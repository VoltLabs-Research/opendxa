#include <volt/pipeline/burgers_loop_builder.h>
#include <volt/pipeline/interface_mesh.h>
#include <volt/helpers/triangle_intersection.h>
#include <volt/core/progress.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/blocked_range.h>
#include <tbb/spin_mutex.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <ranges>
#include <array>
#include <atomic>
#include <algorithm>
#include <utility>

namespace Volt{

namespace{

void sweepFacet(InterfaceMesh::Face* facet, BurgersCircuit* circuit, bool markSwept){
    facet->circuit = circuit;
    if(markSwept) facet->setFlag(InterfaceMeshFace::kSweptBeforeCapFlag);
}

}

BurgersCircuit* BurgersLoopBuilder::allocateCircuit(){
    if(auto* circuit = std::exchange(_unusedCircuit, nullptr)) return circuit;
    return _circuitPool.construct();
}

void BurgersLoopBuilder::traceDislocationSegments(){
    mesh().clearFaceFlag(InterfaceMeshFace::kSweptBeforeCapFlag);

    if(_markCoreAtoms){
        const double alpha = DefectCellIndex::kAlphaScale
            * _mesh.structureAnalysis().maximumNeighborDistance();
        _defectCells.emplace(_mesh.tessellation(), alpha);
        _coreCellClaims = std::vector<CoreCellClaim>(_defectCells->cellCount());
    }

    // One step per circuit length, from 3 up to the stretched maximum.
    ProgressReporter progress("Tracing Burgers circuits",
        static_cast<std::uint64_t>(std::max(0, _maxExtendedBurgersCircuitSize - 2)));
    for(int L : std::views::iota(3, _maxExtendedBurgersCircuitSize + 1)){
        auto dangling = _danglingNodes;

        if(!dangling.empty()){
            for(auto* node : dangling){
                traceSegment(*node, L, L <= _maxBurgersCircuitSize);
            }
        }

        if((L & 1) && L <= _maxBurgersCircuitSize){
            findPrimarySegments(L);
        }
        joinSegments(L);

        if(L >= _maxBurgersCircuitSize && !dangling.empty()){
            tbb::parallel_for_each(dangling.begin(), dangling.end(), [&](DislocationNode* node){
                auto* C = node->circuit;
                if(C->isDangling && C->segmentMeshCap.empty()){
                    C->storeCircuit();
                    C->numPreliminaryPoints = 0;
                }
            });
        }
        progress.update(static_cast<std::uint64_t>(L - 2));
    }
}

void BurgersLoopBuilder::discardCircuit(BurgersCircuit* circuit){
    assert(_unusedCircuit == nullptr);
    _unusedCircuit = circuit;
}

namespace{

// Drops the points each end traced after its circuit stored its cap, so that a dangling line stops at its cap.
void trimPreliminaryPoints(DislocationSegment& segment){
    segment.trimEnds(segment.backwardNode().circuit->numPreliminaryPoints, segment.forwardNode().circuit->numPreliminaryPoints);
}

// Expresses b in a cluster of the reference structure when one transition links the segment's cluster to such a cluster.
void expressInReferenceStructure(DislocationSegment& segment, std::string_view referenceTopologyName){
    auto* cluster = segment.burgersVector.cluster();
    if(cluster->topologyName == referenceTopologyName) return;

    for(auto* t = cluster->transitions; t && t->distance <= 1; t = t->next){
        if(t->cluster2->topologyName == referenceTopologyName){
            segment.burgersVector = ClusterVector(
                t->transform(segment.burgersVector.localVec()),
                t->cluster2
            );
            return;
        }
    }
}

// Runs every line towards +x, +y or +z, along the axis its end-to-end vector is longest on, whichever way it was
// traced. Closed loops, whose ends coincide, keep their direction.
void orientCanonically(DislocationSegment& segment){
    const auto& line = segment.line;
    if(line.empty()) return;
    const Vector3 dir = line.back() - line.front();
    if(dir.isZero(CA_ATOM_VECTOR_EPSILON)) return;

    const auto absx = std::abs(dir.x()), absy = std::abs(dir.y()), absz = std::abs(dir.z());
    if((absx >= absy && absx >= absz && dir.x() < 0) ||
       (absy >= absx && absy >= absz && dir.y() < 0) ||
       (absz >= absx && absz >= absy && dir.z() < 0))
    {
        segment.flipOrientation();
    }
}

}

void BurgersLoopBuilder::finishDislocationSegments(std::string_view referenceTopologyName){
    const auto& segments = network().segments();
    tbb::parallel_for(size_t{0}, segments.size(), [&](size_t i){
        DislocationSegment& segment = *segments[i];
        segment.id = static_cast<int>(i);
        trimPreliminaryPoints(segment);
        expressInReferenceStructure(segment, referenceTopologyName);
        orientCanonically(segment);
    });
}

namespace{

using SearchNode = BurgersLoopSearchNode;

// The search walks only edges that no circuit runs along and whose facet no circuit has swept.
bool isFreeEdge(const InterfaceMesh::Edge* edge){
    return !edge->nextCircuitEdge && !edge->face()->circuit;
}

bool hasFreeEdge(const InterfaceMesh::Vertex* vertex){
    for(auto* edge = vertex->edges(); edge != nullptr; edge = edge->nextVertexEdge()){
        if(isFreeEdge(edge)) return true;
    }
    return false;
}

SearchNode* makeRoot(MemoryPool<SearchNode>& pool, InterfaceMesh::Vertex* vertex){
    auto* root = pool.construct();
    root->vertex = vertex;
    root->coord = Point3::Origin();
    root->tm.setIdentity();
    root->depth = 0;
    root->viaEdge = nullptr;
    root->parent = nullptr;
    return root;
}

SearchNode* makeChild(MemoryPool<SearchNode>& pool, SearchNode* parent, InterfaceMesh::Edge* edge, const Point3& coord){
    auto* child = pool.construct();
    child->vertex = edge->vertex2();
    child->coord = coord;
    child->tm = edge->clusterTransition->isSelfTransition()
        ? parent->tm
        : parent->tm * edge->clusterTransition->reverse->tm;
    child->depth = parent->depth + 1;
    child->viaEdge = edge;
    child->parent = parent;
    return child;
}

// Arriving at a vertex the search already reached closes a circuit. It encloses a dislocation when the two paths to
// the vertex agree on the lattice orientation but not on the lattice coordinate.
bool closesNonZeroCircuit(const SearchNode* from, const InterfaceMesh::Edge* edge, const Point3& coord, const SearchNode* reached){
    const Vector3 b = reached->coord - coord;
    if(b.isZero(CA_LATTICE_VECTOR_EPSILON)) return false;
    const Matrix3 R = from->tm * edge->clusterTransition->reverse->tm;
    return R.equals(reached->tm, CA_TRANSITION_MATRIX_EPSILON);
}

// Breadth-first search over free edges, at most maxDepth edges away from its start vertex.
class CircuitSearch{
public:
    // Calls onClosure(edge, from, reached) for each edge that closes a circuit with a non-zero Burgers vector, and
    // stops when it returns true.
    template<typename OnClosure>
    void run(InterfaceMesh::Vertex* start, int maxDepth, OnClosure&& onClosure){
        _queue.clear();
        _reached.clear();
        _pool.clear(true);
        enqueue(makeRoot(_pool, start));

        for(size_t i = 0; i < _queue.size(); ++i){
            SearchNode* from = _queue[i];
            for(auto* edge = from->vertex->edges(); edge != nullptr; edge = edge->nextVertexEdge()){
                if(!isFreeEdge(edge)) continue;

                const Point3 coord = from->coord + from->tm * edge->clusterVector;
                if(const auto it = _reached.find(edge->vertex2()); it != _reached.end()){
                    if(closesNonZeroCircuit(from, edge, coord, it->second) && onClosure(edge, from, it->second)) return;
                }else if(from->depth < maxDepth){
                    enqueue(makeChild(_pool, from, edge, coord));
                }
            }
        }
    }

private:
    void enqueue(SearchNode* node){
        _queue.push_back(node);
        _reached[node->vertex] = node;
    }

    MemoryPool<SearchNode> _pool;
    std::vector<SearchNode*> _queue;
    std::unordered_map<InterfaceMesh::Vertex*, SearchNode*> _reached;
};

// The first closing edge found from each vertex, sorted so that the serial pass creates circuits in a fixed order
// rather than in the order the threads found them: each circuit it creates changes which edges are free.
std::vector<InterfaceMesh::Edge*> findClosingEdges(const std::vector<InterfaceMesh::Vertex*>& vertices, int searchDepth){
    std::vector<InterfaceMesh::Edge*> closingEdges;
    tbb::spin_mutex closingEdgesMutex;

    tbb::parallel_for(tbb::blocked_range<size_t>(0, vertices.size(), 16), [&](const tbb::blocked_range<size_t>& range){
        CircuitSearch search;
        for(size_t i = range.begin(); i < range.end(); ++i){
            if(!hasFreeEdge(vertices[i])) continue;
            search.run(vertices[i], searchDepth, [&](InterfaceMesh::Edge* edge, SearchNode*, SearchNode*){
                tbb::spin_mutex::scoped_lock lock(closingEdgesMutex);
                closingEdges.push_back(edge);
                return true;
            });
        }
    });

    std::sort(closingEdges.begin(), closingEdges.end(), [](const InterfaceMesh::Edge* a, const InterfaceMesh::Edge* b){
        const int a1 = a->vertex1()->index();
        const int b1 = b->vertex1()->index();
        if(a1 != b1) return a1 < b1;
        return a->vertex2()->index() < b->vertex2()->index();
    });
    return closingEdges;
}

// Frees the edges of a rejected circuit for later searches.
void unlinkCircuit(BurgersCircuit* circuit){
    InterfaceMesh::Edge* edge = circuit->firstEdge;
    do{
        InterfaceMesh::Edge* nextEdge = edge->nextCircuitEdge;
        edge->nextCircuitEdge = nullptr;
        edge->circuit = nullptr;
        edge = nextEdge;
    }while(edge != circuit->firstEdge);
}

// Bounds the walk around a vertex, should its fan of facets not close.
constexpr int kMaxFanSteps = 1000;

// Circuit A enters a vertex over aIn and leaves it over aOut, circuit B over bIn and bOut. Walks the facets around the
// vertex from bIn towards bOut: A crosses B there when it leaves into that side but did not enter from it.
bool leavesOtherCircuit(InterfaceMesh::Edge* aIn, InterfaceMesh::Edge* aOut, InterfaceMesh::Edge* bIn, InterfaceMesh::Edge* bOut){
    const InterfaceMesh::Vertex* vertex = bOut->vertex1();
    assert(aOut->vertex1() == vertex);
    assert(aIn->vertex2() == vertex);
    assert(bIn->vertex2() == vertex);

    bool entersFromSide = false;
    bool leavesIntoSide = false;
    InterfaceMesh::Edge* edge = bIn;
    for(int steps = 0; ; ){
        InterfaceMesh::Edge* nextEdge = edge->nextFaceEdge();
        if(nextEdge == bOut) break;
        InterfaceMesh::Edge* oppositeEdge = nextEdge->oppositeEdge();
        if(oppositeEdge->vertex2() != vertex) break;

        edge = oppositeEdge;
        if(edge == aIn) entersFromSide = true;
        if(nextEdge == aOut) leavesIntoSide = true;
        if(++steps > kMaxFanSteps) break;
    }
    return leavesIntoSide && !entersFromSide;
}

// Whether a circuit through the vertex where circuit edge edge1 ends and edge2 starts crosses that circuit there.
bool crossedAtVertex(InterfaceMesh::Edge* edge1, InterfaceMesh::Edge* edge2){
    InterfaceMesh::Edge* sentinel = edge1->oppositeEdge();
    InterfaceMesh::Edge* edge = sentinel;
    do{
        InterfaceMesh::Edge* incoming = edge->prevFaceEdge();
        const bool circuitPassesThrough = incoming->circuit && incoming->nextCircuitEdge
            && incoming->nextCircuitEdge->vertex1() == incoming->vertex2();
        if(circuitPassesThrough && edge2->oppositeEdge()
           && leavesOtherCircuit(edge2->oppositeEdge(), edge1->oppositeEdge(), incoming, incoming->nextCircuitEdge)){
            return true;
        }
        edge = incoming->oppositeEdge();
    }while(edge != sentinel);
    return false;
}

bool intersectsOtherCircuits(const BurgersCircuit* circuit){
    InterfaceMesh::Edge* edge1 = circuit->firstEdge;
    do{
        InterfaceMesh::Edge* edge2 = edge1->nextCircuitEdge;
        if(edge1 != edge2->oppositeEdge() && crossedAtVertex(edge1, edge2)) return true;
        edge1 = edge2;
    }while(edge1 != circuit->firstEdge);
    return false;
}

}

void BurgersLoopBuilder::findPrimarySegments(int maxBurgersCircuitSize){
    const int searchDepth = (maxBurgersCircuitSize - 1) / 2;
    assert(searchDepth >= 1);

    CircuitSearch search;
    for(InterfaceMesh::Edge* closingEdge : findClosingEdges(mesh().vertices(), searchDepth)){
        if(closingEdge->circuit != nullptr) continue;
        search.run(closingEdge->vertex1(), searchDepth, [&](InterfaceMesh::Edge* edge, SearchNode* from, SearchNode* reached){
            return createPrimaryCircuit(edge, from, reached, maxBurgersCircuitSize);
        });
    }
}

// Builds the circuit that edge closes, out of the search paths that reach its two ends, and starts a segment from it if
// it encloses a dislocation. Returns whether to stop searching from this start vertex: once a segment was started, or
// when the circuit crosses an existing one, whose dislocation is already being traced.
bool BurgersLoopBuilder::createPrimaryCircuit(InterfaceMesh::Edge* edge, SearchNode* from, SearchNode* reached, int maxBurgersCircuitSize){
	if(edge->circuit != nullptr) return false;
	assert(from != reached);

	BurgersCircuit* forwardCircuit = allocateCircuit();
	forwardCircuit->edgeCount = 1;
	forwardCircuit->firstEdge = forwardCircuit->lastEdge = edge->oppositeEdge();
	assert(forwardCircuit->firstEdge->circuit == nullptr);
	forwardCircuit->firstEdge->circuit = forwardCircuit;

	std::unordered_set<InterfaceMesh::Vertex*> fromPath;
	for(SearchNode* a = from; a != nullptr; a = a->parent){
		fromPath.insert(a->vertex);
	}

	SearchNode* meet = reached;
	for(; !fromPath.count(meet->vertex); meet = meet->parent){
		assert(meet->viaEdge != nullptr);
		assert(meet->viaEdge->circuit == nullptr);
		meet->viaEdge->nextCircuitEdge = forwardCircuit->firstEdge;
		forwardCircuit->firstEdge = meet->viaEdge;
		forwardCircuit->edgeCount++;
		forwardCircuit->firstEdge->circuit = forwardCircuit;
	}

	for(SearchNode* a = from; a->vertex != meet->vertex; a = a->parent){
		assert(a->viaEdge != nullptr);
		assert(a->viaEdge->oppositeEdge()->circuit == nullptr);
		forwardCircuit->lastEdge->nextCircuitEdge = a->viaEdge->oppositeEdge();
		forwardCircuit->lastEdge = forwardCircuit->lastEdge->nextCircuitEdge;
		forwardCircuit->edgeCount++;
		forwardCircuit->lastEdge->circuit = forwardCircuit;
	}

	forwardCircuit->lastEdge->nextCircuitEdge = forwardCircuit->firstEdge;
	assert(forwardCircuit->firstEdge != forwardCircuit->firstEdge->nextCircuitEdge);
	assert(forwardCircuit->countEdges() == forwardCircuit->edgeCount);
	assert(forwardCircuit->edgeCount >= 3);

	InterfaceMesh::Edge* e = forwardCircuit->firstEdge;
	Vector3 edgeSum = Vector3::Zero();
	Matrix3 frankRotation = Matrix3::Identity();
	Vector3 b = Vector3::Zero();
	do{
		edgeSum += e->physicalVector;
		b += frankRotation * e->clusterVector;
		if(!e->clusterTransition->isSelfTransition())
			frankRotation = frankRotation * e->clusterTransition->reverse->tm;
		e = e->nextCircuitEdge;
	}while(e != forwardCircuit->firstEdge);
	assert(frankRotation.equals(Matrix3::Identity(), CA_TRANSITION_MATRIX_EPSILON));

	bool intersects = intersectsOtherCircuits(forwardCircuit);
	const bool zeroBurgers = b.isZero(CA_LATTICE_VECTOR_EPSILON);
	const bool badEdgeSum = !edgeSum.isZero(CA_ATOM_VECTOR_EPSILON);
	if(zeroBurgers || badEdgeSum || intersects){
		unlinkCircuit(forwardCircuit);
		discardCircuit(forwardCircuit);
		return intersects;
	}

	createAndTraceSegment(ClusterVector(b, forwardCircuit->firstEdge->clusterTransition->cluster1), forwardCircuit, maxBurgersCircuitSize);
	return true;
}

void BurgersLoopBuilder::createAndTraceSegment(const ClusterVector& burgersVector, BurgersCircuit* forwardCircuit, int maxCircuitLength){
	BurgersCircuit* backwardCircuit = buildReverseCircuit(forwardCircuit);

	DislocationSegment* segment = network().createSegment(burgersVector);
	segment->forwardNode().circuit = forwardCircuit;
	segment->backwardNode().circuit = backwardCircuit;
	forwardCircuit->dislocationNode = &segment->forwardNode();
	backwardCircuit->dislocationNode = &segment->backwardNode();

	_danglingNodes.push_back(&segment->forwardNode());
	_danglingNodes.push_back(&segment->backwardNode());

	segment->addPoint(true, backwardCircuit->calculateCenter(), static_cast<int>(backwardCircuit->countEdges()));

	appendLinePoint(segment->forwardNode());

	traceSegment(segment->forwardNode(), maxCircuitLength, true);
	traceSegment(segment->backwardNode(), maxCircuitLength, true);
}

BurgersCircuit* BurgersLoopBuilder::buildReverseCircuit(BurgersCircuit* forwardCircuit){
	BurgersCircuit* backwardCircuit = allocateCircuit();

	backwardCircuit->edgeCount = 0;
	backwardCircuit->firstEdge = nullptr;
	backwardCircuit->lastEdge = nullptr;
	InterfaceMesh::Edge* edge1 = forwardCircuit->firstEdge;
	do{
		InterfaceMesh::Edge* edge2 = edge1->nextCircuitEdge;
		InterfaceMesh::Edge* oppositeEdge1 = edge1->oppositeEdge();
		InterfaceMesh::Edge* oppositeEdge2 = edge2->oppositeEdge();
		InterfaceMesh::Face* facet1 = oppositeEdge1->face();
		InterfaceMesh::Face* facet2 = oppositeEdge2->face();
		assert(facet1 != nullptr && facet2 != nullptr);
		assert(facet1->circuit == nullptr || facet1->circuit == backwardCircuit);
		assert(facet2->circuit == nullptr || facet2->circuit == backwardCircuit);
		assert(edge1->vertex2() == edge2->vertex1());
		assert((edge1->clusterVector + oppositeEdge1->clusterTransition->tm * oppositeEdge1->clusterVector).isZero(CA_LATTICE_VECTOR_EPSILON));
		assert((edge2->clusterVector + oppositeEdge2->clusterTransition->tm * oppositeEdge2->clusterVector).isZero(CA_LATTICE_VECTOR_EPSILON));

		if(facet1 != facet2){
			InterfaceMesh::Edge* innerEdge1 = oppositeEdge1->prevFaceEdge()->oppositeEdge();
			InterfaceMesh::Edge* innerEdge2 = oppositeEdge2->nextFaceEdge()->oppositeEdge();
			assert(innerEdge1 != nullptr && innerEdge2 != nullptr);
			assert(innerEdge1->vertex1() == edge1->vertex2());
			assert(innerEdge2->vertex2() == edge1->vertex2());
			assert(innerEdge1->vertex1() == innerEdge2->vertex2());
			assert(innerEdge1->circuit == nullptr || innerEdge1->circuit == backwardCircuit);
			assert(innerEdge2->circuit == nullptr || innerEdge2->circuit == backwardCircuit);
			sweepFacet(facet1, backwardCircuit, true);
			sweepFacet(facet2, backwardCircuit, true);
			innerEdge1->circuit = backwardCircuit;
			innerEdge2->circuit = backwardCircuit;
			innerEdge2->nextCircuitEdge = innerEdge1;

			if(backwardCircuit->lastEdge == nullptr){
				assert(backwardCircuit->firstEdge == nullptr);
				assert(innerEdge1->nextCircuitEdge == nullptr);
				backwardCircuit->lastEdge = innerEdge1;
				backwardCircuit->firstEdge = innerEdge2;
				backwardCircuit->edgeCount += 2;
			}else if(backwardCircuit->lastEdge != innerEdge2){
				if(innerEdge1 != backwardCircuit->firstEdge){
					innerEdge1->nextCircuitEdge = backwardCircuit->firstEdge;
					backwardCircuit->edgeCount += 2;
				}else{
					backwardCircuit->edgeCount += 1;
				}
				backwardCircuit->firstEdge = innerEdge2;
			}else if(backwardCircuit->firstEdge != innerEdge1){
				innerEdge1->nextCircuitEdge = backwardCircuit->firstEdge;
				backwardCircuit->firstEdge = innerEdge1;
				backwardCircuit->edgeCount += 1;
			}

			assert(innerEdge1->vertex1() != innerEdge1->vertex2());
			assert(innerEdge2->vertex1() != innerEdge2->vertex2());
		}

		edge1 = edge2;
	}while(edge1 != forwardCircuit->firstEdge);

	assert(backwardCircuit->lastEdge->vertex2() == backwardCircuit->firstEdge->vertex1());
	assert(backwardCircuit->lastEdge->nextCircuitEdge == nullptr || backwardCircuit->lastEdge->nextCircuitEdge == backwardCircuit->firstEdge);

	backwardCircuit->lastEdge->nextCircuitEdge = backwardCircuit->firstEdge;

	assert(backwardCircuit->firstEdge != backwardCircuit->firstEdge->nextCircuitEdge);
	assert(backwardCircuit->countEdges() == backwardCircuit->edgeCount);
	assert(backwardCircuit->edgeCount >= 3);
	assert(!backwardCircuit->calculateBurgersVector().localVec().isZero(CA_LATTICE_VECTOR_EPSILON));

	return backwardCircuit;
}

// Sweeps the circuit forward until it stops moving, then lengthens it by one edge and sweeps again, until it is
// maxCircuitLength edges long or cannot be lengthened. Every move appends a point to the node's line.
void BurgersLoopBuilder::traceSegment(DislocationNode& node, int maxCircuitLength, bool markSwept){
    const BurgersCircuit& circuit = *node.circuit;
    assert(circuit.countEdges() == circuit.edgeCount);
    assert(circuit.isDangling);

    do{
        sweepCircuit(node, markSwept);
    }while(circuit.edgeCount < maxCircuitLength && lengthenCircuit(node, markSwept));
}

// Each lap starts at a different edge, cycling through the edges of every circuit traced so far.
InterfaceMesh::Edge* BurgersLoopBuilder::nextStartEdge(const BurgersCircuit& circuit){
    return circuit.getEdge(_edgeStartIndex++ % circuit.edgeCount);
}

// Laps the circuit, moving it wherever one of the sweep moves applies, until a full lap makes no move.
void BurgersLoopBuilder::sweepCircuit(DislocationNode& node, bool markSwept){
    BurgersCircuit& circuit = *node.circuit;
    InterfaceMesh::Edge* edge0 = nextStartEdge(circuit);
    InterfaceMesh::Edge* edge1 = edge0->nextCircuitEdge;
    InterfaceMesh::Edge* edge2 = edge1->nextCircuitEdge;
    assert(edge1->circuit == &circuit);
    int counter = 0;
    do{
        assert(circuit.edgeCount >= 3);
        // The moves keep the Burgers vector, and both circuit creators reject a zero one.
        assert(!circuit.calculateBurgersVector().localVec().isZero(CA_LATTICE_VECTOR_EPSILON));
        assert(circuit.countEdges() == circuit.edgeCount);
        assert(edge0->circuit == &circuit && edge1->circuit == &circuit && edge2->circuit == &circuit);

        const bool moved = tryRemoveTwoCircuitEdges(edge0, edge1, edge2)
            || tryRemoveThreeCircuitEdges(edge0, edge1, edge2, markSwept)
            || tryRemoveOneCircuitEdge(edge0, edge1, edge2, markSwept)
            || trySweepTwoFacets(edge0, edge1, edge2, markSwept);
        if(moved){
            appendLinePoint(node);
            counter = -1;
        }

        edge0 = edge1;
        edge1 = edge2;
        edge2 = edge2->nextCircuitEdge;
        counter++;
    }while(counter <= circuit.edgeCount);
    assert(circuit.edgeCount >= 3);
    assert(circuit.countEdges() == circuit.edgeCount);
}

// Lengthens the circuit by one edge over the first facet ahead of it that allows that.
bool BurgersLoopBuilder::lengthenCircuit(DislocationNode& node, bool markSwept){
    BurgersCircuit& circuit = *node.circuit;
    InterfaceMesh::Edge* firstEdge = nextStartEdge(circuit);
    InterfaceMesh::Edge* edge0 = firstEdge;
    InterfaceMesh::Edge* edge1 = firstEdge->nextCircuitEdge;
    do{
        if(tryInsertOneCircuitEdge(edge0, edge1, markSwept)){
            appendLinePoint(node);
            return true;
        }

        edge0 = edge1;
        edge1 = edge1->nextCircuitEdge;
    }while(edge0 != firstEdge);
    return false;
}

bool BurgersLoopBuilder::tryRemoveTwoCircuitEdges(InterfaceMesh::Edge*& edge0, InterfaceMesh::Edge*& edge1, InterfaceMesh::Edge*& edge2){
	if(edge1 != edge2->oppositeEdge()) return false;

	BurgersCircuit* circuit = edge0->circuit;
	assert(circuit->edgeCount >= 4);
	edge0->nextCircuitEdge = edge2->nextCircuitEdge;
	
	if(edge0 == circuit->lastEdge){
		circuit->firstEdge = circuit->lastEdge->nextCircuitEdge;
	}else if(edge1 == circuit->lastEdge){
		circuit->lastEdge = edge0;
		circuit->firstEdge = edge0->nextCircuitEdge;
	}else if(edge2 == circuit->lastEdge){
		circuit->lastEdge = edge0;
	}

	circuit->edgeCount -= 2;
	assert(circuit->edgeCount >= 0);

	edge1 = edge0->nextCircuitEdge;
	edge2 = edge1->nextCircuitEdge;
	return true;
}

bool BurgersLoopBuilder::tryRemoveThreeCircuitEdges(
	InterfaceMesh::Edge*& edge0, 
	InterfaceMesh::Edge*& edge1, 
	InterfaceMesh::Edge*& edge2,
	bool markSwept
){
	InterfaceMesh::Face* facet1 = edge1->face();
	InterfaceMesh::Face* facet2 = edge2->face();

	if(facet2 != facet1 || facet1->circuit != nullptr) return false;

	BurgersCircuit* circuit = edge0->circuit;
	assert(circuit->edgeCount > 2);
	InterfaceMesh::Edge* edge3 = edge2->nextCircuitEdge;

	if(edge3->face() != facet1) return false;
	assert(circuit->edgeCount > 4);

	edge0->nextCircuitEdge = edge3->nextCircuitEdge;

	if(edge2 == circuit->firstEdge || edge3 == circuit->firstEdge){
		circuit->firstEdge = edge3->nextCircuitEdge;
		circuit->lastEdge = edge0;
	}else if(edge1 == circuit->firstEdge){
		circuit->firstEdge = edge3->nextCircuitEdge;
		assert(circuit->lastEdge == edge0);
	}else if(edge3 == circuit->lastEdge){
		circuit->lastEdge = edge0;
	}

	circuit->edgeCount -= 3;
	edge1 = edge3->nextCircuitEdge;
	edge2 = edge1->nextCircuitEdge;

	sweepFacet(facet1, circuit, markSwept);

	return true;
}

bool BurgersLoopBuilder::tryRemoveOneCircuitEdge(
	InterfaceMesh::Edge*& edge0, 
	InterfaceMesh::Edge*& edge1, 
	InterfaceMesh::Edge*& edge2, 
	bool markSwept
){
	InterfaceMesh::Face* facet1 = edge1->face();
	InterfaceMesh::Face* facet2 = edge2->face();
	if(facet2 != facet1 || facet1->circuit != nullptr) return false;

	BurgersCircuit* circuit = edge0->circuit;
	assert(circuit->edgeCount > 2);

	if(edge0->face() == facet1) return false;

	InterfaceMesh::Edge* shortEdge = edge1->prevFaceEdge()->oppositeEdge();
	assert(shortEdge->vertex1() == edge1->vertex1());
	assert(shortEdge->vertex2() == edge2->vertex2());

	if(shortEdge->circuit != nullptr) return false;

	assert(shortEdge->nextCircuitEdge == nullptr);
	shortEdge->nextCircuitEdge = edge2->nextCircuitEdge;
	assert(shortEdge != edge2->nextCircuitEdge->oppositeEdge());
	assert(shortEdge != edge0->oppositeEdge());
	edge0->nextCircuitEdge = shortEdge;
	if(edge0 == circuit->lastEdge){
		assert(circuit->lastEdge != edge2);
		assert(circuit->firstEdge == edge1);
		assert(shortEdge != circuit->lastEdge->oppositeEdge());
		circuit->firstEdge = shortEdge;
	}

	if(edge2 == circuit->lastEdge){
		circuit->lastEdge = shortEdge;
	}else if(edge2 == circuit->firstEdge){
		circuit->firstEdge = shortEdge->nextCircuitEdge;
		circuit->lastEdge = shortEdge;
	}

	circuit->edgeCount -= 1;
	edge1 = shortEdge;
	edge2 = shortEdge->nextCircuitEdge;
	shortEdge->circuit = circuit;

	sweepFacet(facet1, circuit, markSwept);

	return true;
}

bool BurgersLoopBuilder::trySweepTwoFacets(
	InterfaceMesh::Edge*& edge0, 
	InterfaceMesh::Edge*& edge1, 
	InterfaceMesh::Edge*& edge2, 
	bool markSwept
){
	InterfaceMesh::Face* facet1 = edge1->face();
	InterfaceMesh::Face* facet2 = edge2->face();

	if(facet1->circuit != nullptr || facet2->circuit != nullptr) return false;

	BurgersCircuit* circuit = edge0->circuit;
	if(facet1 == facet2 || circuit->edgeCount <= 2) return false;

	InterfaceMesh::Edge* outerEdge1 = edge1->prevFaceEdge()->oppositeEdge();
	InterfaceMesh::Edge* innerEdge1 = edge1->nextFaceEdge();
	InterfaceMesh::Edge* outerEdge2 = edge2->nextFaceEdge()->oppositeEdge();
	InterfaceMesh::Edge* innerEdge2 = edge2->prevFaceEdge();

	if(innerEdge1 != innerEdge2->oppositeEdge() || outerEdge1->circuit != nullptr || outerEdge2->circuit != nullptr)
		return false;

	assert(outerEdge1->nextCircuitEdge == nullptr);
	assert(outerEdge2->nextCircuitEdge == nullptr);
	outerEdge1->nextCircuitEdge = outerEdge2;
	outerEdge2->nextCircuitEdge = edge2->nextCircuitEdge;
	edge0->nextCircuitEdge = outerEdge1;

	if(edge0 == circuit->lastEdge){
		circuit->firstEdge = outerEdge1;
	}else if(edge1 == circuit->lastEdge){
		circuit->lastEdge = outerEdge1;
		circuit->firstEdge = outerEdge2;
	}else if(edge2 == circuit->lastEdge){
		circuit->lastEdge = outerEdge2;
	}

	outerEdge1->circuit = circuit;
	outerEdge2->circuit = circuit;

	sweepFacet(facet1, circuit, markSwept);
	sweepFacet(facet2, circuit, markSwept);

	edge0 = outerEdge1;
	edge1 = outerEdge2;
	edge2 = edge1->nextCircuitEdge;

	return true;
}

bool BurgersLoopBuilder::tryInsertOneCircuitEdge(
	InterfaceMesh::Edge*& edge0, 
	InterfaceMesh::Edge*& edge1, 
	bool markSwept
){
	assert(edge0 != edge1->oppositeEdge());

	InterfaceMesh::Face* facet = edge1->face();
	if(facet->circuit != nullptr) return false;

	InterfaceMesh::Edge* insertEdge1 = edge1->prevFaceEdge()->oppositeEdge();
	if(insertEdge1->circuit != nullptr) return false;

	InterfaceMesh::Edge* insertEdge2 = edge1->nextFaceEdge()->oppositeEdge();
	if(insertEdge2->circuit != nullptr) return false;

	assert(insertEdge1->nextCircuitEdge == nullptr);
	assert(insertEdge2->nextCircuitEdge == nullptr);

	BurgersCircuit* circuit = edge0->circuit;
	
	insertEdge1->nextCircuitEdge = insertEdge2;
	insertEdge2->nextCircuitEdge = edge1->nextCircuitEdge;
	
	edge0->nextCircuitEdge = insertEdge1;
	if(edge0 == circuit->lastEdge){
		circuit->firstEdge = insertEdge1;
	}else if(edge1 == circuit->lastEdge){
		circuit->lastEdge = insertEdge2;
	}

	insertEdge1->circuit = circuit;
	insertEdge2->circuit = circuit;
	circuit->edgeCount++;

	assert(circuit->countEdges() == circuit->edgeCount);

	sweepFacet(facet, circuit, markSwept);

	return true;
}

void BurgersLoopBuilder::appendLinePoint(DislocationNode& node){
	DislocationSegment& segment = *node.segment;
	assert(!segment.line.empty());

	const Point3& lastPoint = node.position();
	const Point3 newPoint = lastPoint + cell().wrapVector(node.circuit->calculateCenter() - lastPoint);
	segment.addPoint(node.isForwardNode(), newPoint, static_cast<int>(node.circuit->edgeCount));

	node.circuit->numPreliminaryPoints++;

	if(_markCoreAtoms){
		markCoreCells(node, newPoint);
	}
}

void BurgersLoopBuilder::markCoreCells(DislocationNode& node, const Point3& center){
	if(node.circuit->firstEdge == nullptr) return;

	const DelaunayTessellation& tessellation = _mesh.tessellation();
	const bool capStored = !node.circuit->segmentMeshCap.empty();

	std::vector<TriangleIntersection::Triangle3> cap;
	cap.reserve(node.circuit->edgeCount);

	Box3 capBounds;
	auto* edge = node.circuit->firstEdge;
	do{
		const TriangleIntersection::Triangle3 triangle{
			center + cell().wrapVector(edge->vertex1()->pos() - center),
			center + cell().wrapVector(edge->vertex2()->pos() - center),
			center
		};
		for(const Point3& corner : triangle){
			capBounds.addPoint(corner);
		}
		cap.push_back(triangle);
		edge = edge->nextCircuitEdge;
	}while(edge != nullptr && edge != node.circuit->firstEdge);

	std::vector<DelaunayTessellation::CellHandle> candidates;
	_defectCells->queryOverlapping(capBounds, candidates);

	for(auto candidate : candidates){
		const int slot = tessellation.getUserField(candidate);
		if(slot < 0 || static_cast<std::size_t>(slot) >= _coreCellClaims.size()) continue;

		CoreCellClaim& claim = _coreCellClaims[slot];
		if(claim.node != nullptr) continue;

		std::array<Point3, 4> tetrahedron;
		for(int corner = 0; corner < 4; ++corner){
			tetrahedron[corner] = tessellation.vertexPosition(tessellation.cellVertex(candidate, corner));
		}

		const bool intersects = std::ranges::any_of(cap, [&](const auto& triangle){
			for(int face = 0; face < 4; ++face){
				const TriangleIntersection::Triangle3 facet{
					tetrahedron[DelaunayTessellation::cellFacetVertexIndex(face, 0)],
					tetrahedron[DelaunayTessellation::cellFacetVertexIndex(face, 1)],
					tetrahedron[DelaunayTessellation::cellFacetVertexIndex(face, 2)]
				};
				if(TriangleIntersection::overlap(facet, triangle)) return true;
			}
			return false;
		});

		if(!intersects) continue;

		claim.node = &node;
		claim.claimedWhileCapStored = capStored;
	}
}

std::vector<int> BurgersLoopBuilder::coreAtomDislocationIds(std::size_t particleCount) const{
	std::vector<int> ids(particleCount, -1);
	if(_coreCellClaims.empty()) return ids;

	const DelaunayTessellation& tessellation = _mesh.tessellation();

	for(auto cellHandle : tessellation.cells()){
		const int slot = tessellation.getUserField(cellHandle);
		if(slot < 0 || static_cast<std::size_t>(slot) >= _coreCellClaims.size()) continue;

		const CoreCellClaim& claim = _coreCellClaims[slot];
		DislocationNode* node = claim.node;
		if(node == nullptr) continue;
		if(node->isDangling() && claim.claimedWhileCapStored) continue;

		const DislocationSegment* segment = node->segment;
		while(segment->replacedWith != nullptr){
			segment = segment->replacedWith;
		}

		for(int corner = 0; corner < 4; ++corner){
			const int particle = tessellation.vertexIndex(tessellation.cellVertex(cellHandle, corner));
			if(particle >= 0 && static_cast<std::size_t>(particle) < ids.size()){
				ids[particle] = segment->id;
			}
		}
	}

	return ids;
}

void BurgersLoopBuilder::joinSegments(int maxCircuitLength){
	createSecondarySegments(maxCircuitLength);
	linkAdjacentCircuits();
	for(DislocationNode* node : danglingNodes()){
		joinArms(node);
	}

	std::erase_if(_danglingNodes, [](DislocationNode* node){
		return !node->isDangling();
	});
}

// Tries a secondary segment wherever a dangling circuit borders edges no circuit holds, as where a line branches. Also
// visits the nodes of the segments it creates.
void BurgersLoopBuilder::createSecondarySegments(int maxCircuitLength){
	for(size_t nodeIndex = 0; nodeIndex < danglingNodes().size(); nodeIndex++){
		BurgersCircuit* circuit = danglingNodes()[nodeIndex]->circuit;
		assert(circuit->isDangling);

		InterfaceMesh::Edge* edge = circuit->firstEdge;
		do{
			assert(edge->circuit == circuit);
			BurgersCircuit* oppositeCircuit = edge->oppositeEdge()->circuit;
			if(oppositeCircuit == nullptr){
				assert(edge->oppositeEdge()->nextCircuitEdge == nullptr);

				createSecondarySegment(edge, circuit, maxCircuitLength);

				while(edge->oppositeEdge()->circuit == nullptr && edge != circuit->firstEdge){
					edge = edge->nextCircuitEdge;
				}
			}else{
				edge = edge->nextCircuitEdge;
			}
		}
		while(edge != circuit->firstEdge);
	}
}

// Puts the nodes whose circuits share edges into one junction ring, and marks a circuit completely blocked when other
// circuits hold all the edges beyond it.
void BurgersLoopBuilder::linkAdjacentCircuits(){
	for(DislocationNode* node : danglingNodes()){
		BurgersCircuit* circuit = node->circuit;
		assert(circuit->isDangling);

		circuit->isCompletelyBlocked = true;
		InterfaceMesh::Edge* edge = circuit->firstEdge;
		do{
			assert(edge->circuit == circuit);
			BurgersCircuit* adjacentCircuit = edge->oppositeEdge()->circuit;
			if(adjacentCircuit == nullptr){
				circuit->isCompletelyBlocked = false;
				break;
			}
			if(adjacentCircuit != circuit){
				assert(adjacentCircuit->isDangling);
				DislocationNode* adjacentNode = adjacentCircuit->dislocationNode;
				if(!node->formsJunctionWith(adjacentNode)){
					node->connectNodes(adjacentNode);
				}
			}
			edge = edge->nextCircuitEdge;
		}while(edge != circuit->firstEdge);
	}
}

// Joins the arms of node's junction ring once all their circuits are completely blocked: three or more arms meet at a
// junction, and two arms are one line, which closes a loop or merges two segments. A ring with an arm that can still
// move is dissolved; the next length links it again.
void BurgersLoopBuilder::joinArms(DislocationNode* node){
	if(!node->circuit->isDangling) return;
	if(!node->circuit->isCompletelyBlocked){
		node->dissolveJunction();
		return;
	}
	if(node->junctionRing == node) return;
	assert(node->segment->replacedWith == nullptr);

	Vector3 centerOfMassVector = Vector3::Zero();
	const Point3 basePoint = node->position();
	int armCount = 1;
	for(DislocationNode* armNode = node->junctionRing; armNode != node; armNode = armNode->junctionRing){
		assert(armNode->segment->replacedWith == nullptr);
		assert(armNode->circuit->isDangling);
		if(!armNode->circuit->isCompletelyBlocked){
			node->dissolveJunction();
			return;
		}
		armCount++;
		centerOfMassVector += cell().wrapVector(armNode->position() - basePoint);
	}
	assert(armCount >= 2);

	if(armCount >= 3){
		formJunction(node, basePoint + centerOfMassVector / armCount);
		return;
	}

	DislocationNode* node2 = node->junctionRing;
	assert(node != node2);
	assert(node2->junctionRing == node);
	for(DislocationNode* end : {node, node2}){
		end->circuit->isDangling = false;
		end->circuit->numPreliminaryPoints = 0;
	}

	if(node->oppositeNode == node2){
		closeLoop(node, node2);
	}else{
		mergeSegments(node, node2);
	}
}

// Ends every arm's line at the arms' common center of mass.
void BurgersLoopBuilder::formJunction(DislocationNode* node, const Point3& centerOfMass){
	DislocationNode* armNode = node;
	do{
		armNode->circuit->isDangling = false;
		assert(armNode != armNode->junctionRing);

		DislocationSegment& segment = *armNode->segment;
		const bool atForwardEnd = armNode->isForwardNode();
		const Point3& end = armNode->position();
		const int endCoreSize = atForwardEnd ? segment.coreSize.back() : segment.coreSize.front();
		segment.addPoint(atForwardEnd, end + cell().wrapVector(centerOfMass - end), endCoreSize);

		armNode->circuit->numPreliminaryPoints = 0;
		armNode = armNode->junctionRing;
	}while(armNode != node);
}

// The two ends of one segment met. Unless they already coincide, the line gets a last point at the periodic image of
// its first one nearest its end, which closes it.
void BurgersLoopBuilder::closeLoop(DislocationNode* node1, DislocationNode* node2){
	assert(node1->segment == node2->segment);
	DislocationSegment* loop = node1->segment;
	assert(loop->isClosedLoop());

	if(!cell().wrapVector(node1->position() - node2->position()).isZero(CA_ATOM_VECTOR_EPSILON)){
		loop->addPoint(true, loop->line.back() + cell().wrapVector(loop->line.front() - loop->line.back()), loop->coreSize.back());
		assert(cell().wrapVector(node1->position() - node2->position()).isZero(CA_ATOM_VECTOR_EPSILON));
	}

	assert(loop->line.size() >= 3);
}

// Ends node1 and node2 of two segments met: moves segment2's line onto segment1's end at node1, shifted by the
// periodic image vector between the two ends, and discards segment2.
void BurgersLoopBuilder::mergeSegments(DislocationNode* node1, DislocationNode* node2){
	assert(node1->segment != node2->segment);

	DislocationNode* farEnd1 = node1->oppositeNode;
	DislocationNode* farEnd2 = node2->oppositeNode;
	DislocationSegment* segment1 = node1->segment;
	DislocationSegment* segment2 = node2->segment;

	// segment2's points and core sizes, from the one after node2's to farEnd2's
	std::vector<Point3> points;
	std::vector<int> coreSizes;
	if(node2->isBackwardNode()){
		points.assign(segment2->line.begin() + 1, segment2->line.end());
		coreSizes.assign(segment2->coreSize.begin() + 1, segment2->coreSize.end());
	}else{
		points.assign(segment2->line.rbegin() + 1, segment2->line.rend());
		coreSizes.assign(segment2->coreSize.rbegin() + 1, segment2->coreSize.rend());
	}

	const Vector3 shiftVector = calculateShiftVector(node1->position(), node2->position());
	if(shiftVector != Vector3::Zero()){
		for(Point3& point : points){
			point -= shiftVector;
		}
	}

	if(node1->isForwardNode()){
		segment1->line.insert(segment1->line.end(), points.begin(), points.end());
		segment1->coreSize.insert(segment1->coreSize.end(), coreSizes.begin(), coreSizes.end());
		segment1->nodes[0] = farEnd2;
	}else{
		segment1->line.insert(segment1->line.begin(), points.rbegin(), points.rend());
		segment1->coreSize.insert(segment1->coreSize.begin(), coreSizes.rbegin(), coreSizes.rend());
		segment1->nodes[1] = farEnd2;
	}

	farEnd2->segment = segment1;
	farEnd2->oppositeNode = farEnd1;
	farEnd1->oppositeNode = farEnd2;
	node1->oppositeNode = node2;
	node2->oppositeNode = node1;
	segment2->replacedWith = segment1;
	network().discardSegment(segment2);
}

void BurgersLoopBuilder::createSecondarySegment(InterfaceMesh::Edge* firstEdge, BurgersCircuit* outerCircuit, int maxCircuitLength){
	assert(firstEdge->circuit == outerCircuit);

	int edgeCount = 1;
	Vector3 burgersVector = Vector3::Zero();
	Vector3 edgeSum = Vector3::Zero();
	Cluster* baseCluster = nullptr;
	Matrix3 frankRotation = Matrix3::Identity();
	int numCircuits = 1;
	InterfaceMesh::Edge* circuitStart = firstEdge->oppositeEdge();
	InterfaceMesh::Edge* circuitEnd = circuitStart;
	InterfaceMesh::Edge* edge = circuitStart;
	for(;;){
		for(;;){
			assert(edge->circuit == nullptr);
			InterfaceMesh::Edge* oppositeEdge = edge->oppositeEdge();
			InterfaceMesh::Edge* nextEdge = oppositeEdge->prevFaceEdge();
			assert(nextEdge->vertex2() == oppositeEdge->vertex1());
			assert(nextEdge->vertex2() == edge->vertex2());
			if(nextEdge->circuit != nullptr){
				if(nextEdge->circuit != outerCircuit){
					outerCircuit = nextEdge->circuit;
					numCircuits++;
				}
				
				edge = nextEdge->oppositeEdge();
				break;
			}
			edge = nextEdge;
		}

		circuitEnd->nextCircuitEdge = edge;
		edgeSum += edge->physicalVector;
		burgersVector += frankRotation * edge->clusterVector;
		if(!baseCluster) baseCluster = edge->clusterTransition->cluster1;
		if(!edge->clusterTransition->isSelfTransition()){
			frankRotation = frankRotation * edge->clusterTransition->reverse->tm;
		}

		if(edge == circuitStart) break;
		circuitEnd = edge;
		edgeCount++;

		if(edgeCount > maxCircuitLength) break;
	}

	const bool rejectNumCircuits = numCircuits == 1;
	const bool rejectTooLong = edgeCount > maxCircuitLength;
	const bool rejectZeroBurgers = burgersVector.isZero(CA_LATTICE_VECTOR_EPSILON);
	const bool rejectEdgeSum = edgeSum.isZero(CA_ATOM_VECTOR_EPSILON) == false;
	const bool rejectFrank = !frankRotation.equals(Matrix3::Identity(), CA_TRANSITION_MATRIX_EPSILON);
	if(rejectNumCircuits || rejectTooLong || rejectZeroBurgers || rejectEdgeSum || rejectFrank){

		edge = circuitStart;
		for(;;){
			assert(edge->circuit == nullptr);
			InterfaceMesh::Edge* nextEdge = edge->nextCircuitEdge;
			edge->nextCircuitEdge = nullptr;
			if(edge == circuitEnd) break;
			edge = nextEdge;
		}
		return;
	}
	assert(circuitStart != circuitEnd);

	BurgersCircuit* forwardCircuit = allocateCircuit();
	forwardCircuit->firstEdge = circuitStart;
	forwardCircuit->lastEdge = circuitEnd;
	forwardCircuit->edgeCount = edgeCount;
	edge = circuitStart;
	do{
		assert(edge->circuit == nullptr);
		edge->circuit = forwardCircuit;
		edge = edge->nextCircuitEdge;
	}while(edge != circuitStart);
	
	assert(forwardCircuit->countEdges() == forwardCircuit->edgeCount);

	createAndTraceSegment(ClusterVector(burgersVector, baseCluster), forwardCircuit, maxCircuitLength);
}

}
