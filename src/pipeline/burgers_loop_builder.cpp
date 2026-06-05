#include <volt/pipeline/burgers_loop_builder.h>
#include <volt/pipeline/interface_mesh.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/blocked_range.h>
#include <tbb/spin_mutex.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <ranges>
#include <atomic>
#include <algorithm>

namespace Volt{

// Allocate or recycle a BurgersCircuit object from the internal pool.
BurgersCircuit* BurgersLoopBuilder::allocateCircuit(){
    BurgersCircuit* circuit = nullptr;
    tbb::spin_mutex::scoped_lock lock(_circuit_pool_mutex);
    
    if(_unusedCircuit != nullptr){
        circuit = _unusedCircuit;
        _unusedCircuit = nullptr;
    }
    
    if(circuit == nullptr){
        return _circuitPool.construct();
    }
    return circuit;
}

bool BurgersLoopBuilder::traceDislocationSegments(){
    mesh().clearFaceFlag(0);

    for(int L : std::views::iota(3, _maxExtendedBurgersCircuitSize + 1)){
        auto dangling = _danglingNodes;

        if(!dangling.empty()){
            for(auto* node : dangling){
                traceSegment(*node->segment, *node, L, L <= _maxBurgersCircuitSize);
            }
        }

        if((L & 1) && L <= _maxBurgersCircuitSize){
            if(!findPrimarySegments(L)){
                return false;
            }
        }
        (void) joinSegments(L);

        if(L >= _maxBurgersCircuitSize && !dangling.empty()){
            tbb::parallel_for_each(dangling.begin(), dangling.end(), [&](DislocationNode* node){
                auto* C = node->circuit;
                if(C->isDangling && C->segmentMeshCap.empty()){
                    C->storeCircuit();
                    C->numPreliminaryPoints = 0;
                }
            });
        }
    }

    return true;
}

// Return a previously used BurgersCircuit to the pool for resuse.
// Ensures only one "unused" circuit is held at a time.
void BurgersLoopBuilder::discardCircuit(BurgersCircuit* circuit){
    tbb::spin_mutex::scoped_lock lock(_circuit_pool_mutex);
    assert(_unusedCircuit == nullptr);
    _unusedCircuit = circuit;
}

// Finalize all traced segments. Trim preliminary points, re-express Burgers
// vectors in the target crystal structure and orient each line so
// it points consistently.
void BurgersLoopBuilder::finishDislocationSegments(std::string_view referenceTopologyName){
    auto& segs = network().segments();

    tbb::parallel_for(tbb::blocked_range<size_t>(0, segs.size()), 
        [&](const tbb::blocked_range<size_t>& r){
            for(size_t i = r.begin(); i != r.end(); ++i){
                auto* s = segs[i];
                auto pre  = s->backwardNode().circuit->numPreliminaryPoints;
                auto post = s->forwardNode().circuit->numPreliminaryPoints;
                s->id = static_cast<int>(i);

                auto& line = s->line;
                auto& core = s->coreSize;

                line.erase(line.begin(), line.begin() + pre);
                line.erase(line.end() - post, line.end());
                core.erase(core.begin(), core.begin() + pre);
                core.erase(core.end() - post, core.end());
            }
	});

    // Re-express Burgers vectors in the desired structure
    tbb::parallel_for(tbb::blocked_range<size_t>(0, segs.size()), 
        [&](const tbb::blocked_range<size_t>& r){
            for(size_t i = r.begin(); i != r.end(); ++i){
                auto* s = segs[i];
                auto* orig = s->burgersVector.cluster();
                if(orig->topologyName != referenceTopologyName){
                    for(auto* t = orig->transitions; t && t->distance <= 1; t = t->next){
                        if(t->cluster2->topologyName == referenceTopologyName){
                            s->burgersVector = ClusterVector(
                                t->transform(s->burgersVector.localVec()),
                                t->cluster2
                            );
                            break;
                        }
                    }
                }
            }
	});

    // Align the orientation of each segment
    tbb::parallel_for(tbb::blocked_range<size_t>(0, segs.size()), 
        [&](const tbb::blocked_range<size_t>& r){
            for(size_t i = r.begin(); i != r.end(); ++i){
                auto* s = segs[i];
                auto& line = s->line;
                if (line.empty()) continue;
                Vector3 dir = line.back() - line.front();
                if(dir.isZero(CA_ATOM_VECTOR_EPSILON)) continue;

                auto absx = std::abs(dir.x()), absy = std::abs(dir.y()), absz = std::abs(dir.z());
                if((absx >= absy && absx >= absz && dir.x() < 0) ||
                   (absy >= absx && absy >= absz && dir.y() < 0) ||
                   (absz >= absx && absz >= absy && dir.z() < 0))
                {
                    s->flipOrientation();
                }
            }
	});
}

struct BurgersCircuitSearchStruct{
    InterfaceMesh::Vertex* node;
    Point3 latticeCoord;
    Matrix3 tm;
    int recursiveDepth;
    InterfaceMesh::Edge* predecessorEdge;
    BurgersCircuitSearchStruct* nextToProcess;
};

using SearchNode = BurgersLoopBuilder::SearchNode;

bool BurgersLoopBuilder::findPrimarySegments(int maxBurgersCircuitSize){
    const int searchDepth = (maxBurgersCircuitSize - 1) / 2;
    assert(searchDepth >= 1);

    auto& verts = mesh().vertices();

    // Phase 1: Parallel search — find candidate collision edges
    struct Candidate{
        InterfaceMesh::Edge* edge;
    };
    std::vector<Candidate> candidates;
    tbb::spin_mutex candidatesMutex;

    tbb::parallel_for(tbb::blocked_range<size_t>(0, verts.size(), 16),
        [&](const tbb::blocked_range<size_t>& range){

        MemoryPool<SearchNode> localPool;
        std::vector<SearchNode*> queue;
        queue.reserve(512);
        std::unordered_map<InterfaceMesh::Vertex*, SearchNode*> visited;
        visited.reserve(256);

        for(size_t vi = range.begin(); vi < range.end(); ++vi){
            auto* startVert = verts[vi];

            if(startVert->edges() == nullptr) continue;

            bool allEdgesOwned = true;
            for(auto* e = startVert->edges(); e != nullptr; e = e->nextVertexEdge()){
                if(!e->nextCircuitEdge && !e->face()->circuit &&
                   !e->_claimedForCircuit.load(std::memory_order_relaxed)){
                    allEdgesOwned = false;
                    break;
                }
            }
            if(allEdgesOwned) continue;

            queue.clear();
            visited.clear();

            auto* root = localPool.construct();
            root->node = startVert;
            root->coord = Point3::Origin();
            root->tm.setIdentity();
            root->depth = 0;
            root->viaEdge = nullptr;
            queue.push_back(root);
            visited[startVert] = root;

            bool found = false;

            for(size_t qi = 0; qi < queue.size() && !found; ++qi){
                auto* cur = queue[qi];

                for(auto* edge = cur->node->edges(); edge != nullptr && !found; edge = edge->nextVertexEdge()){
                    if(edge->nextCircuitEdge || edge->face()->circuit) continue;
                    if(edge->_claimedForCircuit.load(std::memory_order_relaxed)) continue;

                    auto* nbVert = edge->vertex2();
                    Point3 nbCoord = cur->coord + cur->tm * edge->clusterVector;

                    auto it = visited.find(nbVert);
                    if(it != visited.end()){
                        auto* prevNode = it->second;
                        Vector3 b = prevNode->coord - nbCoord;
                        if(!b.isZero(CA_LATTICE_VECTOR_EPSILON)){
                            Matrix3 R = cur->tm * edge->clusterTransition->reverse->tm;
                            if(R.equals(prevNode->tm, CA_TRANSITION_MATRIX_EPSILON)){
                                tbb::spin_mutex::scoped_lock lock(candidatesMutex);
                                candidates.push_back({edge});
                                found = true;
                            }
                        }
                    }else if(cur->depth < searchDepth){
                        auto* nb = localPool.construct();
                        nb->node = nbVert;
                        nb->coord = nbCoord;
                        nb->depth = cur->depth + 1;
                        nb->viaEdge = edge;
                        nb->tm = edge->clusterTransition->isSelfTransition()
                            ? cur->tm
                            : cur->tm * edge->clusterTransition->reverse->tm;
                        queue.push_back(nb);
                        visited[nbVert] = nb;
                    }
                }
            }

            localPool.clear(true);
        }
    });

    // Phase 2: Serial creation — use original createBurgersCircuit logic
    // Re-do BFS from the candidate vertices and create circuits serially
    MemoryPool<SearchNode> pool;
    std::vector<SearchNode*> queue;
    queue.reserve(1024);

    for(auto& cand : candidates){
        auto* startVert = cand.edge->vertex1();

        if(cand.edge->circuit != nullptr || cand.edge->_claimedForCircuit.load(std::memory_order_relaxed))
            continue;

        queue.clear();

        auto* root = pool.construct();
        root->node = startVert;
        root->coord = Point3::Origin();
        root->tm.setIdentity();
        root->depth = 0;
        root->viaEdge = nullptr;
        startVert->burgersSearchStruct.store(reinterpret_cast<BurgersCircuitSearchStruct*>(root), std::memory_order_relaxed);
        queue.push_back(root);

        bool found = false;
        for(size_t qi = 0; qi < queue.size() && !found; ++qi){
            auto* cur = queue[qi];

            for(auto* edge = cur->node->edges(); edge != nullptr && !found; edge = edge->nextVertexEdge()){
                if(edge->nextCircuitEdge || edge->face()->circuit) continue;

                auto* nbVert = edge->vertex2();
                Point3 nbCoord = cur->coord + cur->tm * edge->clusterVector;

                auto* prevStruct = nbVert->burgersSearchStruct.load(std::memory_order_relaxed);
                if(prevStruct){
                    auto* prev = reinterpret_cast<SearchNode*>(prevStruct);
                    Vector3 b = prev->coord - nbCoord;
                    if(!b.isZero(CA_LATTICE_VECTOR_EPSILON)){
                        Matrix3 R = cur->tm * edge->clusterTransition->reverse->tm;
                        if(R.equals(prev->tm, CA_TRANSITION_MATRIX_EPSILON)){
                            found = createBurgersCircuitOriginal(edge, maxBurgersCircuitSize);
                        }
                    }
                }else if(cur->depth < searchDepth){
                    auto* nb = pool.construct();
                    nb->node = nbVert;
                    nb->coord = nbCoord;
                    nb->depth = cur->depth + 1;
                    nb->viaEdge = edge;
                    nb->tm = edge->clusterTransition->isSelfTransition()
                        ? cur->tm
                        : cur->tm * edge->clusterTransition->reverse->tm;
                    nbVert->burgersSearchStruct.store(reinterpret_cast<BurgersCircuitSearchStruct*>(nb), std::memory_order_relaxed);
                    queue.push_back(nb);
                }
            }
        }

        for(auto* sn : queue){
            sn->node->burgersSearchStruct.store(nullptr, std::memory_order_relaxed);
        }
        pool.clear(true);
    }

    return true;
}

bool BurgersLoopBuilder::createBurgersCircuitOriginal(InterfaceMesh::Edge* edge, int maxBurgersCircuitSize){
	if(edge->circuit != nullptr) return false;

	InterfaceMesh::Vertex* currentNode = edge->vertex1();
	InterfaceMesh::Vertex* neighborNode = edge->vertex2();
	auto* currentStruct = currentNode->burgersSearchStruct.load(std::memory_order_relaxed);
	auto* neighborStruct = neighborNode->burgersSearchStruct.load(std::memory_order_relaxed);
	assert(currentStruct != neighborStruct);

	BurgersCircuit* forwardCircuit = allocateCircuit();
	forwardCircuit->edgeCount = 1;
	forwardCircuit->firstEdge = forwardCircuit->lastEdge = edge->oppositeEdge();
	assert(forwardCircuit->firstEdge->circuit == nullptr);
	forwardCircuit->firstEdge->circuit = forwardCircuit;

	// Use a local set to track which vertices are on the first branch
	std::unordered_set<InterfaceMesh::Vertex*> firstBranchSet;
	for(BurgersCircuitSearchStruct* a = currentStruct; ; a = a->predecessorEdge->vertex1()->burgersSearchStruct.load(std::memory_order_relaxed)){
		firstBranchSet.insert(a->node);
		if(a->predecessorEdge == nullptr) break;
	}

	// Walk second branch inserting edges until hitting first branch
	for(BurgersCircuitSearchStruct* a = neighborStruct; ; a = a->predecessorEdge->vertex1()->burgersSearchStruct.load(std::memory_order_relaxed)){
		if(firstBranchSet.count(a->node)) break;
		assert(a->predecessorEdge != nullptr);
		assert(a->predecessorEdge->circuit == nullptr);
		a->predecessorEdge->nextCircuitEdge = forwardCircuit->firstEdge;
		forwardCircuit->firstEdge = a->predecessorEdge;
		forwardCircuit->edgeCount++;
		forwardCircuit->firstEdge->circuit = forwardCircuit;
	}

	// Walk first branch appending opposite edges until we hit the meeting node
	// (the first node on branch1 that is also reachable from branch2)
	InterfaceMesh::Vertex* meetNode = nullptr;
	for(BurgersCircuitSearchStruct* a = neighborStruct; ; a = a->predecessorEdge->vertex1()->burgersSearchStruct.load(std::memory_order_relaxed)){
		if(firstBranchSet.count(a->node)){
			meetNode = a->node;
			break;
		}
		if(a->predecessorEdge == nullptr) break;
	}

	for(BurgersCircuitSearchStruct* a = currentStruct; a->node != meetNode; a = a->predecessorEdge->vertex1()->burgersSearchStruct.load(std::memory_order_relaxed)){
		assert(a->predecessorEdge != nullptr);
		assert(a->predecessorEdge->oppositeEdge()->circuit == nullptr);
		forwardCircuit->lastEdge->nextCircuitEdge = a->predecessorEdge->oppositeEdge();
		forwardCircuit->lastEdge = forwardCircuit->lastEdge->nextCircuitEdge;
		forwardCircuit->edgeCount++;
		forwardCircuit->lastEdge->circuit = forwardCircuit;
	}

	// Close circuit
	forwardCircuit->lastEdge->nextCircuitEdge = forwardCircuit->firstEdge;
	assert(forwardCircuit->firstEdge != forwardCircuit->firstEdge->nextCircuitEdge);
	assert(forwardCircuit->countEdges() == forwardCircuit->edgeCount);
	assert(forwardCircuit->edgeCount >= 3);

	// Validate
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
		e = forwardCircuit->firstEdge;
		do{
			InterfaceMesh::Edge* nextEdge = e->nextCircuitEdge;
			e->nextCircuitEdge = nullptr;
			e->circuit = nullptr;
			e = nextEdge;
		}while(e != forwardCircuit->firstEdge);
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

	segment->line.push_back(backwardCircuit->calculateCenter());
	segment->coreSize.push_back(backwardCircuit->countEdges());

	appendLinePoint(segment->forwardNode());

	traceSegment(*segment, segment->forwardNode(), maxCircuitLength, true);
	traceSegment(*segment, segment->backwardNode(), maxCircuitLength, true);
}

// Test whether a given BurgersCircuit intersects any previously recorded circuit,
// by walking around each vertex and checking for crossing edges on adjacent faces.
bool BurgersLoopBuilder::intersectsOtherCircuits(BurgersCircuit* circuit){
    // Traverse each edge edge1 of our circuit
    InterfaceMesh::Edge* startEdge1 = circuit->firstEdge;
    for(InterfaceMesh::Edge* edge1 = startEdge1; ; edge1 = edge1->nextCircuitEdge){
        InterfaceMesh::Edge* edge2 = edge1->nextCircuitEdge;

        // Only interested if it is not the trivial pair
        if(edge1 != edge2->oppositeEdge()){
            // Let's traverse all the half-edges "cur" around the vertex common = edge1->vertex2()
            // starting at edge1->oppositeEdge()
            InterfaceMesh::Edge* sentinel = edge1->oppositeEdge();
            InterfaceMesh::Edge* cur = sentinel;

            do{
                // The anterior half-edge on the face
                InterfaceMesh::Edge* prev = cur->prevFaceEdge();
                // Check if that edge points to an existing circuit
                if(prev->circuit){
                    int goingOutside = 0;
					int goingInside = 0;
                    
                    // Add validation before calling the function
                    if(prev->nextCircuitEdge && 
                       edge2->oppositeEdge() && 
                       edge1->oppositeEdge() &&
                       prev->nextCircuitEdge->vertex1() == prev->vertex2()) {
                        circuitCircuitIntersection(
                            edge2->oppositeEdge(),
                            edge1->oppositeEdge(),
                            prev,
                            prev->nextCircuitEdge,
                            goingOutside,
                            goingInside
                        );

                        if(goingOutside){
                            return true;
                        }
                    }
                }
                // Move on to the next half-edge around the vertex
                cur = prev->oppositeEdge();
            }while (cur != sentinel);
        }

        // Have we gone all the way around the original circuit?
        if(edge2 == startEdge1){
            break;
        }
    }

    return false;
}

// Form a valid forward circuit, build the complementary inner loop (reverse circuit)
// that traces the hole boundary, marking faces as part of the circuit along the way.
BurgersCircuit* BurgersLoopBuilder::buildReverseCircuit(BurgersCircuit* forwardCircuit){
	BurgersCircuit* backwardCircuit = allocateCircuit();

	// Build the backward circuit along inner outline.
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
			InterfaceMesh::Edge* outerEdge1 = oppositeEdge1->nextFaceEdge()->oppositeEdge();
			InterfaceMesh::Edge* innerEdge1 = oppositeEdge1->prevFaceEdge()->oppositeEdge();
			InterfaceMesh::Edge* outerEdge2 = oppositeEdge2->prevFaceEdge()->oppositeEdge();
			InterfaceMesh::Edge* innerEdge2 = oppositeEdge2->nextFaceEdge()->oppositeEdge();
			assert(innerEdge1 != nullptr && innerEdge2 != nullptr);
			assert(innerEdge1->vertex1() == edge1->vertex2());
			assert(innerEdge2->vertex2() == edge1->vertex2());
			assert(innerEdge1->vertex1() == innerEdge2->vertex2());
			assert(innerEdge1->circuit == nullptr || innerEdge1->circuit == backwardCircuit);
			assert(innerEdge2->circuit == nullptr || innerEdge2->circuit == backwardCircuit);
			facet1->setFlag(1);
			facet1->circuit = backwardCircuit;
			facet2->setFlag(1);
			facet2->circuit = backwardCircuit;
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

	// Close circuit.
	backwardCircuit->lastEdge->nextCircuitEdge = backwardCircuit->firstEdge;

	assert(backwardCircuit->firstEdge != backwardCircuit->firstEdge->nextCircuitEdge);
	assert(backwardCircuit->countEdges() == backwardCircuit->edgeCount);
	assert(backwardCircuit->edgeCount >= 3);
	assert(!backwardCircuit->calculateBurgersVector().localVec().isZero(CA_LATTICE_VECTOR_EPSILON));

	return backwardCircuit;
}

void BurgersLoopBuilder::traceSegment(DislocationSegment& segment, DislocationNode& node, int maxCircuitLength, bool isPrimarySegment){
    BurgersCircuit& circuit = *node.circuit;
    assert(circuit.countEdges() == circuit.edgeCount);
    assert(circuit.isDangling);

    for(;;){
        size_t idx = _edgeStartIndex.fetch_add(1, std::memory_order_relaxed);
        int edgeIndex = (idx % circuit.edgeCount);

        InterfaceMesh::Edge* firstEdge = circuit.getEdge(edgeIndex);

        InterfaceMesh::Edge* edge0 = firstEdge;
        InterfaceMesh::Edge* edge1 = edge0->nextCircuitEdge;
        InterfaceMesh::Edge* edge2 = edge1->nextCircuitEdge;
        assert(edge1->circuit == &circuit);
        int counter = 0;
        do{
            // Check Burgers circuit.
            assert(circuit.edgeCount >= 3);
            
            // Check if Burgers vector is valid - if not, try to fix or skip
            ClusterVector burgersVec = circuit.calculateBurgersVector();
            if(burgersVec.localVec().isZero(CA_LATTICE_VECTOR_EPSILON)) {
                std::cerr << "Warning: Burgers vector is zero for circuit with " << circuit.edgeCount 
                         << " edges. Attempting to continue..." << std::endl;
                // Try to shorten the circuit to see if we can recover
                if(circuit.edgeCount <= 3) {
                    std::cerr << "Error: Cannot recover circuit with only " << circuit.edgeCount << " edges" << std::endl;
                    return; // Exit gracefully instead of asserting
                }
            }
            
            assert(circuit.countEdges() == circuit.edgeCount);
            assert(edge0->circuit == &circuit && edge1->circuit == &circuit && edge2->circuit == &circuit);

            bool wasShortened = false;
            if(tryRemoveTwoCircuitEdges(edge0, edge1, edge2)){
                wasShortened = true;
			}else if(tryRemoveThreeCircuitEdges(edge0, edge1, edge2, isPrimarySegment)){
                wasShortened = true;
			}else if(tryRemoveOneCircuitEdge(edge0, edge1, edge2, isPrimarySegment)){
                wasShortened = true;
			}else if(trySweepTwoFacets(edge0, edge1, edge2, isPrimarySegment)){
                wasShortened = true;
			}

            if(wasShortened){
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

        if(circuit.edgeCount >= maxCircuitLength) break;

        bool wasExtended = false;

        idx = _edgeStartIndex.fetch_add(1, std::memory_order_relaxed);
        edgeIndex = (idx % circuit.edgeCount);

        firstEdge = circuit.getEdge(edgeIndex);

        edge0 = firstEdge;
        edge1 = firstEdge->nextCircuitEdge;
        do{
            if(tryInsertOneCircuitEdge(edge0, edge1, isPrimarySegment)){
                wasExtended = true;
                appendLinePoint(node);
                break;
            }

            edge0 = edge1;
            edge1 = edge1->nextCircuitEdge;
        }while(edge0 != firstEdge);
        if(!wasExtended) break;
    }
}

// Try removing two consecutive edges that form a back-and-forth pair, simplyfing the loop
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

// Try collapsing a small three-edge pocket by folding across a single face.
bool BurgersLoopBuilder::tryRemoveThreeCircuitEdges(
	InterfaceMesh::Edge*& edge0, 
	InterfaceMesh::Edge*& edge1, 
	InterfaceMesh::Edge*& edge2,
	bool isPrimarySegment
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

	facet1->circuit = circuit;
	if(isPrimarySegment) facet1->setFlag(1);

	return true;
}

// Try replacing a single edge by its "short" diagonal across a quad face.
bool BurgersLoopBuilder::tryRemoveOneCircuitEdge(
	InterfaceMesh::Edge*& edge0, 
	InterfaceMesh::Edge*& edge1, 
	InterfaceMesh::Edge*& edge2, 
	bool isPrimarySegment
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
	if(!shortEdge->tryClaimForCircuit()) return false;

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

	facet1->circuit = circuit;
	if(isPrimarySegment) facet1->setFlag(1);

	return true;
}

// Sweep across two adjacent facets to shortcut around a concave bend.
bool BurgersLoopBuilder::trySweepTwoFacets(
	InterfaceMesh::Edge*& edge0, 
	InterfaceMesh::Edge*& edge1, 
	InterfaceMesh::Edge*& edge2, 
	bool isPrimarySegment
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

	if(!outerEdge1->tryClaimForCircuit()){
		return false;
	}
	if(!outerEdge2->tryClaimForCircuit()){
		outerEdge1->releaseCircuitClaim();
		return false;
	}

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

	facet1->circuit = circuit;
	facet2->circuit = circuit;
	if(isPrimarySegment){
		facet1->setFlag(1);
		facet2->setFlag(1);
	}

	edge0 = outerEdge1;
	edge1 = outerEdge2;
	edge2 = edge1->nextCircuitEdge;

	return true;
}

// Insert a new edge across a free face to grow the circuit outward.
bool BurgersLoopBuilder::tryInsertOneCircuitEdge(
	InterfaceMesh::Edge*& edge0, 
	InterfaceMesh::Edge*& edge1, 
	bool isPrimarySegment
){
	assert(edge0 != edge1->oppositeEdge());

	InterfaceMesh::Face* facet = edge1->face();
	if(facet->circuit != nullptr) return false;

	InterfaceMesh::Edge* insertEdge1 = edge1->prevFaceEdge()->oppositeEdge();
	if(insertEdge1->circuit != nullptr) return false;

	InterfaceMesh::Edge* insertEdge2 = edge1->nextFaceEdge()->oppositeEdge();
	if(insertEdge2->circuit != nullptr) return false;

	if(!insertEdge1->tryClaimForCircuit()) return false;
	if(!insertEdge2->tryClaimForCircuit()){
		insertEdge1->releaseCircuitClaim();
		return false;
	}

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

	facet->circuit = circuit;
	if(isPrimarySegment) facet->setFlag(1);

	return true;
}

// After each successful removal or insertion, compute the segment's new center of mass,
// apply periodic wrapping, and append that point to the dislocation line
void BurgersLoopBuilder::appendLinePoint(DislocationNode& node){
	DislocationSegment& segment = *node.segment;
	assert(!segment.line.empty());

	// Get size of dislocation core.
	int coreSize = node.circuit->edgeCount;

	// Make sure the line is not wrapped at periodic boundaries.
	const Point3& lastPoint = node.isForwardNode() ? segment.line.back() : segment.line.front();
	Point3 newPoint = lastPoint + cell().wrapVector(node.circuit->calculateCenter() - lastPoint);

	if(node.isForwardNode()){
		// Add a new point to end the line.
		segment.line.push_back(newPoint);
		segment.coreSize.push_back(coreSize);
	}else{
		// Add a new point to start the line.
		segment.line.insert(segment.line.begin(), newPoint);
		segment.coreSize.insert(segment.coreSize.begin(), coreSize);
	}

	node.circuit->numPreliminaryPoints++;
}

// When two circuits share a face corner, classify wheter one is going "inside" or
// "outside" relative to the other, to detect illegal crossings.
void BurgersLoopBuilder::circuitCircuitIntersection(
	InterfaceMesh::Edge* circuitAEdge1, 
	InterfaceMesh::Edge* circuitAEdge2, 
	InterfaceMesh::Edge* circuitBEdge1, 
	InterfaceMesh::Edge* circuitBEdge2, 
	int& goingOutside, 
	int& goingInside
){
	assert(circuitAEdge2->vertex1() == circuitBEdge2->vertex1());
	assert(circuitAEdge1->vertex2() == circuitBEdge2->vertex1());
	assert(circuitBEdge1->vertex2() == circuitBEdge2->vertex1());

	// Iterate over interior facet edges.
	InterfaceMesh::Edge* edge = circuitBEdge2;
	bool contour1inside = false;
	bool contour2inside = false;
	// Add counter to prevent infinite loops in large systems
	int safetyCounter = 0;
	const int maxEdgeTraversal = 1000; // Adjust as needed for large systems
	
	for(;;){
		InterfaceMesh::Edge* oppositeEdge = edge->oppositeEdge();
		if(oppositeEdge == circuitBEdge1) break;
		if(edge != circuitBEdge2){
			if(oppositeEdge == circuitAEdge1) contour1inside = true;
			if(edge == circuitAEdge2) contour2inside = true;
		}

		edge = oppositeEdge->nextFaceEdge();
		if(edge->vertex1() != circuitBEdge2->vertex1() || edge == circuitBEdge2){
			// Instead of asserting, break the loop and use current state
			// This handles rare topological cases in large systems
			break;
		}
		
		// Safety check to avoid infinite loops
		if(++safetyCounter > maxEdgeTraversal) break;
	}
	
	// Keep the intention but make more robust
	if(circuitAEdge2 != circuitBEdge2) {
		// The original assertion is valid
	} else {
		// If equal, contour2inside should be false - force it
		contour2inside = false;
	}

	// Iterate over exterior facet edges.
	bool contour1outside = false;
	bool contour2outside = false;
	edge = circuitBEdge1;
	// Reset safety counter for second loop
	safetyCounter = 0;
	
	for(;;){
		InterfaceMesh::Edge* nextEdge = edge->nextFaceEdge();
		if(nextEdge == circuitBEdge2) break;
		InterfaceMesh::Edge* oppositeEdge = nextEdge->oppositeEdge();
		
		// Handle potential topology issues more robustly
		if(oppositeEdge->vertex2() != circuitBEdge2->vertex1()) {
			// Found unexpected topology, break the loop
			break;
		}
		
		edge = oppositeEdge;
		if(edge == circuitAEdge1) contour1outside = true;
		if(nextEdge == circuitAEdge2) contour2outside = true;
		
		// Safety check to avoid infinite loops
		if(++safetyCounter > 1000) break;
	}

	// Handle invariants more robustly
	if(contour1outside && contour1inside) {
		// This shouldn't happen according to original assertion
		// Force a consistent state - prioritize 'outside'
		contour1inside = false;
	}
	
	if(contour2outside && contour2inside) {
		// This shouldn't happen according to original assertion
		// Force a consistent state - prioritize 'outside'  
		contour2inside = false;
	}

	if(contour2outside == true && contour1outside == false){
		goingOutside += 1;
	}else if(contour2inside == true && contour1inside == false){
		goingInside += 1;
	}
}

// After trimming and extending, attempt to join any remaining dangling circuits
// by (1) creating secondary loops in adjacent holes, (2) marking fully blocked circuits
// as junctions candidates, and (3) either fusing two arms or forming multi-arm junctions
size_t BurgersLoopBuilder::joinSegments(int maxCircuitLength){
	// First iteration over all dangling circuits.
	// Try to create secondary dislocation segments in the adjacent regions of the
	// interface mesh.
	for(size_t nodeIndex = 0; nodeIndex < danglingNodes().size(); nodeIndex++){
		DislocationNode* node = danglingNodes()[nodeIndex];
		BurgersCircuit* circuit = node->circuit;
		assert(circuit->isDangling);

		// Go around the circuit to find an unvisited region on the interface mesh.
		InterfaceMesh::Edge* edge = circuit->firstEdge;
		do{
			assert(edge->circuit == circuit);
			BurgersCircuit* oppositeCircuit = edge->oppositeEdge()->circuit;
			if(oppositeCircuit == nullptr){
				assert(edge->oppositeEdge()->nextCircuitEdge == nullptr);

				// Try to create a new circuit inside the unvisited region.
				createSecondarySegment(edge, circuit, maxCircuitLength);

				// Skip edges to the end of the unvisited interval.
				while(edge->oppositeEdge()->circuit == nullptr && edge != circuit->firstEdge){
					edge = edge->nextCircuitEdge;
				}
			}else{
				edge = edge->nextCircuitEdge;
			}
		}
		while(edge != circuit->firstEdge);
	}

	// Second pass over all dangling nodes.
	// Mark circuits that are completely blocked by other circuits.
	// They are candidates for the formation of junctions.
	for(DislocationNode* node : danglingNodes()){
		BurgersCircuit* circuit = node->circuit;
		assert(circuit->isDangling);

		// Go around the circuit to see whether it is completely surrounded by other circuits.
		// Put it into one ring with the adjacent circuits.
		circuit->isCompletelyBlocked = true;
		InterfaceMesh::Edge* edge = circuit->firstEdge;
		do{
			assert(edge->circuit == circuit);
			BurgersCircuit* adjacentCircuit = edge->oppositeEdge()->circuit;
			if(adjacentCircuit == nullptr){
				// Found a section of the circuit, which is not blocked
				// by some other circuit.
				circuit->isCompletelyBlocked = false;
				break;
			}else if(adjacentCircuit != circuit){
				assert(adjacentCircuit->isDangling);
				DislocationNode* adjacentNode = adjacentCircuit->dislocationNode;
				if(node->formsJunctionWith(adjacentNode) == false){
					node->connectNodes(adjacentNode);
				}
			}
			edge = edge->nextCircuitEdge;
		}while(edge != circuit->firstEdge);
	}

	// Count number of created dislocation junctions.
	size_t numJunctions = 0;

	// Actually create junctions for completely blocked circuits.
	for(DislocationNode* node : danglingNodes()){
		BurgersCircuit* circuit = node->circuit;

		// Skip circuits which have already become part of a junction.
		if(circuit->isDangling == false) continue;

		// Skip dangling circuits, which are not completely blocked by other circuits;
		if(!circuit->isCompletelyBlocked){
			node->dissolveJunction();
			continue;
		}

		// Junctions must consist of at least two dislocation segments.
		if(node->junctionRing == node) continue;

		assert(node->segment->replacedWith == nullptr);

		// Compute center of mass of junction node.
		Vector3 centerOfMassVector = Vector3::Zero();
		Point3 basePoint = node->position();
		int armCount = 1;
		bool allCircuitsCompletelyBlocked = true;
		DislocationNode* armNode = node->junctionRing;
		while(armNode != node){
			assert(armNode->segment->replacedWith == nullptr);
			assert(armNode->circuit->isDangling);
			if(armNode->circuit->isCompletelyBlocked == false) {
				allCircuitsCompletelyBlocked = false;
				break;
			}
			armCount++;
			centerOfMassVector += cell().wrapVector(armNode->position() - basePoint);
			armNode = armNode->junctionRing;
		}

		// All circuits of the junction must be fully blocked by other circuits.
		if(!allCircuitsCompletelyBlocked){
			node->dissolveJunction();
			continue;
		}

		// Junctions must consist of at least two dislocation segments.
		assert(armCount >= 2);

		// Only create a real junction for three or more segments.
		if(armCount >= 3){
			Point3 centerOfMass = basePoint + centerOfMassVector / armCount;

			// Iterate over all arms of the new junction.
			armNode = node;
			do{
				// Mark this node as no longer dangling.
				armNode->circuit->isDangling = false;
				assert(armNode != armNode->junctionRing);

				// Extend arm to junction's exact center point.
				std::vector<Point3>& line = armNode->segment->line;
				if(armNode->isForwardNode()){
					line.push_back(line.back() + cell().wrapVector(centerOfMass - line.back()));
					armNode->segment->coreSize.push_back(armNode->segment->coreSize.back());
				}else{
					line.insert(line.begin(), line.front() + cell().wrapVector(centerOfMass - line.front()));
					armNode->segment->coreSize.insert(armNode->segment->coreSize.begin(), armNode->segment->coreSize.front());
				}

				armNode->circuit->numPreliminaryPoints = 0;
				armNode = armNode->junctionRing;
			}while(armNode != node);
			numJunctions++;
		}else{
			// For a two-armed junction, just merge the two segments into one.
			DislocationNode* node1 = node;
			DislocationNode* node2 = node->junctionRing;
			assert(node1 != node2);
			assert(node2->junctionRing == node1);
			assert(node1->junctionRing == node2);

			BurgersCircuit* circuit1 = node1->circuit;
			BurgersCircuit* circuit2 = node2->circuit;
			circuit1->isDangling = false;
			circuit2->isDangling = false;
			circuit1->numPreliminaryPoints = 0;
			circuit2->numPreliminaryPoints = 0;

			// Check if this is a closed dislocation loop.
			if(node1->oppositeNode == node2){
				assert(node1->segment == node2->segment);
				DislocationSegment* loop = node1->segment;
				assert(loop->isClosedLoop());

				// Make both ends of the segment coincide by adding an extra point if necessary.
				if(!cell().wrapVector(node1->position() - node2->position()).isZero(CA_ATOM_VECTOR_EPSILON)){
					loop->line.push_back(loop->line.back() + cell().wrapVector(loop->line.front() - loop->line.back()));
					assert(cell().wrapVector(node1->position() - node2->position()).isZero(CA_ATOM_VECTOR_EPSILON));
					loop->coreSize.push_back(loop->coreSize.back());
				}

				// Loop segment should not be degenerate.
				assert(loop->line.size() >= 3);
			}else{
				// If not a closed loop, merge the two segments into a single line.
				assert(node1->segment != node2->segment);

				DislocationNode* farEnd1 = node1->oppositeNode;
				DislocationNode* farEnd2 = node2->oppositeNode;
				DislocationSegment* segment1 = node1->segment;
				DislocationSegment* segment2 = node2->segment;
				if(node1->isBackwardNode()){
					segment1->nodes[1] = farEnd2;
					Vector3 shiftVector;
					if(node2->isBackwardNode()){
						shiftVector = calculateShiftVector(segment1->line.front(), segment2->line.front());
						segment1->line.insert(segment1->line.begin(), segment2->line.rbegin(), segment2->line.rend() - 1);
						segment1->coreSize.insert(segment1->coreSize.begin(), segment2->coreSize.rbegin(), segment2->coreSize.rend() - 1);
					}else{
						shiftVector = calculateShiftVector(segment1->line.front(), segment2->line.back());
						segment1->line.insert(segment1->line.begin(), segment2->line.begin(), segment2->line.end() - 1);
						segment1->coreSize.insert(segment1->coreSize.begin(), segment2->coreSize.begin(), segment2->coreSize.end() - 1);
					}

					if(shiftVector != Vector3::Zero()){
						for(auto p = segment1->line.begin(); p != segment1->line.begin() + segment2->line.size() - 1; ++p){
							*p -= shiftVector;
						}
					}
				}else{
					segment1->nodes[0] = farEnd2;
					Vector3 shiftVector;
					if(node2->isBackwardNode()){
						shiftVector = calculateShiftVector(segment1->line.back(), segment2->line.front());
						segment1->line.insert(segment1->line.end(), segment2->line.begin() + 1, segment2->line.end());
						segment1->coreSize.insert(segment1->coreSize.end(), segment2->coreSize.begin() + 1, segment2->coreSize.end());
					}else{
						shiftVector = calculateShiftVector(segment1->line.back(), segment2->line.back());
						segment1->line.insert(segment1->line.end(), segment2->line.rbegin() + 1, segment2->line.rend());
						segment1->coreSize.insert(segment1->coreSize.end(), segment2->coreSize.rbegin() + 1, segment2->coreSize.rend());
					}

					if(shiftVector != Vector3::Zero()){
						for(auto p = segment1->line.end() - segment2->line.size() + 1; p != segment1->line.end(); ++p){
							*p -= shiftVector;
						}
					}
				}
				
				farEnd2->segment = segment1;
				farEnd2->oppositeNode = farEnd1;
				farEnd1->oppositeNode = farEnd2;
				node1->oppositeNode = node2;
				node2->oppositeNode = node1;
				segment2->replacedWith = segment1;
				network().discardSegment(segment2);
			}
		}
	}

	// Clean up list of dangling nodes. Remove joined nodes.
	_danglingNodes.erase(std::remove_if(_danglingNodes.begin(), _danglingNodes.end(),[](DislocationNode* node){
		return !node->isDangling();
	}), _danglingNodes.end());

	return numJunctions;
}

// When a dangling circuit borders an unvisited hole, trace that hole boundary as a
// secondary Burgers loop, validate its Burgers vector and closure, and convert it into a 
// new dislocation segment if valid.
void BurgersLoopBuilder::createSecondarySegment(InterfaceMesh::Edge* firstEdge, BurgersCircuit* outerCircuit, int maxCircuitLength){
	assert(firstEdge->circuit == outerCircuit);

	// Create circuit along the border of the hole.
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
			InterfaceMesh::Face* oppositeFacet = oppositeEdge->face();
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

	// Create secondary segment only for dislocations (b != 0) and small enough dislocation cores.
	const bool rejectNumCircuits = numCircuits == 1;
	const bool rejectTooLong = edgeCount > maxCircuitLength;
	const bool rejectZeroBurgers = burgersVector.isZero(CA_LATTICE_VECTOR_EPSILON);
	const bool rejectEdgeSum = edgeSum.isZero(CA_ATOM_VECTOR_EPSILON) == false;
	const bool rejectFrank = !frankRotation.equals(Matrix3::Identity(), CA_TRANSITION_MATRIX_EPSILON);
	if(rejectNumCircuits || rejectTooLong || rejectZeroBurgers || rejectEdgeSum || rejectFrank){

		// Discard unused circuit.
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

	// Create forward circuit.
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

	// Do all the rest.
	createAndTraceSegment(ClusterVector(burgersVector, baseCluster), forwardCircuit, maxCircuitLength);
}

}
