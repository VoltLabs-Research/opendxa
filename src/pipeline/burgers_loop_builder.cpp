#include <volt/pipeline/burgers_loop_builder.h>
#include <volt/pipeline/interface_mesh.h>
#include <volt/helpers/triangle_intersection.h>
#include <volt/core/progress.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/blocked_range.h>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <numeric>
#include <optional>
#include <span>
#include <ranges>
#include <array>
#include <algorithm>
#include <utility>

namespace Volt{

// Records the facet's corners when the serial primary pass is collecting the vertices it changes.
void BurgersLoopBuilder::sweepFacet(InterfaceMesh::Face* facet, BurgersCircuit* circuit, bool markSwept){
    facet->circuit = circuit;
    if(markSwept) facet->setFlag(InterfaceMeshFace::kSweptBeforeCapFlag);
    if(_changedVertices){
        auto* edge = facet->edges();
        do{
            _changedVertices->push_back(edge->vertex1());
            edge = edge->nextFaceEdge();
        }while(edge != facet->edges());
    }
}

// From the trial size on, each circuit that still dangles stores its cap once: the defect mesh closes the circuit's
// hole with it, and the line drops the points traced after it.
void BurgersLoopBuilder::storeCaps(const std::vector<DislocationNode*>& nodes){
    tbb::parallel_for_each(nodes.begin(), nodes.end(), [](DislocationNode* node){
        BurgersCircuit* circuit = node->circuit;
        if(circuit->isDangling && circuit->segmentMeshCap.empty()){
            circuit->storeCircuit();
            circuit->numPreliminaryPoints = 0;
        }
    });
}


void BurgersLoopBuilder::traceDislocationSegments(){
    mesh().clearFaceFlag(InterfaceMeshFace::kSweptBeforeCapFlag);
    if(_markCoreAtoms){
        indexDefectCells();
    }

    // One step per circuit length, from 3 up to the stretched maximum.
    ProgressReporter progress("Tracing Burgers circuits",
        static_cast<std::uint64_t>(std::max(0, _maxExtendedCircuitSize - 2)));
    for(int maxCircuitLength : std::views::iota(3, _maxExtendedCircuitSize + 1)){
        const bool withinTrialSize = maxCircuitLength <= _maxTrialCircuitSize;
        const std::vector<DislocationNode*> dangling = _danglingNodes;
        for(DislocationNode* node : dangling){
            traceSegment(*node, maxCircuitLength, withinTrialSize);
        }

        // An even length has the same search depth, (n - 1) / 2, as the odd length before it, so it is skipped.
        if((maxCircuitLength & 1) && withinTrialSize){
            findPrimarySegments(maxCircuitLength);
        }
        joinSegments(maxCircuitLength);

        if(maxCircuitLength >= _maxTrialCircuitSize){
            storeCaps(dangling);
        }
        progress.update(static_cast<std::uint64_t>(maxCircuitLength - 2));
    }
}

void BurgersLoopBuilder::indexDefectCells(){
    const double alpha = DefectCellIndex::kAlphaScale * _mesh.structureAnalysis().maximumNeighborDistance();
    _defectCells.emplace(_mesh.tessellation(), alpha);
    _coreCellClaims = std::vector<CoreCellClaim>(_defectCells->cellCount());
}

// Drops the points each end traced after its circuit stored its cap, so that a dangling line stops at its cap.
void BurgersLoopBuilder::trimPreliminaryPoints(DislocationSegment& segment){
    segment.trimEnds(segment.backwardNode().circuit->numPreliminaryPoints, segment.forwardNode().circuit->numPreliminaryPoints);
}

// Expresses b in a cluster of the reference structure when one transition links the segment's cluster to such a cluster.
void BurgersLoopBuilder::expressInReferenceStructure(DislocationSegment& segment, std::string_view referenceTopologyName){
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
void BurgersLoopBuilder::orientCanonically(DislocationSegment& segment){
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

// The search walks only edges that no circuit runs along and whose facet no circuit has swept.
bool BurgersLoopBuilder::isFreeEdge(const InterfaceMesh::Edge* edge){
    return !edge->nextCircuitEdge && !edge->face()->circuit;
}

bool BurgersLoopBuilder::hasFreeEdge(const InterfaceMesh::Vertex* vertex){
    for(auto* edge = vertex->edges(); edge != nullptr; edge = edge->nextVertexEdge()){
        if(isFreeEdge(edge)) return true;
    }
    return false;
}

BurgersLoopBuilder::SearchNode* BurgersLoopBuilder::CircuitSearch::newNode(){
    const std::size_t block = _usedNodes / kBlockSize;
    if(block == _blocks.size()){
        _blocks.push_back(std::make_unique<SearchNode[]>(kBlockSize));
    }
    return &_blocks[block][_usedNodes++ % kBlockSize];
}

// A new generation empties the reached table without touching it; on the rare wrap-around it is cleared for real.
void BurgersLoopBuilder::CircuitSearch::startSearch(){
    _queue.clear();
    _usedNodes = 0;
    if(++_generation == 0){
        std::ranges::fill(_slots, Slot{});
        _generation = 1;
    }
}

BurgersLoopBuilder::SearchNode* BurgersLoopBuilder::CircuitSearch::reached(int key) const{
    const std::size_t mask = _slots.size() - 1;
    for(std::size_t i = slotOf(key); ; i = (i + 1) & mask){
        const Slot& slot = _slots[i];
        if(slot.generation != _generation) return nullptr;
        if(slot.vertex == key) return slot.node;
    }
}

// Only vertices the search has not reached are enqueued, so a key is never inserted twice.
void BurgersLoopBuilder::CircuitSearch::insertSlot(int vertex, SearchNode* node){
    const std::size_t mask = _slots.size() - 1;
    for(std::size_t i = slotOf(vertex); ; i = (i + 1) & mask){
        Slot& slot = _slots[i];
        if(slot.generation != _generation){
            slot = Slot{_generation, vertex, node};
            return;
        }
    }
}

void BurgersLoopBuilder::CircuitSearch::growSlots(){
    std::vector<Slot> old = std::exchange(_slots, std::vector<Slot>(_slots.size() * 2));
    for(const Slot& slot : old){
        if(slot.generation == _generation) insertSlot(slot.vertex, slot.node);
    }
}

// Keeps the table at most half full.
void BurgersLoopBuilder::CircuitSearch::enqueue(SearchNode* node){
    _queue.push_back(node);
    if(2 * _queue.size() > _slots.size()) growSlots();
    insertSlot(node->vertexIndex, node);
}

void BurgersLoopBuilder::initRoot(SearchNode* root, InterfaceMesh::Vertex* vertex){
    root->vertex = vertex;
    root->vertexIndex = vertex->index();
    root->coord = Point3::Origin();
    root->tm.setIdentity();
    root->depth = 0;
    root->viaEdge = nullptr;
    root->parent = nullptr;
}

void BurgersLoopBuilder::initChild(SearchNode* child, SearchNode* parent, const FreeEdge& edge, const Point3& coord){
    child->vertex = edge.target;
    child->vertexIndex = edge.targetIndex;
    child->coord = coord;
    child->tm = edge.transition->isSelfTransition()
        ? parent->tm
        : parent->tm * edge.transition->reverse->tm;
    child->depth = parent->depth + 1;
    child->viaEdge = edge.edge;
    child->parent = parent;
}

// Arriving at a vertex the search already reached closes a circuit. It encloses a dislocation when the two paths to
// the vertex agree on the lattice orientation but not on the lattice coordinate.
bool BurgersLoopBuilder::closesNonZeroCircuit(const SearchNode* from, const FreeEdge& edge, const Point3& coord, const SearchNode* reached){
    const Vector3 b = reached->coord - coord;
    if(b.isZero(CA_LATTICE_VECTOR_EPSILON)) return false;
    const Matrix3 R = from->tm * edge.transition->reverse->tm;
    return R.equals(reached->tm, CA_TRANSITION_MATRIX_EPSILON);
}

// The free edges leaving a vertex, read from the mesh as it is now.
struct BurgersLoopBuilder::LiveFreeEdges{
    template<typename Visit>
    void forEach(const SearchNode* from, Visit&& visit) const{
        for(auto* edge = from->vertex->edges(); edge != nullptr; edge = edge->nextVertexEdge()){
            if(!isFreeEdge(edge)) continue;
            if(visit(FreeEdge{edge->clusterVector, edge->clusterTransition, edge, edge->vertex2(), edge->vertex2()->index()})) return;
        }
    }
};

// The free edges leaving each vertex, copied once into one array in the order of each vertex's edge list, for a
// phase that searches from every vertex while nothing changes the mesh. A search then reads its edges from contiguous
// memory instead of chasing each edge, its facet and its end vertex through the mesh, with the same values in the
// same order, so it finds what it would find on the mesh.
class BurgersLoopBuilder::FrozenFreeEdges{
public:
    explicit FrozenFreeEdges(const InterfaceMesh& mesh){
        const auto& vertices = mesh.vertices();
        _begin.assign(vertices.size() + 1, 0);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, vertices.size(), 4096), [&](const tbb::blocked_range<size_t>& r){
            for(size_t i = r.begin(); i < r.end(); ++i){
                std::uint32_t count = 0;
                for(auto* edge = vertices[i]->edges(); edge != nullptr; edge = edge->nextVertexEdge()){
                    if(isFreeEdge(edge)) ++count;
                }
                _begin[i + 1] = count;
            }
        });
        std::inclusive_scan(_begin.begin(), _begin.end(), _begin.begin());
        _edges.resize(_begin.back());
        tbb::parallel_for(tbb::blocked_range<size_t>(0, vertices.size(), 4096), [&](const tbb::blocked_range<size_t>& r){
            for(size_t i = r.begin(); i < r.end(); ++i){
                std::size_t slot = _begin[i];
                for(auto* edge = vertices[i]->edges(); edge != nullptr; edge = edge->nextVertexEdge()){
                    if(isFreeEdge(edge)) _edges[slot++] = FreeEdge{edge->clusterVector, edge->clusterTransition, edge, edge->vertex2(), edge->vertex2()->index()};
                }
            }
        });
    }

    bool hasFreeEdge(std::size_t vertex) const{
        return _begin[vertex + 1] > _begin[vertex];
    }

    template<typename Visit>
    void forEach(const SearchNode* from, Visit&& visit) const{
        const auto index = static_cast<std::size_t>(from->vertexIndex);
        for(std::size_t i = _begin[index]; i < _begin[index + 1]; ++i){
            if(visit(_edges[i])) return;
        }
    }

private:
    std::vector<std::size_t> _begin;
    std::vector<FreeEdge> _edges;
};

// Calls onClosure(edge, from, reached) for each edge that closes a circuit with a non-zero Burgers vector, and
// stops when it returns true.
template<typename Edges, typename OnClosure>
void BurgersLoopBuilder::CircuitSearch::run(const Edges& freeEdges, InterfaceMesh::Vertex* start, int maxDepth, OnClosure&& onClosure){
    startSearch();
    SearchNode* root = newNode();
    initRoot(root, start);
    enqueue(root);

    bool stopped = false;
    for(size_t i = 0; i < _queue.size() && !stopped; ++i){
        SearchNode* from = _queue[i];
        freeEdges.forEach(from, [&](const FreeEdge& edge){
            const Point3 coord = from->coord + from->tm * edge.clusterVector;
            if(SearchNode* other = reached(edge.targetIndex)){
                if(closesNonZeroCircuit(from, edge, coord, other) && onClosure(edge.edge, from, other)){
                    stopped = true;
                    return true;
                }
            }else if(from->depth < maxDepth){
                SearchNode* child = newNode();
                initChild(child, from, edge, coord);
                enqueue(child);
            }
            return false;
        });
    }
}

BurgersLoopBuilder::CircuitSearch& BurgersLoopBuilder::threadSearch(){
    auto& search = _searches.local();
    if(!search) search = std::make_unique<CircuitSearch>();
    return *search;
}

// The first closing edge found from each vertex, sorted so that the serial pass creates circuits in a fixed order
// rather than in the order the threads found them: each circuit it creates changes which edges are free. Each vertex
// writes its own slot, and the sort is stable over vertex order, so edges with the same end vertices keep the order
// of the vertices they were found from, whatever the thread count.
std::vector<InterfaceMesh::Edge*> BurgersLoopBuilder::findClosingEdges(const FrozenFreeEdges& freeEdges, int searchDepth){
    const auto& vertices = mesh().vertices();
    std::vector<InterfaceMesh::Edge*> found(vertices.size(), nullptr);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, vertices.size(), 16), [&](const tbb::blocked_range<size_t>& range){
        CircuitSearch& search = threadSearch();
        for(size_t i = range.begin(); i < range.end(); ++i){
            if(!freeEdges.hasFreeEdge(i)) continue;
            search.run(freeEdges, vertices[i], searchDepth, [&](InterfaceMesh::Edge* edge, SearchNode*, SearchNode*){
                found[i] = edge;
                return true;
            });
        }
    });

    std::vector<InterfaceMesh::Edge*> closingEdges;
    for(InterfaceMesh::Edge* edge : found){
        if(edge) closingEdges.push_back(edge);
    }
    std::ranges::stable_sort(closingEdges, [](const InterfaceMesh::Edge* a, const InterfaceMesh::Edge* b){
        const int a1 = a->vertex1()->index();
        const int b1 = b->vertex1()->index();
        if(a1 != b1) return a1 < b1;
        return a->vertex2()->index() < b->vertex2()->index();
    });
    return closingEdges;
}

// How the crossing test sees the circuits on the mesh: as linked into it, with a candidate circuit that is not linked
// laid over them. The candidate's edges are few (at most the trial circuit size), so they are looked up in order.
struct BurgersLoopBuilder::CandidateLinks{
    std::span<InterfaceMesh::Edge* const> candidate;

    std::optional<std::size_t> position(const InterfaceMesh::Edge* edge) const{
        for(std::size_t i = 0; i < candidate.size(); ++i){
            if(candidate[i] == edge) return i;
        }
        return std::nullopt;
    }

    bool onCircuit(const InterfaceMesh::Edge* edge) const{
        return edge->circuit != nullptr || position(edge).has_value();
    }

    InterfaceMesh::Edge* next(const InterfaceMesh::Edge* edge) const{
        if(const auto i = position(edge)) return candidate[(*i + 1) % candidate.size()];
        return edge->nextCircuitEdge;
    }
};

// Vertices within a given number of edges, either way, of a vertex the serial pass changed. A search from a vertex
// reads only the edges and facets of the vertices it reaches, all within its depth of it, so its look-ahead verdict
// still holds as long as its start vertex is outside this region.
class BurgersLoopBuilder::ChangedRegion{
public:
    explicit ChangedRegion(std::size_t vertexCount) : _reach(vertexCount, -1){}

    bool contains(const InterfaceMesh::Vertex* vertex) const{
        return _reach[static_cast<std::size_t>(vertex->index())] >= 0;
    }

    void spread(const std::vector<InterfaceMesh::Vertex*>& changed, int depth){
        _queue.clear();
        for(InterfaceMesh::Vertex* vertex : changed) relax(vertex, depth);
        for(std::size_t i = 0; i < _queue.size(); ++i){
            InterfaceMesh::Vertex* vertex = _queue[i];
            const int remaining = _reach[static_cast<std::size_t>(vertex->index())] - 1;
            if(remaining < 0) continue;
            for(auto* edge = vertex->edges(); edge != nullptr; edge = edge->nextVertexEdge()){
                relax(edge->vertex2(), remaining);
                relax(edge->prevFaceEdge()->vertex1(), remaining);
            }
        }
    }

private:
    void relax(InterfaceMesh::Vertex* vertex, int remaining){
        auto& reach = _reach[static_cast<std::size_t>(vertex->index())];
        if(reach >= remaining) return;
        reach = remaining;
        _queue.push_back(vertex);
    }

    std::vector<int> _reach;
    std::vector<InterfaceMesh::Vertex*> _queue;
};

// Disjoint sets over the dangling nodes that linkAdjacentCircuits keeps equal to their junction rings, so that asking
// whether two nodes share a ring does not walk it. Rings grow to thousands of nodes along a grain boundary full of
// dangling circuits, and walking one for each edge two circuits share made linking quadratic in the ring size.
class BurgersLoopBuilder::RingSets{
public:
    explicit RingSets(std::size_t count) : _parent(count), _size(count, 1){
        std::iota(_parent.begin(), _parent.end(), std::size_t{0});
    }

    std::size_t find(std::size_t i){
        while(_parent[i] != i){
            _parent[i] = _parent[_parent[i]];
            i = _parent[i];
        }
        return i;
    }

    void unite(std::size_t a, std::size_t b){
        a = find(a);
        b = find(b);
        if(a == b) return;
        if(_size[a] < _size[b]) std::swap(a, b);
        _parent[b] = a;
        _size[a] += _size[b];
    }

private:
    std::vector<std::size_t> _parent;
    std::vector<std::size_t> _size;
};

// Circuit A enters a vertex over aIn and leaves it over aOut, circuit B over bIn and bOut. Walks the facets around the
// vertex from bIn towards bOut: A crosses B there when it leaves into that side but did not enter from it.
bool BurgersLoopBuilder::leavesOtherCircuit(InterfaceMesh::Edge* aIn, InterfaceMesh::Edge* aOut, InterfaceMesh::Edge* bIn, InterfaceMesh::Edge* bOut){
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
bool BurgersLoopBuilder::crossedAtVertex(InterfaceMesh::Edge* edge1, InterfaceMesh::Edge* edge2, const CandidateLinks& links){
    InterfaceMesh::Edge* sentinel = edge1->oppositeEdge();
    InterfaceMesh::Edge* edge = sentinel;
    do{
        InterfaceMesh::Edge* incoming = edge->prevFaceEdge();
        InterfaceMesh::Edge* incomingNext = links.onCircuit(incoming) ? links.next(incoming) : nullptr;
        const bool circuitPassesThrough = incomingNext && incomingNext->vertex1() == incoming->vertex2();
        if(circuitPassesThrough && edge2->oppositeEdge()
           && leavesOtherCircuit(edge2->oppositeEdge(), edge1->oppositeEdge(), incoming, incomingNext)){
            return true;
        }
        edge = incoming->oppositeEdge();
    }while(edge != sentinel);
    return false;
}

// Whether the candidate circuit crosses a circuit on the mesh, or itself, as if it were linked into the mesh.
bool BurgersLoopBuilder::intersectsOtherCircuits(std::span<InterfaceMesh::Edge* const> candidate){
    const CandidateLinks links{candidate};
    for(std::size_t i = 0; i < candidate.size(); ++i){
        InterfaceMesh::Edge* edge1 = candidate[i];
        InterfaceMesh::Edge* edge2 = candidate[(i + 1) % candidate.size()];
        if(edge1 != edge2->oppositeEdge() && crossedAtVertex(edge1, edge2, links)) return true;
    }
    return false;
}

// The serial pass creates the primary circuits in the order of the closing edges, and each circuit it creates and
// traces changes which edges are free for the searches after it. Most searches create nothing, though: their circuits
// enclose no dislocation, or cross one already traced. So the searches are first run in parallel, as a look-ahead on
// the mesh as it is before the pass, and the serial pass repeats only those that start a segment there, or that start
// close enough to a change the pass made (within the search depth) to read it. Every other search would find on the
// mesh what the look-ahead found, and change nothing, so the result is the serial one whatever the thread count.
void BurgersLoopBuilder::findPrimarySegments(int maxCircuitLength){
    const int searchDepth = (maxCircuitLength - 1) / 2;
    assert(searchDepth >= 1);

    const FrozenFreeEdges frozenEdges(mesh());
    const std::vector<InterfaceMesh::Edge*> closingEdges = findClosingEdges(frozenEdges, searchDepth);

    std::vector<std::uint8_t> startsSegment(closingEdges.size(), 0);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, closingEdges.size(), 16), [&](const tbb::blocked_range<size_t>& range){
        CircuitSearch& search = threadSearch();
        for(size_t i = range.begin(); i < range.end(); ++i){
            if(closingEdges[i]->circuit != nullptr) continue;
            search.run(frozenEdges, closingEdges[i]->vertex1(), searchDepth, [&](InterfaceMesh::Edge* edge, SearchNode* from, SearchNode* reached){
                const PrimaryVerdict verdict = judgePrimaryCircuit(edge, from, reached, search.candidate());
                if(verdict == PrimaryVerdict::StartsSegment) startsSegment[i] = 1;
                return verdict != PrimaryVerdict::KeepSearching;
            });
        }
    });

    ChangedRegion changed(mesh().vertexCount());
    std::vector<InterfaceMesh::Vertex*> changedVertices;
    CircuitSearch& search = threadSearch();
    for(size_t i = 0; i < closingEdges.size(); ++i){
        InterfaceMesh::Edge* closingEdge = closingEdges[i];
        if(closingEdge->circuit != nullptr) continue;
        InterfaceMesh::Vertex* start = closingEdge->vertex1();
        if(!startsSegment[i] && !changed.contains(start)) continue;

        changedVertices.clear();
        _changedVertices = &changedVertices;
        search.run(LiveFreeEdges{}, start, searchDepth, [&](InterfaceMesh::Edge* edge, SearchNode* from, SearchNode* reached){
            return createPrimaryCircuit(edge, from, reached, maxCircuitLength, search.candidate());
        });
        _changedVertices = nullptr;
        if(!changedVertices.empty()) changed.spread(changedVertices, searchDepth);
    }
}

// The circuit that edge closes, out of the search paths that reach its two ends: down the path to reached from where
// the two paths meet, back over edge, and up the path to from.
void BurgersLoopBuilder::collectPrimaryCircuit(InterfaceMesh::Edge* edge, const SearchNode* from, const SearchNode* reached, std::vector<InterfaceMesh::Edge*>& circuit){
    // The paths are at most the search depth long, so walking from's path for each vertex costs less than a set.
    const auto onFromPath = [from](const InterfaceMesh::Vertex* vertex){
        for(const SearchNode* a = from; a != nullptr; a = a->parent){
            if(a->vertex == vertex) return true;
        }
        return false;
    };

    circuit.clear();
    const SearchNode* meet = reached;
    for(; !onFromPath(meet->vertex); meet = meet->parent){
        assert(meet->viaEdge != nullptr);
        circuit.push_back(meet->viaEdge);
    }
    std::ranges::reverse(circuit);
    circuit.push_back(edge->oppositeEdge());
    for(const SearchNode* a = from; a->vertex != meet->vertex; a = a->parent){
        assert(a->viaEdge != nullptr);
        circuit.push_back(a->viaEdge->oppositeEdge());
    }
    assert(circuit.size() >= 3);
}

// Whether the circuit that edge closes encloses a dislocation and starts a segment, or the search goes on (zero Burgers
// vector, or edges that do not add up to zero), or stops there because the circuit crosses an existing one, whose
// dislocation is already being traced. Reads the mesh only.
BurgersLoopBuilder::PrimaryVerdict BurgersLoopBuilder::judgePrimaryCircuit(
    InterfaceMesh::Edge* edge,
    const SearchNode* from,
    const SearchNode* reached,
    std::vector<InterfaceMesh::Edge*>& circuit,
    Vector3* burgersVector
){
    if(edge->circuit != nullptr) return PrimaryVerdict::KeepSearching;
    assert(from != reached);
    collectPrimaryCircuit(edge, from, reached, circuit);

    Vector3 edgeSum = Vector3::Zero();
    Matrix3 frankRotation = Matrix3::Identity();
    Vector3 b = Vector3::Zero();
    for(InterfaceMesh::Edge* e : circuit){
        assert(e->circuit == nullptr && e->nextCircuitEdge == nullptr);
        edgeSum += e->physicalVector;
        b += frankRotation * e->clusterVector;
        if(!e->clusterTransition->isSelfTransition())
            frankRotation = frankRotation * e->clusterTransition->reverse->tm;
    }
    assert(frankRotation.equals(Matrix3::Identity(), CA_TRANSITION_MATRIX_EPSILON));

    if(intersectsOtherCircuits(circuit)) return PrimaryVerdict::Stops;
    if(b.isZero(CA_LATTICE_VECTOR_EPSILON) || !edgeSum.isZero(CA_ATOM_VECTOR_EPSILON)) return PrimaryVerdict::KeepSearching;
    if(burgersVector) *burgersVector = b;
    return PrimaryVerdict::StartsSegment;
}

// Starts a segment from the circuit that edge closes if it encloses a dislocation. Returns whether to stop searching
// from this start vertex: once a segment was started, or when the circuit crosses an existing one.
bool BurgersLoopBuilder::createPrimaryCircuit(InterfaceMesh::Edge* edge, SearchNode* from, SearchNode* reached, int maxCircuitLength, std::vector<InterfaceMesh::Edge*>& circuit){
	Vector3 b;
	const PrimaryVerdict verdict = judgePrimaryCircuit(edge, from, reached, circuit, &b);
	if(verdict != PrimaryVerdict::StartsSegment) return verdict == PrimaryVerdict::Stops;

	BurgersCircuit* forwardCircuit = _circuitPool.construct();
	forwardCircuit->firstEdge = circuit.front();
	forwardCircuit->lastEdge = circuit.back();
	forwardCircuit->edgeCount = circuit.size();
	for(std::size_t i = 0; i < circuit.size(); ++i){
		circuit[i]->circuit = forwardCircuit;
		circuit[i]->nextCircuitEdge = circuit[(i + 1) % circuit.size()];
		if(_changedVertices) _changedVertices->push_back(circuit[i]->vertex1());
	}
	assert(forwardCircuit->countEdges() == forwardCircuit->edgeCount);

	createAndTraceSegment(ClusterVector(b, forwardCircuit->firstEdge->clusterTransition->cluster1), forwardCircuit, maxCircuitLength);
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
	BurgersCircuit* backwardCircuit = _circuitPool.construct();

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

	for(auto candidate : _defectCells->queryOverlapping(capBounds)){
		const auto slot = claimSlot(candidate);
		if(!slot) continue;

		CoreCellClaim& claim = _coreCellClaims[*slot];
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

// The cell's index into _coreCellClaims; DefectCellIndex gave the cells it left out a negative one.
std::optional<std::size_t> BurgersLoopBuilder::claimSlot(DelaunayTessellation::CellHandle cell) const{
	const int slot = _mesh.tessellation().getUserField(cell);
	if(slot < 0 || static_cast<std::size_t>(slot) >= _coreCellClaims.size()) return std::nullopt;
	return static_cast<std::size_t>(slot);
}

std::vector<int> BurgersLoopBuilder::coreAtomDislocationIds(std::size_t particleCount) const{
	std::vector<int> ids(particleCount, -1);
	if(_coreCellClaims.empty()) return ids;

	const DelaunayTessellation& tessellation = _mesh.tessellation();

	for(auto cellHandle : tessellation.cells()){
		const auto slot = claimSlot(cellHandle);
		if(!slot) continue;

		const CoreCellClaim& claim = _coreCellClaims[*slot];
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
// circuits hold all the edges beyond it. Whether two nodes already share a ring is asked of disjoint sets that follow
// the rings, which are built once from the rings as they are; a node outside the dangling list, which a dangling
// circuit should never border, sends every later question back to the ring itself.
void BurgersLoopBuilder::linkAdjacentCircuits(){
	const std::vector<DislocationNode*>& nodes = danglingNodes();
	for(std::size_t i = 0; i < nodes.size(); ++i){
		nodes[i]->circuit->linkIndex = i;
	}
	const auto indexOf = [&](const DislocationNode* node) -> std::optional<std::size_t>{
		const std::size_t i = node->circuit->linkIndex;
		if(i < nodes.size() && nodes[i] == node) return i;
		return std::nullopt;
	};

	RingSets rings(nodes.size());
	bool ringsFollowed = true;
	std::vector<std::uint8_t> seen(nodes.size(), 0);
	for(std::size_t i = 0; i < nodes.size(); ++i){
		if(seen[i]) continue;
		seen[i] = 1;
		for(DislocationNode* member = nodes[i]->junctionRing; member != nodes[i]; member = member->junctionRing){
			const auto j = indexOf(member);
			if(!j){
				ringsFollowed = false;
				continue;
			}
			seen[*j] = 1;
			rings.unite(i, *j);
		}
	}

	for(std::size_t i = 0; i < nodes.size(); ++i){
		DislocationNode* node = nodes[i];
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
				const auto j = ringsFollowed ? indexOf(adjacentNode) : std::nullopt;
				if(!j) ringsFollowed = false;

				const bool sameRing = j ? rings.find(i) == rings.find(*j) : node->formsJunctionWith(adjacentNode);
				if(!sameRing){
					node->connectNodes(adjacentNode);
					if(j) rings.unite(i, *j);
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

	BurgersCircuit* forwardCircuit = _circuitPool.construct();
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
