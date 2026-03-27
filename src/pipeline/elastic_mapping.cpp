#include <volt/core/volt.h>
#include <volt/pipeline/burgers_loop_builder.h>
#include <volt/helpers/crystal_path_finder.h>
#include <volt/pipeline/elastic_mapping.h>
#include <tbb/parallel_for.h>
#include <mutex>
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <queue>

namespace Volt{

static constexpr std::array<std::pair<int, int>, 6> tetraEdgeVertices{{
    {0, 1}, {0, 2}, {0, 3},
    {1, 2}, {1, 3}, {2, 3}
}};

// In order to measure how each tetrahedron in our Delaunay mesh connects
// two atoms across the grain boundary, we talk every tetrahedral cell and record each
// of its six edges exactly once. We skip any "ghost" cells that lie outside the
// real simulation box. For each real edge, we look up the two vertex IDs (v1, v2),
// skip degenerate or wrapped edges, and then build a TessellationEdge object if that
// connection has not already been recorded. Each new edge is inserted into two linked
// lists: one at its source vertex (edges leaving) and one at its destination vertex
// (edges arriving), so that we can later traverse all edges adjacent to any given vertex.
void ElasticMapping::generateTessellationEdges(){
    const auto &simCell = structureAnalysis().context().simCell;
    _edges.clear();
    for(auto& heads : _vertexEdges){
        heads = {-1, -1};
    }

    std::vector<std::vector<uint64_t>> sortedEdgeChunks;
    std::mutex chunksMutex;

    tbb::parallel_for(tbb::blocked_range<size_t>(0, tessellation().numberOfTetrahedra(), 32'768), 
        [&](const tbb::blocked_range<size_t>& r) {
        std::vector<uint64_t> edgeKeys;
        edgeKeys.reserve((r.size() * tetraEdgeVertices.size()) / 2);

        std::array<decltype(tessellation().cellVertex(0, 0)), 4> vertices{};
        std::array<int, 4> vertexIndices{};
        std::array<Point3, 4> vertexPositions{};
        for(size_t cellIdx = r.begin(); cellIdx != r.end(); ++cellIdx){
            if(tessellation().isGhostCell(cellIdx)) continue;

            for(int localIdx = 0; localIdx < 4; ++localIdx){
                vertices[localIdx] = tessellation().cellVertex(cellIdx, localIdx);
                vertexIndices[localIdx] = tessellation().vertexIndex(vertices[localIdx]);
                vertexPositions[localIdx] = tessellation().vertexPosition(vertices[localIdx]);
            }

            for(auto [vi, vj] : tetraEdgeVertices){
                int v1 = vertexIndices[vi];
                int v2 = vertexIndices[vj];

                if(v1 == v2) continue;

                if(simCell.isWrappedVector(vertexPositions[vi] - vertexPositions[vj])) continue;

                int minV = std::min(v1, v2);
                int maxV = std::max(v1, v2);
                // Pack (minV, maxV) into one sortable key.
                const uint64_t key = (uint64_t{static_cast<uint32_t>(minV)} << 32) |
                                     uint64_t{static_cast<uint32_t>(maxV)};
                edgeKeys.push_back(key);
            }
        }

        if(edgeKeys.empty()){
            return;
        }

        std::sort(edgeKeys.begin(), edgeKeys.end());
        edgeKeys.erase(std::unique(edgeKeys.begin(), edgeKeys.end()), edgeKeys.end());
        if(edgeKeys.empty()){
            return;
        }

        std::lock_guard<std::mutex> lock(chunksMutex);
        sortedEdgeChunks.emplace_back(std::move(edgeKeys));
    });

    if(sortedEdgeChunks.empty()){
        _edgeCount = 0;
        return;
    }

    size_t maxUniqueEdges = 0;
    for(const auto& chunk : sortedEdgeChunks){
        maxUniqueEdges += chunk.size();
    }
    _edges.reserve(maxUniqueEdges);

    struct HeapEntry{
        uint64_t key;
        size_t chunkIndex;
        size_t offset;
    };
    auto compareEntries = [](const HeapEntry& a, const HeapEntry& b){
        return a.key > b.key;
    };
    std::priority_queue<
        HeapEntry,
        std::vector<HeapEntry>,
        decltype(compareEntries)
    > minHeap(compareEntries);

    for(size_t i = 0; i < sortedEdgeChunks.size(); ++i){
        if(!sortedEdgeChunks[i].empty()){
            minHeap.push({sortedEdgeChunks[i][0], i, 0});
        }
    }

    _edgeCount = 0;
    bool hasPreviousKey = false;
    uint64_t previousKey = 0;
    while(!minHeap.empty()){
        HeapEntry entry = minHeap.top();
        minHeap.pop();

        if(!hasPreviousKey || entry.key != previousKey){
            hasPreviousKey = true;
            previousKey = entry.key;
            const int v1 = static_cast<int>(entry.key >> 32);
            const int v2 = static_cast<int>(entry.key & 0xFFFFFFFFu);
            const int edgeIndex = static_cast<int>(_edges.size());
            _edges.emplace_back(v1, v2);
            auto& edge = _edges.back();
            edge.nextLeavingEdge = _vertexEdges[v1].first;
            _vertexEdges[v1].first = edgeIndex;

            edge.nextArrivingEdge = _vertexEdges[v2].second;
            _vertexEdges[v2].second = edgeIndex;
            ++_edgeCount;
        }

        auto& keys = sortedEdgeChunks[entry.chunkIndex];
        const size_t nextOffset = entry.offset + 1;
        if(nextOffset < keys.size()){
            minHeap.push({keys[nextOffset], entry.chunkIndex, nextOffset});
        }else{
            std::vector<uint64_t>().swap(keys);
        }
    }

    std::vector<std::vector<uint64_t>>().swap(sortedEdgeChunks);

    _edges.shrink_to_fit();
}

// Once we have a graph of edges connecting mesh vertices, we need to assign
// each vertex to the grain (cluster) it belongs to.
// Initially, vertices that coincide exactly with an atomic cluster
// center get that cluster's ID; other vertices start with zero.
// In a simple propagation loop, we look at each unassigned vertex and check
// its neighboring vertices (both edges leaving and arriving).
// As soon as it touches a vertex already assigned to a nonzero cluster,
// we adopt that cluster ID. We repeat the scan until no changes occur,
// so that every vertex on the interface inherits the grain identity
// from at lesat one of its neighbors.
void ElasticMapping::assignVerticesToClusters(){
    const size_t vertex_count = _vertexClusters.size();
    
    // Initial assignment (can be parallel since each vertex is independent)
    tbb::parallel_for(tbb::blocked_range<size_t>(0, vertex_count),
        [this](const tbb::blocked_range<size_t>& r){
            for(size_t i = r.begin(); i < r.end(); ++i){
                _vertexClusters[i] = structureAnalysis().atomCluster(int(i));
            }
        });

    std::vector<Cluster*> nextClusters(vertex_count, nullptr);

    bool changed;
    do{
        std::atomic<bool> anyChanged{false};
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, vertex_count),
            [&](const tbb::blocked_range<size_t>& r){
                for(size_t idx = r.begin(); idx < r.end(); ++idx){
                    Cluster* currentCluster = _vertexClusters[idx];
                    Cluster* assignedCluster = currentCluster;
                    if(currentCluster->id == 0){
                        for(int edgeIdx = _vertexEdges[idx].first; edgeIdx >= 0; edgeIdx = _edges[edgeIdx].nextLeavingEdge){
                            auto const* e = &_edges[edgeIdx];
                            Cluster* neighborCluster = _vertexClusters[e->vertex2];
                            if(neighborCluster->id != 0){
                                assignedCluster = neighborCluster;
                                break;
                            }
                        }

                        if(assignedCluster->id == 0){
                            for(int edgeIdx = _vertexEdges[idx].second; edgeIdx >= 0; edgeIdx = _edges[edgeIdx].nextArrivingEdge){
                                auto const* e = &_edges[edgeIdx];
                                Cluster* neighborCluster = _vertexClusters[e->vertex1];
                                if(neighborCluster->id != 0){
                                    assignedCluster = neighborCluster;
                                    break;
                                }
                            }
                        }
                    }

                    nextClusters[idx] = assignedCluster;
                    if(assignedCluster != currentCluster){
                        anyChanged.store(true, std::memory_order_relaxed);
                    }
                }
            }
        );

        changed = anyChanged.load(std::memory_order_relaxed);
        if(changed){
            _vertexClusters.swap(nextClusters);
        }
    }while(changed);
}

// With every mesh edge now knowing the grain ID of its two endpoints, we want to compute
// an "ideal" Burgers vector on each edge so that when we later trace dislocation loops
// we know how the lattice would distort ideally between the two grains.
// We instantiate a helper class, CrystalPatFinder, which can find a lattice-aligned path
// between two atomic sites. For each edge that doesn't already have a vector, we check
// that both its vertices belong to valid clusters. We the ask the path finder for the 
// ideal vector and the grain where that vector originates. If necessary, we apply the inverse
// of the transition that brough that vector from its source cluster into the first cluster,
// then compose with the transition from the first cluster to the second. The final result
// is stored on the edge so that later elastic compatibility checks can 
// verify closed-loops balances.
void ElasticMapping::assignIdealVectorsToEdges(bool reconstructEdgeVectors, int crystalPathSteps){
    tbb::parallel_for(tbb::blocked_range<size_t>(0, _vertexEdges.size()), [&](const tbb::blocked_range<size_t>& r){
        (void)reconstructEdgeVectors;
        CrystalPathFinder pathFinder{ structureAnalysis(), crystalPathSteps };
        for(size_t headIdx = r.begin(); headIdx != r.end(); ++headIdx){
            for(int edgeIdx = _vertexEdges[headIdx].first; edgeIdx >= 0; edgeIdx = _edges[edgeIdx].nextLeavingEdge){
                auto* edge = &_edges[edgeIdx];
                if(edge->hasClusterVector()) { continue; }

                Cluster* c1 = clusterOfVertex(edge->vertex1);
                Cluster* c2 = clusterOfVertex(edge->vertex2);
                
                if(c1->id == 0 || c2->id == 0) continue;

                if(auto optCv = pathFinder.findPath(edge->vertex1, edge->vertex2)){
                    Vector3 localVec = optCv->localVec();
                    Cluster* srcCl = optCv->cluster();

                    Vector3 vecInC1;
                    if(srcCl == c1){
                        vecInC1 = localVec;
                    }else if(auto* tr = clusterGraph().determineClusterTransition(srcCl, c1)){
                        vecInC1 = tr->transform(localVec);
                    }else{
                        continue;
                    }

                    if(auto* tr12 = clusterGraph().determineClusterTransition(c1, c2)){
                        edge->assignClusterVector(vecInC1, tr12);
                    }
                }
            }
        }
    });
}

// Before accepting the elastic mapping as valid for simulation or further analysis
// we must confirm that every tetrahedron's six edges close consistently. That means each
// triangular face in the tetrahedron must satisfy both:
// 1) The sum of Burgers vectors around that triange is zero.
// 2) Any lattice symmetry transitions across those edges combine to the identity rotation.
// We extract each of the four unique 3-edge circuits on the tetrahedron, transform
// and sum their stored vectors, and test for zero magnitude.
// We likewise compose their transition matrices and confirm no net rotation.
// If any check fails, the mapping is compatible and we return false.
bool ElasticMapping::isElasticMappingCompatible(DelaunayTessellation::CellHandle cell) const{
    if(!tessellation().isValidCell(cell)) return false;

    // Gather each of the six edge's vector and transition
    std::array<std::pair<Vector3, ClusterTransition*>, 6> edgeVecs;
    for(int i = 0; i < 6; ++i){
        auto [vi, vj] = tetraEdgeVertices[i];
        int v1 = tessellation().vertexIndex(tessellation().cellVertex(cell, vi));
        int v2 = tessellation().vertexIndex(tessellation().cellVertex(cell, vj));
        auto* te = findEdge(v1, v2);

        // Every edge must exist and have a stored vector
        if(!te || !te->hasClusterVector()){
            return false;
		}

        // Make sure the vector is always oriented from v1 to v2
        if(te->vertex1 == v1){
            edgeVecs[i] = { te->clusterVector, te->clusterTransition };
        }else{
            // If reversed, invert the vector and use the reverse transition
            edgeVecs[i] = {
                te->clusterTransition->transform(-te->clusterVector),
                te->clusterTransition->reverse
            };
        }
    }

    // Define four triangular loops on a tetrahedron by edge indices
    static constexpr std::array<std::array<int, 3>, 4> circuits{{
        {{0, 4, 2}}, {{1, 5, 2}}, {{0, 3, 1}}, {{3, 5, 4}}
    }};

    // Check that the vector sum around each triangle is zero
    for(auto const& c : circuits){
        Vector3 B = edgeVecs[c[0]].first
                  + edgeVecs[c[0]].second->reverseTransform(edgeVecs[c[1]].first)
                  - edgeVecs[c[2]].first;

		if(!B.isZero(CA_LATTICE_VECTOR_EPSILON)) return false;
    }

    // Check that the combine rotations around each triangle are identity
    for(auto const& c : circuits){
        auto* t1 = edgeVecs[c[0]].second;
        auto* t2 = edgeVecs[c[1]].second;
        auto* t3 = edgeVecs[c[2]].second;
        if(!t1->isSelfTransition() ||
            !t2->isSelfTransition() ||
            !t3->isSelfTransition()
		){
            Matrix3 R = t3->reverse->tm * t2->tm * t1->tm;
		if(!R.equals(Matrix3::Identity(), CA_TRANSITION_MATRIX_EPSILON)) return false;
        }
    }

    return true;
}

void ElasticMapping::releaseCaches() noexcept{
    std::vector<TessellationEdge>().swap(_edges);
    _edgeCount = 0;
    std::vector<std::pair<int, int>>().swap(_vertexEdges);
    std::vector<Cluster*>().swap(_vertexClusters);
}

}
