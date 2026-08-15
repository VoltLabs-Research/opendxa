#include <volt/core/volt.h>
#include <volt/pipeline/burgers_loop_builder.h>
#include <volt/helpers/crystal_path_finder.h>
#include <volt/pipeline/elastic_mapping.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
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

    size_t totalKeys = 0;
    for(const auto& chunk : sortedEdgeChunks)
        totalKeys += chunk.size();

    std::vector<uint64_t> allKeys;
    allKeys.reserve(totalKeys);
    for(auto& chunk : sortedEdgeChunks){
        allKeys.insert(allKeys.end(), chunk.begin(), chunk.end());
        std::vector<uint64_t>().swap(chunk);
    }
    std::vector<std::vector<uint64_t>>().swap(sortedEdgeChunks);

    tbb::parallel_sort(allKeys.begin(), allKeys.end());
    allKeys.erase(std::unique(allKeys.begin(), allKeys.end()), allKeys.end());

    _edges.reserve(allKeys.size());
    _edgeCount = 0;
    for(uint64_t key : allKeys){
        const int v1 = static_cast<int>(key >> 32);
        const int v2 = static_cast<int>(key & 0xFFFFFFFFu);
        const int edgeIndex = static_cast<int>(_edges.size());
        _edges.emplace_back(v1, v2);
        auto& edge = _edges.back();
        edge.nextLeavingEdge = _vertexEdges[v1].first;
        _vertexEdges[v1].first = edgeIndex;
        edge.nextArrivingEdge = _vertexEdges[v2].second;
        _vertexEdges[v2].second = edgeIndex;
        ++_edgeCount;
    }
}

void ElasticMapping::assignVerticesToClusters(){
    const size_t vertex_count = _vertexClusters.size();
    
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

bool ElasticMapping::isElasticMappingCompatible(DelaunayTessellation::CellHandle cell) const{
    if(!tessellation().isValidCell(cell)) return false;

    std::array<std::pair<Vector3, ClusterTransition*>, 6> edgeVecs;
    for(int i = 0; i < 6; ++i){
        auto [vi, vj] = tetraEdgeVertices[i];
        int v1 = tessellation().vertexIndex(tessellation().cellVertex(cell, vi));
        int v2 = tessellation().vertexIndex(tessellation().cellVertex(cell, vj));
        auto* te = findEdge(v1, v2);

        if(!te || !te->hasClusterVector()){
            return false;
		}

        if(te->vertex1 == v1){
            edgeVecs[i] = { te->clusterVector, te->clusterTransition };
        }else{
            edgeVecs[i] = {
                te->clusterTransition->transform(-te->clusterVector),
                te->clusterTransition->reverse
            };
        }
    }

    static constexpr std::array<std::array<int, 3>, 4> circuits{{
        {{0, 4, 2}}, {{1, 5, 2}}, {{0, 3, 1}}, {{3, 5, 4}}
    }};

    for(auto const& c : circuits){
        Vector3 B = edgeVecs[c[0]].first
                  + edgeVecs[c[0]].second->reverseTransform(edgeVecs[c[1]].first)
                  - edgeVecs[c[2]].first;

		if(!B.isZero(CA_LATTICE_VECTOR_EPSILON)) return false;
    }

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
