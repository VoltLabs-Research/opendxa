#pragma once

#include <volt/core/volt.h>
#include <volt/core/simulation_cell.h>
#include <volt/core/particle_property.h>
#include <volt/pipeline/delaunay_tessellation.h>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>

#include <boost/functional/hash.hpp>
#include <unordered_map>
#include <array>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <ranges>

namespace Volt{

// Builds the surface between the solid tetrahedra of an alpha shape and the rest of the tessellation: one face
// per facet that separates two regions, with each half-edge linked to its opposite across the surface.
template<class HalfEdgeStructureType>
class ManifoldConstructionHelper{
public:
    using CellHandle = DelaunayTessellation::CellHandle;
    using Face = typename HalfEdgeStructureType::Face;

    struct DefaultPrepareMeshFaceFunc{
        void operator()(
            Face*,
            const std::array<int, 3>&,
            const std::array<DelaunayTessellation::VertexHandle, 3>&,
            CellHandle
        ){}
    };

    ManifoldConstructionHelper(DelaunayTessellation& tessellation, HalfEdgeStructureType& outputMesh, double alpha, ParticleProperty* positions)
        : _tessellation(tessellation), _mesh(outputMesh), _alpha(alpha), _positions(positions){}

    template<typename CellRegionFunc, typename PrepareMeshFaceFunc = DefaultPrepareMeshFaceFunc>
    bool construct(
        CellRegionFunc&& determineCellRegion,
        PrepareMeshFaceFunc&& prepareMeshFaceFunc = PrepareMeshFaceFunc()
    ){
        classifyTetrahedra(std::move(determineCellRegion));
        createInterfaceFacets(std::move(prepareMeshFaceFunc));
        linkHalfedges();
        return true;
    }

private:
    // A degenerate tetrahedron has no circumsphere to test, so it takes the verdict of its neighbours:
    // it is filled unless one of them is infinite or fails the alpha test.
    bool isFilled(CellHandle cell) const{
        if(!_tessellation.isValidCell(cell)) return false;
        if(auto result = _tessellation.alphaTest(cell, _alpha)) return *result;
        return std::ranges::all_of(std::views::iota(0, 4), [&](int f){
            const CellHandle neighbor = _tessellation.mirrorFacet(cell, f).first;
            if(!_tessellation.isValidCell(neighbor)) return false;
            const auto neighborResult = _tessellation.alphaTest(neighbor, _alpha);
            return !neighborResult.has_value() || *neighborResult;
        });
    }

    // Stores each cell's region in its user field and numbers the primary solid cells, in cell order, through
    // the cell index; _solidCells maps that number back to the cell.
    template<typename CellRegionFunc>
    void classifyTetrahedra(CellRegionFunc&& determineCellRegion){
        const size_t cellCount = _tessellation.numberOfTetrahedra();

        std::vector<uint8_t> filled(cellCount, 0);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, cellCount, 512),
            [&](const tbb::blocked_range<size_t>& r){
            for(size_t i = r.begin(); i < r.end(); ++i){
                filled[i] = isFilled(static_cast<CellHandle>(i)) ? 1 : 0;
            }
        });

        std::vector<int> regions(cellCount, 0);
        tbb::parallel_for(tbb::blocked_range<size_t>(0, cellCount, 4096),
            [&](const tbb::blocked_range<size_t>& r){
            for(size_t i = r.begin(); i < r.end(); ++i){
                if(filled[i]) regions[i] = determineCellRegion(static_cast<CellHandle>(i));
            }
        });

        _solidCells.clear();
        for(size_t i = 0; i < cellCount; ++i){
            const auto cell = static_cast<CellHandle>(i);
            _tessellation.setUserField(cell, regions[i]);
            if(regions[i] != 0 && !_tessellation.isGhostCell(cell)){
                _tessellation.setCellIndex(cell, static_cast<int>(_solidCells.size()));
                _solidCells.push_back(cell);
            }else{
                _tessellation.setCellIndex(cell, -1);
            }
        }
    }

    template<typename PrepareMeshFaceFunc>
    void createInterfaceFacets(PrepareMeshFaceFunc&& prepareMeshFaceFunc){
        std::vector<int> vertexMapIndices(_positions->size(), -1);

        _tetrahedraFaceList.clear();
        _tetrahedraFaceList.resize(_solidCells.size(), { nullptr, nullptr, nullptr, nullptr });
        _faceLookupMap.clear();

        struct FacetCandidate{
            CellHandle cell;
            int internalIdx;
            int faceIdx;
            std::array<DelaunayTessellation::VertexHandle, 3> vertexHandles;
            std::array<int, 3> vertexIndices;
        };

        tbb::concurrent_vector<FacetCandidate> allCandidates;
        allCandidates.reserve(_solidCells.size());

        tbb::parallel_for(tbb::blocked_range<size_t>(0, _solidCells.size(), 512),
            [&](const tbb::blocked_range<size_t>& r){
            for(size_t i = r.begin(); i < r.end(); ++i){
                const CellHandle cell = _solidCells[i];
                const int solidRegion = _tessellation.getUserField(cell);

                for(int f = 0; f < 4; ++f){
                    const CellHandle adjacentCell = _tessellation.mirrorFacet(cell, f).first;
                    if(_tessellation.getUserField(adjacentCell) == solidRegion) continue;

                    FacetCandidate candidate;
                    candidate.cell = cell;
                    candidate.internalIdx = static_cast<int>(i);
                    candidate.faceIdx = f;
                    for(int v = 0; v < 3; ++v){
                        candidate.vertexHandles[v] = _tessellation.cellVertex(cell, DelaunayTessellation::cellFacetVertexIndex(f, v));
                        candidate.vertexIndices[v] = _tessellation.vertexIndex(candidate.vertexHandles[v]);
                    }
                    allCandidates.push_back(candidate);
                }
            }
        });

        // Faces are created serially in candidate order: the mesh's vertex and face numbering follows it.
        _faceLookupMap.reserve(allCandidates.size());
        for(const auto& candidate : allCandidates){
            std::array<typename HalfEdgeStructureType::Vertex*, 3> facetVertices{};
            for(int v = 0; v < 3; ++v){
                const int idx = candidate.vertexIndices[v];
                if(vertexMapIndices[idx] < 0){
                    vertexMapIndices[idx] = _mesh.createVertex(_positions->getPoint3(idx))->index();
                }
                facetVertices[v] = _mesh.vertex(vertexMapIndices[idx]);
            }

            auto* face = _mesh.createFace(facetVertices.begin(), facetVertices.end());
            prepareMeshFaceFunc(face, candidate.vertexIndices, candidate.vertexHandles, candidate.cell);

            auto orderedIndices = candidate.vertexIndices;
            reorderFaceVertices(orderedIndices);
            _faceLookupMap.insert({orderedIndices, face});
            _tetrahedraFaceList[static_cast<size_t>(candidate.internalIdx)][static_cast<size_t>(candidate.faceIdx)] = face;
        }
    }

    // Walks around the tetrahedron edge e of facet f until it leaves the cell's region; the facet found there,
    // seen from the other side, is the face that shares this edge.
    Face* findAdjacentFace(CellHandle cell, int f, int e){
        const int v1 = DelaunayTessellation::cellFacetVertexIndex(f, (e + 1) % 3);
        const int v2 = DelaunayTessellation::cellFacetVertexIndex(f, e);

        auto circ = _tessellation.incident_facets(cell, v1, v2, cell, f);
        --circ;

        const int region = _tessellation.getUserField(cell);
        while(_tessellation.getUserField((*circ).first) == region){
            --circ;
        }

        return findCellFace(_tessellation.mirrorFacet(*circ));
    }

    void linkHalfedges(){
        tbb::parallel_for(tbb::blocked_range<size_t>(0, _solidCells.size(), 256),
            [&](const tbb::blocked_range<size_t>& r){
            for(size_t solidIndex = r.begin(); solidIndex < r.end(); ++solidIndex){
                const CellHandle cell = _solidCells[solidIndex];
                for(int f = 0; f < 4; f++){
                    auto* facet = _tetrahedraFaceList[solidIndex][f];
                    if(!facet) continue;

                    auto* edge = facet->edges();
                    for(int e = 0; e < 3; ++e, edge = edge->nextFaceEdge()){
                        if(edge->oppositeEdge()) continue;
                        auto* oppFace = findAdjacentFace(cell, f, e);
                        if(!oppFace) continue;
                        auto* oppEdge = oppFace->findEdge(edge->vertex2(), edge->vertex1());
                        if(oppEdge) edge->linkToOppositeEdge(oppEdge);
                    }
                }
            }
        });
    }

    Face* findCellFace(const DelaunayTessellation::Facet& facet){
        const CellHandle cell = facet.first;
        if(_tessellation.getCellIndex(cell) != -1){
            return _tetrahedraFaceList[_tessellation.getCellIndex(cell)][facet.second];
        }
        std::array<int, 3> faceVerts;
        for(int i = 0; i < 3; ++i){
            const int idx = DelaunayTessellation::cellFacetVertexIndex(facet.second, i);
            faceVerts[i] = _tessellation.vertexIndex(_tessellation.cellVertex(cell, idx));
        }
        reorderFaceVertices(faceVerts);
        auto it = _faceLookupMap.find(faceVerts);
        return (it != _faceLookupMap.end()) ? it->second : nullptr;
    }

    // Rotates the smallest index to the front, so that a face has the same key whichever vertex it starts at.
    static void reorderFaceVertices(std::array<int, 3>& vertexIndices){
        std::rotate(
            vertexIndices.begin(),
            std::min_element(vertexIndices.begin(), vertexIndices.end()),
            vertexIndices.end()
        );
    }

    DelaunayTessellation& _tessellation;
    HalfEdgeStructureType& _mesh;
    double _alpha;
    ParticleProperty* _positions;
    std::vector<CellHandle> _solidCells;
    std::vector<std::array<Face*, 4>> _tetrahedraFaceList;
    std::unordered_map<std::array<int, 3>, Face*, boost::hash<std::array<int, 3>>> _faceLookupMap;
};

}
