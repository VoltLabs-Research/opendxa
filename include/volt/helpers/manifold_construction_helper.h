#pragma once

#include <volt/core/volt.h>
#include <volt/core/simulation_cell.h>
#include <volt/core/particle_property.h>
#include <volt/pipeline/delaunay_tessellation.h>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <atomic>

#include <boost/functional/hash.hpp>
#include <type_traits>
#include <unordered_map>
#include <array>
#include <vector>
#include <cstdint>
#include <algorithm>

namespace Volt{

template<class HalfEdgeStructureType, bool FlipOrientation = false, bool CreateTwoSidedMesh = false>
class ManifoldConstructionHelper{
public:
	struct DefaultPrepareMeshFaceFunc{
		void operator()(
			typename HalfEdgeStructureType::Face*, 
			const std::array<int, 3>&, 
			const std::array<DelaunayTessellation::VertexHandle, 3>&, 
			DelaunayTessellation::CellHandle
		){}
	};

	struct DefaultLinkManifoldsFunc{
		void operator()(
			typename HalfEdgeStructureType::Edge*, 
			typename HalfEdgeStructureType::Edge*
		){}
	};

	ManifoldConstructionHelper(DelaunayTessellation& tessellation, HalfEdgeStructureType& outputMesh, double alpha, ParticleProperty* positions)
		: _tessellation(tessellation), _mesh(outputMesh), _alpha(alpha), _positions(positions){}

	template<typename CellRegionFunc, typename PrepareMeshFaceFunc = DefaultPrepareMeshFaceFunc, typename LinkManifoldsFunc = DefaultLinkManifoldsFunc>
	bool construct(
		CellRegionFunc&& determineCellRegion,
		PrepareMeshFaceFunc&& prepareMeshFaceFunc = PrepareMeshFaceFunc(),
		LinkManifoldsFunc&& linkManifoldsFunc = LinkManifoldsFunc()
	){
		if(!classifyTetrahedra(std::move(determineCellRegion))) return false;
		if(!createInterfaceFacets(std::move(prepareMeshFaceFunc))) return false;
		if(!linkHalfedges(std::move(linkManifoldsFunc))) return false;
		return true;
	}

	int spaceFillingRegion() const{
		return _spaceFillingRegion;
	}

private:
	template<typename CellRegionFunc>
	bool classifyTetrahedra(CellRegionFunc&& determineCellRegion){
		_numSolidCells = 0;
		_spaceFillingRegion = -2;

		std::vector<DelaunayTessellation::CellHandle> cells;
		for(auto cell : _tessellation.cells()){
			cells.push_back(cell);
		}

		std::vector<uint8_t> isFilled(cells.size(), 0);
		tbb::parallel_for(tbb::blocked_range<size_t>(0, cells.size(), 512),
			[&](const tbb::blocked_range<size_t>& r){
			for(size_t i = r.begin(); i < r.end(); ++i){
				auto cell = cells[i];
				bool filled = false;
				if(_tessellation.isValidCell(cell)){
					if(auto res = _tessellation.alphaTest(cell, _alpha)){
						filled = *res;
					}else{
						// sliver test
						int f = 0;
						for(; f < 4; ++f){
							auto nbr = _tessellation.mirrorFacet(cell, f).first;
							if(!_tessellation.isValidCell(nbr)) break;
							auto nr = _tessellation.alphaTest(nbr, _alpha);
							if(nr.has_value() && !nr.value()) break;
						}
						if(f == 4) filled = true;
					}
				}
				isFilled[i] = filled ? 1 : 0;
			}
		});

		// Parallel region assignment
		std::vector<int> regions(cells.size(), 0);
		tbb::parallel_for(tbb::blocked_range<size_t>(0, cells.size(), 4096),
			[&](const tbb::blocked_range<size_t>& r){
			for(size_t i = r.begin(); i < r.end(); ++i){
				if(isFilled[i])
					regions[i] = determineCellRegion(cells[i]);
			}
		});

		// Serial commit (cheap: only writes fields, no computation)
		for(size_t i = 0; i < cells.size(); ++i){
			auto cell = cells[i];
			_tessellation.setUserField(cell, regions[i]);

			if(!_tessellation.isGhostCell(cell)){
				if(_spaceFillingRegion == -2){
					_spaceFillingRegion = regions[i];
				}else if(_spaceFillingRegion != regions[i]){
					_spaceFillingRegion = -1;
				}
			}
		}

		int cellIndex = 0;
		for(size_t i = 0; i < cells.size(); ++i){
			auto cell = cells[i];
			if(regions[i] != 0 && !_tessellation.isGhostCell(cell)){
				_tessellation.setCellIndex(cell, cellIndex++);
				_numSolidCells++;
			}else{
				_tessellation.setCellIndex(cell, -1);
			}
		}

		if(_spaceFillingRegion == -2) _spaceFillingRegion = 0;
		
		return true;
	}

	template<typename PrepareMeshFaceFunc>
	bool createInterfaceFacets(PrepareMeshFaceFunc&& prepareMeshFaceFunc){
		std::vector<int> vertexMapIndices(_positions->size(), -1);

		_tetrahedraFaceList.clear();
        _tetrahedraFaceList.resize(_numSolidCells, { nullptr, nullptr, nullptr, nullptr });
		_faceLookupMap.clear();

		struct FacetCandidate{
			DelaunayTessellation::CellHandle cell;
			int internalIdx;
			int faceIdx;
			std::array<DelaunayTessellation::VertexHandle, 3> vertexHandles;
			std::array<int, 3> vertexIndices;
		};

		// Phase 1: Collect active cells (parallel filter)
		const size_t totalCells = _tessellation.numberOfTetrahedra();

		std::vector<uint8_t> isActive(totalCells, 0);
		tbb::parallel_for(tbb::blocked_range<size_t>(0, totalCells, 8192),
			[&](const tbb::blocked_range<size_t>& r){
			for(size_t cellIdx = r.begin(); cellIdx < r.end(); ++cellIdx){
				if(_tessellation.getCellIndex(static_cast<DelaunayTessellation::CellHandle>(cellIdx)) != -1)
					isActive[cellIdx] = 1;
			}
		});

		std::vector<uint32_t> activeCellIndices;
		activeCellIndices.reserve(totalCells / 8);
		for(size_t cellIdx = 0; cellIdx < totalCells; ++cellIdx){
			if(isActive[cellIdx])
				activeCellIndices.push_back(static_cast<uint32_t>(cellIdx));
		}

		// Phase 2: Find all interface facets (one big parallel pass)
		tbb::concurrent_vector<FacetCandidate> allCandidates;
		allCandidates.reserve(activeCellIndices.size());

		tbb::parallel_for(tbb::blocked_range<size_t>(0, activeCellIndices.size(), 512),
			[&](const tbb::blocked_range<size_t>& r){
			for(size_t i = r.begin(); i < r.end(); ++i){
				auto cell = static_cast<DelaunayTessellation::CellHandle>(activeCellIndices[i]);
				int solidRegion = _tessellation.getUserField(cell);
				int internalIdx = _tessellation.getCellIndex(cell);

				for(int f = 0; f < 4; ++f){
					auto mirrorFacet = _tessellation.mirrorFacet(cell, f);
					auto adjacentCell = mirrorFacet.first;
					if(_tessellation.getUserField(adjacentCell) == solidRegion) continue;

					FacetCandidate candidate;
					candidate.cell = cell;
					candidate.internalIdx = internalIdx;
					candidate.faceIdx = f;
					for(int v = 0; v < 3; ++v){
						candidate.vertexHandles[v] = _tessellation.cellVertex(
							cell,
							DelaunayTessellation::cellFacetVertexIndex(f, FlipOrientation ? (2 - v) : v)
						);
						candidate.vertexIndices[v] = _tessellation.vertexIndex(candidate.vertexHandles[v]);
					}
					allCandidates.push_back(candidate);
				}
			}
		});

		// Phase 3: Serial commit (only interface faces, not all 52M cells)
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
			if constexpr(!std::is_same_v<PrepareMeshFaceFunc, std::nullptr_t>){
				prepareMeshFaceFunc(face, candidate.vertexIndices, candidate.vertexHandles, candidate.cell);
			}

			auto orderedIndices = candidate.vertexIndices;
			reorderFaceVertices(orderedIndices);
			_faceLookupMap.insert({orderedIndices, face});
			_tetrahedraFaceList[static_cast<size_t>(candidate.internalIdx)][static_cast<size_t>(candidate.faceIdx)] = face;
		}

		return true;
	}

	typename HalfEdgeStructureType::Face* findAdjacentFace(DelaunayTessellation::CellHandle cell, int f, int e){
		int v1 = FlipOrientation ? DelaunayTessellation::cellFacetVertexIndex(f, 2-e) : DelaunayTessellation::cellFacetVertexIndex(f, (e+1)%3);
		int v2 = FlipOrientation ? DelaunayTessellation::cellFacetVertexIndex(f, (4-e)%3) : DelaunayTessellation::cellFacetVertexIndex(f, e);

		auto start = _tessellation.incident_facets(cell, v1, v2, cell, f);
		auto circ = start;
		--circ;

		int region = _tessellation.getUserField(cell);
		while(_tessellation.getUserField((*circ).first) == region){
			--circ;
		}

		auto mirror = _tessellation.mirrorFacet(*circ);
		return findCellFace(mirror);
	}

	template<typename LinkManifoldsFunc>
	bool linkHalfedges(LinkManifoldsFunc&& linkManifoldsFunc){
		std::vector<std::pair<DelaunayTessellation::CellHandle, size_t>> activeCells;
		activeCells.reserve(static_cast<size_t>(_numSolidCells));
		size_t tetIdx = 0;
		for(DelaunayTessellation::CellHandle cell : _tessellation.cells()){
			if(_tessellation.getCellIndex(cell) == -1) continue;
			activeCells.push_back({cell, tetIdx});
			++tetIdx;
		}

		tbb::parallel_for(tbb::blocked_range<size_t>(0, activeCells.size(), 256),
			[&](const tbb::blocked_range<size_t>& r){
			for(size_t ci = r.begin(); ci < r.end(); ++ci){
				auto [cell, ti] = activeCells[ci];
				for(int f = 0; f < 4; f++){
					auto* facet = _tetrahedraFaceList[ti][f];
					if(!facet) continue;

					auto* edge = facet->edges();
					for(int e = 0; e < 3; ++e, edge = edge->nextFaceEdge()){
						if(edge->oppositeEdge()) continue;
						auto* oppFace = findAdjacentFace(cell, f, e);
						if(oppFace){
							auto* oppEdge = oppFace->findEdge(edge->vertex2(), edge->vertex1());
							if(oppEdge) edge->linkToOppositeEdge(oppEdge);
						}
					}
				}
			}
		});

		if constexpr(CreateTwoSidedMesh){
			for(auto [cell, ti] : activeCells){
				for(int f = 0; f < 4; f++){
					auto* facet = _tetrahedraFaceList[ti][f];
					if(!facet) continue;

					auto oppFacet = _tessellation.mirrorFacet(cell, f);
					auto* outerFacet = findCellFace(oppFacet);
					if(!outerFacet) continue;

					auto* edge1 = facet->edges();
					for(int i = 0; i < 3; ++i, edge1 = edge1->nextFaceEdge()){
						for(auto* edge2 = outerFacet->edges(); ; edge2 = edge2->nextFaceEdge()) {
							if(edge2->vertex1() == edge1->vertex2()) {
								linkManifoldsFunc(edge1, edge2);
								break;
							}
						}
					}

					if(_tessellation.getUserField(oppFacet.first) == 0){
						auto* edge = outerFacet->edges();
						for(int e = 0; e < 3; ++e, edge = edge->nextFaceEdge()){
							if(edge->oppositeEdge()) continue;
							auto* oppFace = findAdjacentFace(oppFacet.first, oppFacet.second, e);
							if(!oppFace) continue;
							auto* oppEdge = oppFace->findEdge(edge->vertex2(), edge->vertex1());
							if(oppEdge) edge->linkToOppositeEdge(oppEdge);
						}
					}
				}
			}
		}
		return true;
	}

	typename HalfEdgeStructureType::Face* findCellFace(const std::pair<DelaunayTessellation::CellHandle,int>& facet){
		auto cell = facet.first;
		if(_tessellation.getCellIndex(cell) != -1){
			return _tetrahedraFaceList[_tessellation.getCellIndex(cell)][facet.second];
		}
		std::array<int,3> faceVerts;
		for(std::size_t i = 0; i < 3; ++i){
			int idx = DelaunayTessellation::cellFacetVertexIndex(facet.second, FlipOrientation ? (2-i) : i);
			faceVerts[i] = _tessellation.vertexIndex(_tessellation.cellVertex(cell, idx));
		}
		reorderFaceVertices(faceVerts);
		auto it = _faceLookupMap.find(faceVerts);
		return (it != _faceLookupMap.end()) ? it->second : nullptr;
	}

	static void reorderFaceVertices(std::array<int,3>& vertexIndices){
		std::rotate(
			vertexIndices.begin(), 
			std::min_element(vertexIndices.begin(), vertexIndices.end()), 
			vertexIndices.end()
		);
	}

	DelaunayTessellation& _tessellation;
	double _alpha;
	int _numSolidCells = 0;
	int _spaceFillingRegion = -1;
	ParticleProperty* _positions;
	HalfEdgeStructureType& _mesh;
	std::vector<std::array<typename HalfEdgeStructureType::Face*, 4>> _tetrahedraFaceList;
    std::unordered_map<std::array<int,3>, typename HalfEdgeStructureType::Face*, boost::hash<std::array<int, 3>>> _faceLookupMap;
};

}
