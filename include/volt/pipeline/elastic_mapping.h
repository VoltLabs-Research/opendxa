#pragma once

#include <volt/core/volt.h>
#include <volt/pipeline/delaunay_tessellation.h>
#include <volt/structures/cluster.h>
#include <volt/structures/cluster_graph.h>
#include <volt/analysis/structure_analysis.h>

#include <optional>
#include <utility>

namespace Volt {

class ElasticMapping{
    struct TessellationEdge{
        int vertex1;
        int vertex2;
        int nextLeavingEdge = -1;
        int nextArrivingEdge = -1;
        ClusterTransition* clusterTransition = nullptr;
        Vector3 clusterVector{};

        TessellationEdge(int v1, int v2) noexcept : vertex1(v1), vertex2(v2){}

        [[nodiscard]] bool hasClusterVector() const noexcept{
			return clusterTransition != nullptr;
		}

        void assignClusterVector(Vector3 const& v, ClusterTransition* t) noexcept{
            clusterVector = v;
            clusterTransition = t;
        }
    };

public:
    explicit ElasticMapping(StructureAnalysis& sa, DelaunayTessellation& tess) noexcept
        : _structureAnalysis(sa)
        , _tessellation(tess)
        , _clusterGraph(sa.clusterGraph())
        , _vertexEdges(sa.context().atomCount(), {-1, -1})
        , _vertexClusters(sa.context().atomCount(), nullptr){}

    [[nodiscard]] auto structureAnalysis() const noexcept -> StructureAnalysis& {
		return _structureAnalysis;
	}

    [[nodiscard]] auto tessellation() noexcept -> DelaunayTessellation& {
		return _tessellation;
	}

    [[nodiscard]] auto tessellation() const noexcept -> DelaunayTessellation const& {
		return _tessellation;
	}

    [[nodiscard]] auto clusterGraph() noexcept -> ClusterGraph& {
		return _clusterGraph;
	}

    [[nodiscard]] auto clusterGraph() const noexcept -> ClusterGraph const& {
		return _clusterGraph;
	}

	void generateTessellationEdges();
    void assignVerticesToClusters();
    void assignIdealVectorsToEdges(int crystalPathSteps);
    [[nodiscard]] auto isElasticMappingCompatible(DelaunayTessellation::CellHandle cell) const -> bool;
    void releaseCaches() noexcept;

    [[nodiscard]] auto clusterOfVertex(int idx) const noexcept -> Cluster*{
		assert(idx < (int)_vertexClusters.size());
		return _vertexClusters[idx];
	}

    [[nodiscard]] auto getEdgeClusterVector(int v1, int v2) const -> std::pair<Vector3, ClusterTransition*>{
        const auto oriented = orientedClusterVector(v1, v2);
        assert(oriented);
        return *oriented;
    }

    [[nodiscard]] auto edgeCount() const noexcept -> int {
		return static_cast<int>(_edges.size());
	}

private:
    // The edge's ideal vector and transition read from v1 to v2, when the edge exists and has one.
    [[nodiscard]] auto orientedClusterVector(int v1, int v2) const -> std::optional<std::pair<Vector3, ClusterTransition*>>{
        const auto* e = findEdge(v1, v2);
        if(!e || !e->hasClusterVector()) return std::nullopt;
        if(e->vertex1 == v1) return std::pair{e->clusterVector, e->clusterTransition};
        return std::pair{e->clusterTransition->transform(-e->clusterVector), e->clusterTransition->reverse};
    }

    [[nodiscard]] auto firstClusteredNeighbor(int vertex) const noexcept -> Cluster*;

    [[nodiscard]] auto findEdge(int v1, int v2) const noexcept -> TessellationEdge const* {
        assert(v1 >= 0 && v1 < static_cast<int>(_vertexEdges.size()));
        assert(v2 >= 0 && v2 < static_cast<int>(_vertexEdges.size()));

		for(int edgeIdx = _vertexEdges[v1].first; edgeIdx >= 0; edgeIdx = _edges[edgeIdx].nextLeavingEdge){
            auto const* e = &_edges[edgeIdx];
            if(e->vertex2 == v2) return e;
		}

		for(int edgeIdx = _vertexEdges[v1].second; edgeIdx >= 0; edgeIdx = _edges[edgeIdx].nextArrivingEdge){
            auto const* e = &_edges[edgeIdx];
            if(e->vertex1 == v2) return e;
		}

        return nullptr;
    }

private:
    StructureAnalysis& _structureAnalysis;
    DelaunayTessellation& _tessellation;
    ClusterGraph& _clusterGraph;

    std::vector<TessellationEdge> _edges;
    std::vector<std::pair<int, int>> _vertexEdges;
    std::vector<Cluster*> _vertexClusters;
};

}
