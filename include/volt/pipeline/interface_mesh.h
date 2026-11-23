#pragma once

#include <volt/core/volt.h>
#include <volt/helpers/half_edge_mesh.h>
#include <volt/pipeline/elastic_mapping.h>
#include <atomic>

namespace Volt{

struct BurgersCircuit;
struct BurgersCircuitSearchStruct;
class BurgersLoopBuilder;

struct InterfaceMeshVertex{
    std::atomic<BurgersCircuitSearchStruct*> burgersSearchStruct{nullptr};
};

struct InterfaceMeshFace{
    // Set on a facet a circuit swept before storing its cap. generateDefectMesh drops these, and keeps the facets a
    // still-dangling circuit swept after storing its cap: they lie between the cap and the circuit and close the mesh.
    static constexpr unsigned int kSweptBeforeCapFlag = 1;

    BurgersCircuit* circuit = nullptr;
};

struct InterfaceMeshEdge{
    Vector3 physicalVector;
    Vector3 clusterVector;
    ClusterTransition* clusterTransition = nullptr;
    BurgersCircuit* circuit = nullptr;
    HalfEdgeMesh<InterfaceMeshEdge, InterfaceMeshFace, InterfaceMeshVertex>::Edge* nextCircuitEdge = nullptr;
};

using InterfaceMeshTopology = HalfEdgeMesh<InterfaceMeshEdge, InterfaceMeshFace, InterfaceMeshVertex>;

class InterfaceMesh : public InterfaceMeshTopology{
public:
    explicit InterfaceMesh(ElasticMapping& mapping) noexcept
        : _elasticMapping(mapping){}

    [[nodiscard]] ElasticMapping& elasticMapping() noexcept{
		return _elasticMapping;
	}

    [[nodiscard]] ElasticMapping const& elasticMapping() const noexcept{
		return _elasticMapping;
	}

    [[nodiscard]] DelaunayTessellation& tessellation() noexcept{
		return elasticMapping().tessellation();
	}

    void makeManifold();

    [[nodiscard]] StructureAnalysis const& structureAnalysis() const noexcept{
        return elasticMapping().structureAnalysis();
    }

    void createMesh(double maximumNeighborDistance, double alphaScale = 5.0);

    [[nodiscard]] bool isCompletelyGood() const noexcept{
		return _isCompletelyGood.load();
	}

    [[nodiscard]] bool isCompletelyBad()  const noexcept{
		return _isCompletelyBad.load();
	}

    [[nodiscard]] bool generateDefectMesh(BurgersLoopBuilder const& tracer, InterfaceMeshTopology& defectMesh) const;

private:
    ElasticMapping& _elasticMapping;
    std::atomic<bool> _isCompletelyGood{true};
    std::atomic<bool> _isCompletelyBad{true};
};

}
