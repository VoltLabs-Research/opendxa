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
    BurgersCircuit* circuit = nullptr;
};

struct InterfaceMeshEdge{
    Vector3 physicalVector;
    Vector3 clusterVector;
    ClusterTransition* clusterTransition = nullptr;
    BurgersCircuit* circuit = nullptr;
    HalfEdgeMesh<InterfaceMeshEdge, InterfaceMeshFace, InterfaceMeshVertex>::Edge* nextCircuitEdge = nullptr;
    
    std::atomic<bool> _claimedForCircuit{false};
    
    bool tryClaimForCircuit(){
        bool expected = false;
        return _claimedForCircuit.compare_exchange_strong(expected, true);
    }
    
    void releaseCircuitClaim(){
        _claimedForCircuit.store(false);
    }
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
