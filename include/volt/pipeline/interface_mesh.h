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
    
    // Atomic flag for thread-safe edge claiming during parallel circuit creation
    std::atomic<bool> _claimedForCircuit{false};
    
    // Try to atomically claim this edge for circuit creation
    // Returns true if successfully claimed, false if already claimed
    bool tryClaimForCircuit(){
        bool expected = false;
        return _claimedForCircuit.compare_exchange_strong(expected, true);
    }
    
    void releaseCircuitClaim(){
        _claimedForCircuit.store(false);
    }
};

// The bare half-edge topology carrying the DXA per-element payloads. The interface mesh is
// one instance of it and the defect mesh derived by generateDefectMesh() is another, so code
// that accepts either one (exporters, in particular) takes this base type rather than
// InterfaceMesh: only the interface mesh owns an ElasticMapping and the good/bad flags.
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

    // Derives the defect mesh: this mesh minus every facet the Burgers-circuit tracing claimed
    // for a dislocation, with the holes left behind capped off. Two ordering constraints:
    //  - must run after BurgersLoopBuilder::traceDislocationSegments(), which is what fills
    //    BurgersCircuit::segmentMeshCap and marks the claimed facets;
    //  - must run before any coordinate transform of the dislocation lines (metric rescale,
    //    smoothing), because each cap vertex is placed at a dangling node's line endpoint and
    //    has to land in the same frame as the mesh vertices it gets stitched to.
    // `defectMesh` is expected to be empty. Returns true if the result came out closed; false
    // means some hole was left uncapped and the surface is not watertight.
    [[nodiscard]] bool generateDefectMesh(BurgersLoopBuilder const& tracer, InterfaceMeshTopology& defectMesh) const;

private:
    ElasticMapping& _elasticMapping;
    std::atomic<bool> _isCompletelyGood{true};
    std::atomic<bool> _isCompletelyBad{true};
};

}
