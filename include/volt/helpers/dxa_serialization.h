#pragma once

#include <nlohmann/json.hpp>

#include <optional>

#include <volt/helpers/dislocation_network.h>
#include <volt/pipeline/delaunay_tessellation.h>
#include <volt/pipeline/interface_mesh.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/core/lammps_parser.h>

namespace Volt::DxaSerialization{

using json = nlohmann::json;

struct DislocationsExportOptions{
    bool clipPbcSegments = true;
    bool exportCircuitInformation = true;
    bool exportDislocationNetworkStats = true;
    bool exportJunctions = true;
};

// Burgers-vector family classification in the cluster lattice frame. Public so
// downstream plugins (line-reconstruction-dxa) classify their own segments —
// family knowledge lives in the plugins, not in VOLT or the daemon.
struct BurgersFamily{
    std::string name;   // e.g. "1/2<110>"
    std::string label;  // e.g. "1/2<110> (Perfect)"
};
BurgersFamily classifyBurgersFamily(const Vector3& localBurgers, const std::string& crystalStructure);

// Streaming export — one row per dislocation segment in the standard VOLT line
// entity table (id, points, per-segment property columns), plus a separate
// JSON-payload summary file carrying network statistics and chart data.
void streamDislocationsToFile(
    const std::string& linesFilePath,
    const std::string& summaryFilePath,
    const DislocationNetwork* network,
    const SimulationCell* simulationCell = nullptr,
    const DislocationsExportOptions& options = {}
);

// Whole-mesh classification produced by InterfaceMesh::createMesh. It describes the interface
// mesh only, so the defect mesh — which is a filtered copy of it — omits the block instead of
// repeating a verdict that does not apply to it.
struct InterfaceMeshFlags{
    bool isCompletelyGood;
    bool isCompletelyBad;
};

// Serializes either DXA surface (interface mesh or defect mesh) to the standard mesh payload:
// export.MeshExporter for the 3D renderer, sub_listings for the daemon's tabular listing
// reader, plus an optional topology block. Takes the topology base type so the defect mesh,
// which is not an InterfaceMesh, goes through the same path.
void streamMeshToFile(
    const std::string& filePath,
    const InterfaceMeshTopology& mesh,
    const StructureAnalysis& structureAnalysis,
    bool includeTopologyInfo,
    const std::optional<InterfaceMeshFlags>& interfaceFlags = std::nullopt
);

void streamDelaunayTessellationToFile(
    const std::string& filePath,
    const DelaunayTessellation& tessellation
);

void streamCoherentCrystallineRegionsToFile(
    const std::string& filePath,
    const LammpsParser::Frame& frame,
    const StructureAnalysis& structureAnalysis
);

}
