#pragma once

#include <nlohmann/json.hpp>

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

void streamDefectMeshToFile(
    const std::string& filePath,
    const InterfaceMesh& interfaceMesh,
    const StructureAnalysis& structureAnalysis,
    bool includeTopologyInfo
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
