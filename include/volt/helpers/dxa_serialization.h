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

struct BurgersFamily{
    std::string name;
    std::string label;
};
BurgersFamily classifyBurgersFamily(const Vector3& localBurgers, const std::string& crystalStructure);

void streamDislocationsToFile(
    const std::string& linesFilePath,
    const std::string& summaryFilePath,
    const DislocationNetwork* network,
    const SimulationCell* simulationCell = nullptr,
    const DislocationsExportOptions& options = {},
    const StructureAnalysis* structureAnalysis = nullptr
);

struct InterfaceMeshFlags{
    bool isCompletelyGood;
    bool isCompletelyBad;
};

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
