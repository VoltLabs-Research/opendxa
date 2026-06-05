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

json buildDislocationsJson(
    const DislocationNetwork* network,
    const SimulationCell* simulationCell = nullptr,
    const DislocationsExportOptions& options = {}
);

json buildMeshJson(
    const InterfaceMesh& mesh,
    const StructureAnalysis& structureAnalysis,
    bool includeTopologyInfo,
    const InterfaceMesh* interfaceMeshForTopology
);

json buildDefectMeshJson(
    const InterfaceMesh& interfaceMesh,
    const StructureAnalysis& structureAnalysis,
    bool includeTopologyInfo
);

json buildDelaunayTessellationJson(const DelaunayTessellation& tessellation);

json buildStructureIdentificationJson(
    const LammpsParser::Frame& frame,
    const StructureAnalysis& structureAnalysis
);

json buildCoherentCrystallineRegionsJson(
    const LammpsParser::Frame& frame,
    const StructureAnalysis& structureAnalysis
);

// Streaming export — writes directly to file without building DOM
void streamDislocationsToFile(
    const std::string& filePath,
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

}
