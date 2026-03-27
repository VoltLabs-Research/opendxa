#pragma once

#include <nlohmann/json.hpp>

#include <volt/helpers/dislocation_network.h>
#include <volt/pipeline/interface_mesh.h>
#include <volt/analysis/structure_analysis.h>

namespace Volt::DxaSerialization{

using json = nlohmann::json;

json buildDislocationsJson(
    const DislocationNetwork* network,
    const SimulationCell* simulationCell = nullptr
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

}
