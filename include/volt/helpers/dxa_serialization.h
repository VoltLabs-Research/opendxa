#pragma once

#include <nlohmann/json.hpp>

#include <vector>

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

// Splits a dislocation line where it crosses a periodic face of the cell, shifting the pieces by whole cell vectors
// so that the line stays inside the cell along its periodic directions: where one piece leaves through a face, the
// next one starts on the opposite face.
std::vector<std::vector<Point3>> splitAtPeriodicBoundaries(const std::vector<Point3>& line, const SimulationCell& cell);

void streamDislocationsToFile(
    const std::string& linesFilePath,
    const std::string& summaryFilePath,
    const DislocationNetwork& network,
    const StructureAnalysis& structureAnalysis,
    const DislocationsExportOptions& options
);

void streamMeshToFile(const std::string& filePath, const InterfaceMeshTopology& mesh, const SimulationCell& cell);

// The mesh, plus whether every tetrahedron of the tessellation was elastically compatible or none was.
void streamInterfaceMeshToFile(const std::string& filePath, const InterfaceMesh& interfaceMesh);

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
