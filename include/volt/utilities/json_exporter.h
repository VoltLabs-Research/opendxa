#pragma once
#include <nlohmann/json.hpp>
#include <string>
#include <vector>
#include <limits>
#include <volt/structures/dislocation_network.h>
#include <volt/geometry/interface_mesh.h>
#include <volt/analysis/burgers_circuit.h>
#include <volt/core/lammps_parser.h>
#include <volt/math/lin_alg.h>
#include <volt/analysis/analysis_context.h>

namespace Volt{

using json = nlohmann::json;

class BurgersLoopBuilder;

class DXAJsonExporter{
public:
    explicit DXAJsonExporter() = default;

    void exportForStructureIdentification(
        const LammpsParser::Frame& frame,
        const StructureAnalysis& structureAnalysis,
        const std::string& outputFilename
    );

    json exportDislocationsToJson(const DislocationNetwork* network, const SimulationCell* simulationCell = nullptr);
    static inline uint32_t checked_u32_size(std::size_t n){
        if(n > static_cast<std::size_t>(std::numeric_limits<uint32_t>::max())){
            throw std::runtime_error("JSON container too large for msgpack u32 header.");
        }
        return static_cast<uint32_t>(n);
    }


    void writeDislocationsMsgpackToFile(
        const DislocationNetwork* network,
        const SimulationCell& simulationCell,
        const std::string& filePath
    );

    void exportPTMData(const AnalysisContext& context, const std::vector<int>& ids, const std::string& outputFilename);

    json getNetworkStatistics(const DislocationNetwork* network, double cellVolume);
    json getJunctionInformation(const DislocationNetwork* network);
    json getCircuitInformation(const DislocationNetwork* network);
    json getExtendedSimulationCellInfo(const SimulationCell& cell);
    
    void exportCoreAtoms(
        const LammpsParser::Frame& frame,
        const std::vector<uint8_t>& coreAtomFlags,
        const std::string& outputFilename
    );

    template <typename MeshType>
    json getMeshData(
        const MeshType& mesh,
        const StructureAnalysis& structureAnalysis,
        bool includeTopologyInfo,
        const InterfaceMesh* interfaceMeshForTopology
    );

    template <typename MeshType>
    void writeMeshMsgpackToFile(
        const MeshType& mesh,
        const StructureAnalysis& structureAnalysis,
        bool includeTopologyInfo,
        const InterfaceMesh* interfaceMeshForTopology,
        const std::string& filePath
    );

    void writeDefectMeshMsgpackToFile(
        const InterfaceMesh& interfaceMesh,
        const BurgersLoopBuilder& tracer,
        const StructureAnalysis& structureAnalysis,
        bool includeTopologyInfo,
        const std::string& filePath
    );

private:
    json vectorToJson(const Vector3& vector);
    json affineTransformationToJson(const AffineTransformation& transform);
    json simulationCellToJson(const SimulationCell& cell);
    
    std::string getBurgersVectorString(const Vector3& burgers);
    
    int countJunctions(const DislocationNetwork* network);
    int countDanglingSegments(const DislocationNetwork* network);
    double calculateAngle(const Vector3& a, const Vector3& b);

};

} 
