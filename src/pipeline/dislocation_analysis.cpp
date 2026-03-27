#include <volt/pipeline/dislocation_analysis.h>

#include <volt/analysis/structure_analysis.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/reconstructed_structure.h>
#include <volt/pipeline/burgers_loop_builder.h>
#include <volt/pipeline/delaunay_tessellation.h>
#include <volt/pipeline/elastic_mapping.h>
#include <volt/pipeline/interface_mesh.h>
#include <volt/utilities/json_utils.h>
#include <volt/helpers/dxa_serialization.h>

#include <spdlog/spdlog.h>

#include <memory>
#include <stdexcept>
#include <utility>

namespace Volt{

using namespace Volt::Particles;

DislocationAnalysis::DislocationAnalysis()
    : _referenceStructureLabel(0),
      _maxTrialCircuitSize(14),
      _circuitStretchability(9),
      _lineSmoothingLevel(10),
      _linePointInterval(2.5){}

void DislocationAnalysis::compute(const LammpsParser::Frame& frame, const std::string& outputFile){
    spdlog::info("Processing frame {} with {} atoms", frame.timestep, frame.natoms);

    FrameAdapter::PreparedAnalysisInput prepared;
    std::string frameError;
    if(!FrameAdapter::prepareAnalysisInput(frame, prepared, &frameError)){
        throw std::runtime_error(frameError);
    }

    if(_clustersTablePath.empty() || _clusterTransitionsPath.empty()){
        throw std::runtime_error("OpenDXA requires --clusters-table and --clusters-transitions for reconstruct Cluster Graph");
    }

    std::shared_ptr<ParticleProperty> positions = std::move(prepared.positions);
    
    spdlog::info("Trying structural identification context reconstruction");
    ReconstructedStructureContext context(positions.get(), frame.simulationCell);

    auto structureAnalysis = std::make_unique<StructureAnalysis>(context);
    std::string reconstructionError;
    if(!ReconstructedStructureLoader::load(
        frame,
        {_clustersTablePath, _clusterTransitionsPath},
        *structureAnalysis,
        context,
        &reconstructionError
    )){
        throw std::runtime_error(reconstructionError);
    }

    spdlog::info("Delaunay tessellation");
    DelaunayTessellation tessellation;
    double ghostLayerSize = 3.5 * structureAnalysis->maximumNeighborDistance();
    tessellation.generateTessellation(
        context.simCell,
        context.positions->constDataPoint3(),
        context.atomCount(),
        ghostLayerSize,
        false,
        nullptr
    );

    spdlog::info("Elastic mapping");
    ElasticMapping elasticMap(*structureAnalysis, tessellation);
    elasticMap.generateTessellationEdges();
    elasticMap.assignVerticesToClusters();
    elasticMap.assignIdealVectorsToEdges(false, 4);
    elasticMap.shrinkVertexStorage();

    spdlog::info("Creating interface mesh");
    InterfaceMesh interfaceMesh(elasticMap);
    interfaceMesh.createMesh(structureAnalysis->maximumNeighborDistance());

    elasticMap.releaseCaches();
    tessellation.releaseMemory();

    spdlog::info("Burgers loops construction");
    BurgersLoopBuilder tracer(
        interfaceMesh,
        &structureAnalysis->clusterGraph(),
        static_cast<int>(_maxTrialCircuitSize),
        static_cast<int>(_circuitStretchability)
    );

    tracer.traceDislocationSegments();
    tracer.finishDislocationSegments(_referenceStructureLabel);

    DislocationNetwork& network = tracer.network();
    spdlog::info("Found {} dislocation segments", network.segments().size());
    network.smoothDislocationLines(_lineSmoothingLevel, _linePointInterval);

    if(!outputFile.empty()){
        spdlog::info("Writing defect mesh data");
        JsonUtils::writeJsonMsgpackToFile(
            DxaSerialization::buildDefectMeshJson(interfaceMesh, interfaceMesh.structureAnalysis(), true),
            outputFile + "_defect_mesh.msgpack",
            false
        );

        spdlog::info("Writing dislocations data");
        JsonUtils::writeJsonMsgpackToFile(
            DxaSerialization::buildDislocationsJson(&network, &frame.simulationCell),
            outputFile + "_dislocations.msgpack",
            false
        );

        spdlog::info("Writing mesh data");
        JsonUtils::writeJsonMsgpackToFile(
            DxaSerialization::buildMeshJson(
                interfaceMesh,
                interfaceMesh.structureAnalysis(),
                true,
                &interfaceMesh
            ),
            outputFile + "_interface_mesh.msgpack",
            false
        );
    }

    structureAnalysis.reset();
    positions.reset();
}

}
