#include <volt/pipeline/dislocation_analysis.h>

#include <volt/analysis/structure_analysis.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/reconstructed_structure.h>
#include <volt/pipeline/burgers_loop_builder.h>
#include <volt/pipeline/delaunay_tessellation.h>
#include <volt/pipeline/elastic_mapping.h>
#include <volt/pipeline/interface_mesh.h>
#include <volt/helpers/dxa_serialization.h>
#include <volt/analysis/structure_identification_export.h>

#include <spdlog/spdlog.h>

#include <memory>
#include <stdexcept>
#include <utility>
#include <chrono>

namespace Volt{

using namespace Volt::Particles;

DislocationAnalysis::DislocationAnalysis()
    : _maxTrialCircuitSize(14),
      _circuitStretchability(9),
      _lineSmoothingLevel(1.0),
      _linePointInterval(2.5),
      _ghostLayerScale(3.5),
      _interfaceAlphaScale(5.0),
      _crystalPathSteps(4),
      _exportDefectMesh(true),
      _exportInterfaceMesh(false),
      _exportDelaunayTessellation(false),
      _exportStructureIdentification(false),
      _exportCoherentCrystallineRegions(false),
      _exportDislocations(true),
      _exportCircuitInformation(true),
      _exportDislocationNetworkStats(true),
      _exportJunctions(true),
      _clipPbcSegments(true),
      _coverDomainWithFiniteTets(false){}

void DislocationAnalysis::compute(const LammpsParser::Frame& frame, const std::string& outputFile){
    spdlog::info("Processing frame {} with {} atoms", frame.timestep, frame.natoms);
    auto tTotal = std::chrono::high_resolution_clock::now();
    auto tStep = tTotal;

    auto elapsed = [&](){
        auto now = std::chrono::high_resolution_clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - tStep).count();
        tStep = now;
        return ms;
    };

    FrameAdapter::PreparedAnalysisInput prepared;
    std::string frameError;
    if(!FrameAdapter::prepareAnalysisInput(frame, prepared, &frameError)){
        throw std::runtime_error(frameError);
    }

    if(_clustersTablePath.empty() || _clusterTransitionsPath.empty()){
        throw std::runtime_error("OpenDXA requires --clusters_table and --clusters_transitions for reconstruct Cluster Graph");
    }

    std::shared_ptr<ParticleProperty> positions = std::move(prepared.positions);

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
    spdlog::info("[{:>6}ms] Structure reconstruction", elapsed());

    DelaunayTessellation tessellation;
    double ghostLayerSize = _ghostLayerScale * structureAnalysis->maximumNeighborDistance();
    tessellation.generateTessellation(
        context.simCell,
        context.positions->constDataPoint3(),
        context.atomCount(),
        ghostLayerSize,
        _coverDomainWithFiniteTets,
        nullptr
    );
    spdlog::info("[{:>6}ms] Delaunay tessellation", elapsed());

    ElasticMapping elasticMap(*structureAnalysis, tessellation);
    elasticMap.generateTessellationEdges();
    elasticMap.assignVerticesToClusters();
    elasticMap.assignIdealVectorsToEdges(false, _crystalPathSteps);
    elasticMap.shrinkVertexStorage();
    spdlog::info("[{:>6}ms] Elastic mapping", elapsed());

    InterfaceMesh interfaceMesh(elasticMap);
    interfaceMesh.createMesh(structureAnalysis->maximumNeighborDistance(), _interfaceAlphaScale);
    elasticMap.releaseCaches();
    spdlog::info("[{:>6}ms] Interface mesh", elapsed());

    BurgersLoopBuilder tracer(
        interfaceMesh,
        &structureAnalysis->clusterGraph(),
        static_cast<int>(_maxTrialCircuitSize),
        static_cast<int>(_circuitStretchability)
    );
    tracer.traceDislocationSegments();
    tracer.finishDislocationSegments(_referenceTopologyName);
    spdlog::info("[{:>6}ms] Burgers loop tracing", elapsed());

    DislocationNetwork& network = tracer.network();
    spdlog::info("Found {} dislocation segments", network.segments().size());
    network.smoothDislocationLines(_lineSmoothingLevel, _linePointInterval);
    spdlog::info("[{:>6}ms] Line smoothing", elapsed());

    if(!outputFile.empty()){
        if(_exportDefectMesh){
            spdlog::info("Writing defect mesh data");
            DxaSerialization::streamDefectMeshToFile(
                outputFile + "_defect_mesh.parquet",
                interfaceMesh,
                interfaceMesh.structureAnalysis(),
                true
            );
        }

        if(_exportDislocations){
            spdlog::info("Writing dislocations data");
            const DxaSerialization::DislocationsExportOptions exportOptions{
                _clipPbcSegments,
                _exportCircuitInformation,
                _exportDislocationNetworkStats,
                _exportJunctions
            };
            DxaSerialization::streamDislocationsToFile(
                outputFile + "_dislocations.parquet",
                &network,
                &frame.simulationCell,
                exportOptions
            );
        }

        if(_exportInterfaceMesh){
            spdlog::info("Writing mesh data");
            DxaSerialization::streamDefectMeshToFile(
                outputFile + "_interface_mesh.parquet",
                interfaceMesh,
                interfaceMesh.structureAnalysis(),
                true
            );
        }

        if(_exportDelaunayTessellation){
            spdlog::info("Writing Delaunay tessellation data");
            DxaSerialization::streamDelaunayTessellationToFile(
                outputFile + "_delaunay_tessellation.parquet", tessellation
            );
        }

        if(_exportStructureIdentification){
            spdlog::info("Writing structure identification data");
            StructureIdentificationExport::streamStructureIdentificationToParquet(
                outputFile + "_atoms.parquet", frame, *structureAnalysis
            );
        }

        if(_exportCoherentCrystallineRegions){
            spdlog::info("Writing coherent crystalline region data");
            DxaSerialization::streamCoherentCrystallineRegionsToFile(
                outputFile + "_coherent_crystalline_regions.parquet", frame, *structureAnalysis
            );
        }
    }

    spdlog::info("[{:>6}ms] Export", elapsed());

    auto totalMs = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - tTotal).count();
    spdlog::info("[{:>6}ms] TOTAL", totalMs);

    tessellation.releaseMemory();
    structureAnalysis.reset();
    positions.reset();
}

}
