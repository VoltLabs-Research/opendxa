#include <volt/pipeline/dislocation_analysis.h>

#include <volt/analysis/structure_analysis.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/reconstructed_structure.h>
#include <volt/pipeline/burgers_loop_builder.h>
#include <volt/pipeline/delaunay_tessellation.h>
#include <volt/pipeline/elastic_mapping.h>
#include <volt/pipeline/interface_mesh.h>
#include <volt/helpers/dxa_serialization.h>
#include <volt/helpers/full_crystal_context.h>
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

    const bool fullCrystalMode = !_referenceCrystal.empty();
    if(!fullCrystalMode && (_clustersTablePath.empty() || _clusterTransitionsPath.empty())){
        throw std::runtime_error("OpenDXA requires --clusters_table and --clusters_transitions for reconstruct Cluster Graph");
    }

    std::shared_ptr<ParticleProperty> positions = std::move(prepared.positions);

    SimulationCell analysisCell = frame.simulationCell;
    ReconstructedStructureContext context(positions.get(), analysisCell);
    auto structureAnalysis = std::make_unique<StructureAnalysis>(context);
    
    const bool metricRescaleActive = 
            std::abs(_metricRescaleX - 1.0) > 1e-12 ||
            std::abs(_metricRescaleY - 1.0) > 1e-12 ||
            std::abs(_metricRescaleZ - 1.0) > 1e-12;

    // Full-crystal mode => build the elastic-mapping input from the perfect reference
    // (physical frame; the grain-frame ICP needs physical bonds), then rescale below.
    FullCrystalContextResult fullCrystalCtxResult;
    if(fullCrystalMode){
        FullCrystalContextParams params;
        params.referenceFile = _referenceCrystal;
        params.cationSpecies = _cationSpecies;
        params.bondCutoff = _fullCrystalCutoff;
        params.topologyName = _referenceTopologyName.empty() ? "crystal" : _referenceTopologyName;
        fullCrystalCtxResult = buildFullCrystalContext(frame, context, *structureAnalysis, params);

        if(!fullCrystalCtxResult.ok){
            throw std::runtime_error("Full-crystal context: " + fullCrystalCtxResult.message);
        }

        if(!metricRescaleActive){
            _metricRescaleX = fullCrystalCtxResult.metricRescaleX;
            _metricRescaleY = fullCrystalCtxResult.metricRescaleY;
            _metricRescaleZ = fullCrystalCtxResult.metricRescaleZ;
        }

        spdlog::info("[{:>6}ms] Full-crystal context (cutoff {:.2f} A, grain residual {:.3f} A)",
                elapsed(), fullCrystalCtxResult.selectedCutoff, fullCrystalCtxResult.grainSnapResidual);
    }

    if(metricRescaleActive){
        spdlog::info(
            "Metric isotropization active: analysis frame rescaled by (1/{}, 1/{}, 1/{})",
            _metricRescaleX, _metricRescaleY, _metricRescaleZ
        );  

        // Rescale atom positions in place.
        Point3* pos = positions->dataPoint3();
        const std::size_t n = positions->size();
        for(std::size_t i = 0; i < n; i++){
            pos[i] = Point3(
                pos[i].x() / _metricRescaleX,
                pos[i].y() / _metricRescaleY,
                pos[i].z() / _metricRescaleZ
            );
        }

        // Rescale the simulation cell (columns = cell vectors, column 3 = origin) so the
        // PBC wrapping in the isotropized frame stays consistent
        AffineTransformation cellMatrix = frame.simulationCell.matrix();
        for(int col = 0; col < 4; col++){
            cellMatrix(0, col) /= _metricRescaleX;
            cellMatrix(1, col) /= _metricRescaleY;
            cellMatrix(2, col) /= _metricRescaleZ;
        }

        analysisCell.setMatrix(cellMatrix);
    }

    if(!fullCrystalMode){
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
    }

    spdlog::info("[{:>6}ms] Structure reconstruction", elapsed());

    if(metricRescaleActive){
        // The ideal lattice vectors (neighbor_lattice_*) were produced in the physical
        // (ansiotropic) frame. The tetrahedron compatibility tests (isElasticMappingCompatible)
        // and the Burgers sum use an absolute tolerance (CA_LATTICE_VECTOR_EPSILON) tuned 
        // for ~unit-magnitude ideal vectors; ansiotropic vectors of magnitude ~a, b, c break those
        // checks so every good tetahedron is rejected.
        // Rescale the ideal vectors by (1/sx, 1/sy, 1/sz) into the same isotropized 
        // frame as the positions so the whole circuit machinery is self-consistent. 
        // The traced Burgers vector then comes out in isotropized lattice units; 
        // multiply its components by (sx, sy, sz) to recover physical Angstrom.
        if(structureAnalysis->hasNeighborLatticeVectorOverrides()){
            std::vector<Vector3> overrides = structureAnalysis->neighborLatticeVectorOverrides();
            const std::size_t stride= structureAnalysis->neighborLatticeVectorOverrideStride();
            for(Vector3& v : overrides){
                v = Vector3(
                    v.x() / _metricRescaleX,
                    v.y() / _metricRescaleY,
                    v.z() / _metricRescaleZ
                );
            }

            structureAnalysis->setNeighborLatticeVectorOverrides(std::move(overrides), stride);
        }

        // Full-crystal maxNeighborDistance was set in physical units; isotropize it.
        if(fullCrystalMode){
            const double smax = std::max({_metricRescaleX, _metricRescaleY, _metricRescaleZ});
            context.maximumNeighborDistance /= smax;
        }
    }

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

    if(metricRescaleActive){
        for(DislocationSegment* segment : network.segments()){
            if(!segment) continue;
            for(Point3& p : segment->line){
                p = Point3(
                    p.x() * _metricRescaleX,
                    p.y() * _metricRescaleY,
                    p.z() * _metricRescaleZ
                );
            }
        }
    }

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
                outputFile + "_dislocation_summary.parquet",
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
