#include <volt/pipeline/dislocation_analysis.h>

#include <volt/analysis/structure_analysis.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/reconstructed_structure.h>
#include <volt/analysis/cluster_graph_io.h>
#include <volt/pipeline/burgers_loop_builder.h>
#include <volt/pipeline/delaunay_tessellation.h>
#include <volt/pipeline/elastic_mapping.h>
#include <volt/pipeline/interface_mesh.h>
#include <volt/helpers/dxa_serialization.h>
#include <volt/analysis/structure_identification_export.h>
#include <volt/utilities/parquet_atom_writer.h>

#include <spdlog/spdlog.h>

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
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
      _markCoreAtoms(false),
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
    if(_neighborLatticePath.empty()){
        throw std::runtime_error("OpenDXA requires --neighbor_lattice (per-atom neighbor topology parquet) for reconstruct Cluster Graph");
    }

    std::shared_ptr<ParticleProperty> positions = std::move(prepared.positions);

    SimulationCell analysisCell = frame.simulationCell;
    ReconstructedStructureContext context(positions.get(), analysisCell);
    auto structureAnalysis = std::make_unique<StructureAnalysis>(context);

    const bool metricRescaleActive =
            std::abs(_metricRescaleX - 1.0) > 1e-12 ||
            std::abs(_metricRescaleY - 1.0) > 1e-12 ||
            std::abs(_metricRescaleZ - 1.0) > 1e-12;

    if(metricRescaleActive){
        spdlog::info(
            "Metric isotropization active: analysis frame rescaled by (1/{}, 1/{}, 1/{})",
            _metricRescaleX, _metricRescaleY, _metricRescaleZ
        );  

        Point3* pos = positions->dataPoint3();
        const std::size_t n = positions->size();
        for(std::size_t i = 0; i < n; i++){
            pos[i] = Point3(
                pos[i].x() / _metricRescaleX,
                pos[i].y() / _metricRescaleY,
                pos[i].z() / _metricRescaleZ
            );
        }

        AffineTransformation cellMatrix = frame.simulationCell.matrix();
        for(int col = 0; col < 4; col++){
            cellMatrix(0, col) /= _metricRescaleX;
            cellMatrix(1, col) /= _metricRescaleY;
            cellMatrix(2, col) /= _metricRescaleZ;
        }

        analysisCell.setMatrix(cellMatrix);
    }

    {
        std::string reconstructionError;
        if(!ReconstructedStructureLoader::load(
            frame,
            _neighborLatticePath,
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
    spdlog::info("[{:>6}ms]   ... tessellation edges", elapsed());
    elasticMap.assignVerticesToClusters();
    spdlog::info("[{:>6}ms]   ... vertex-to-cluster labels", elapsed());
    elasticMap.assignIdealVectorsToEdges(false, _crystalPathSteps);
    spdlog::info("[{:>6}ms]   ... ideal vectors on edges", elapsed());
    elasticMap.shrinkVertexStorage();
    spdlog::info("[{:>6}ms] Elastic mapping (total of the four above)", elapsed());

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
    tracer.setMarkCoreAtoms(_markCoreAtoms);
    tracer.traceDislocationSegments();
    tracer.finishDislocationSegments(_referenceTopologyName);
    spdlog::info("[{:>6}ms] Burgers loop tracing", elapsed());

    std::vector<int> coreAtomDislocationIds;
    if(_markCoreAtoms){
        coreAtomDislocationIds = tracer.coreAtomDislocationIds(structureAnalysis->context().atomCount());
        const auto markedAtoms = std::count_if(
            coreAtomDislocationIds.begin(), coreAtomDislocationIds.end(),
            [](int id){ return id >= 0; }
        );
        spdlog::info("[{:>6}ms] Core atom marking ({} atoms marked)", elapsed(), markedAtoms);
    }

    DislocationNetwork& network = tracer.network();
    spdlog::info("Found {} dislocation segments", network.segments().size());

    InterfaceMeshTopology defectMesh;
    const bool needDefectMesh = _exportDefectMesh && !outputFile.empty();
    if(needDefectMesh){
        if(!interfaceMesh.generateDefectMesh(tracer, defectMesh)){
            spdlog::warn("Defect mesh is not watertight: at least one dislocation hole was left uncapped");
        }
        spdlog::info(
            "[{:>6}ms] Defect mesh ({} of {} interface facets kept, {} vertices)",
            elapsed(), defectMesh.faceCount(), interfaceMesh.faceCount(), defectMesh.vertexCount()
        );
    }

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
        if(needDefectMesh){
            spdlog::info("Writing defect mesh data");
            DxaSerialization::streamMeshToFile(
                outputFile + "_defect_mesh.parquet",
                defectMesh,
                interfaceMesh.structureAnalysis(),
                true
            );

            defectMesh.clear();
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
                exportOptions,
                structureAnalysis.get()
            );
        }

        if(_exportInterfaceMesh){
            spdlog::info("Writing interface mesh data");
            DxaSerialization::streamMeshToFile(
                outputFile + "_interface_mesh.parquet",
                interfaceMesh,
                interfaceMesh.structureAnalysis(),
                true,
                DxaSerialization::InterfaceMeshFlags{
                    interfaceMesh.isCompletelyGood(),
                    interfaceMesh.isCompletelyBad()
                }
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

        if(!coreAtomDislocationIds.empty()){
            std::vector<std::size_t> coreAtoms;
            for(std::size_t i = 0; i < coreAtomDislocationIds.size(); ++i){
                if(coreAtomDislocationIds[i] >= 0) coreAtoms.push_back(i);
            }

            LammpsParser::Frame coreFrame;
            coreFrame.timestep = frame.timestep;
            coreFrame.natoms = static_cast<int>(coreAtoms.size());
            coreFrame.simulationCell = frame.simulationCell;
            coreFrame.positions.reserve(coreAtoms.size());
            coreFrame.ids.reserve(coreAtoms.size());
            for(const std::size_t source : coreAtoms){
                coreFrame.positions.push_back(source < frame.positions.size()
                    ? frame.positions[source] : Point3::Origin());
                coreFrame.ids.push_back(source < frame.ids.size()
                    ? frame.ids[source] : static_cast<int>(source));
            }

            spdlog::info("Writing core atom data ({} atoms)", coreAtoms.size());
            streamAtomsToParquet(
                outputFile + "_core_atoms.parquet",
                coreFrame,
                [](std::size_t){ return std::string("Core"); },
                [&coreAtoms, &coreAtomDislocationIds](ColumnarAtomWriter& writer, std::size_t atomIndex){
                    writer.field("dislocation_id", coreAtomDislocationIds[coreAtoms[atomIndex]]);
                }
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
