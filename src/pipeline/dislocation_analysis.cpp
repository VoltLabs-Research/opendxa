#include <volt/pipeline/dislocation_analysis.h>

#include <volt/analysis/cluster_graph_io.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/structure_identification_export.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/progress.h>
#include <volt/core/reconstructed_structure.h>
#include <volt/helpers/dxa_serialization.h>
#include <volt/pipeline/burgers_loop_builder.h>
#include <volt/pipeline/delaunay_tessellation.h>
#include <volt/pipeline/elastic_mapping.h>
#include <volt/pipeline/interface_mesh.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/utilities/parquet_atom_writer.h>

#include <spdlog/spdlog.h>

#include <algorithm>
#include <chrono>
#include <functional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace Volt{

using namespace Volt::Particles;

struct DislocationOutputs{
    bool defectMesh = false;
    bool dislocations = false;
    bool interfaceMesh = false;
    bool tessellation = false;
    bool atoms = false;
    bool coreAtoms = false;
    bool coherentRegions = false;
    DxaSerialization::DislocationsExportOptions dislocationOptions;
};

namespace{

// Relative wall time of each phase on sample runs: how far each one moves the progress VOLT shows.
namespace Cost{
constexpr double Reconstruction = 2.0;
constexpr double Tessellation = 4.0;
constexpr double TessellationEdges = 1.0;
constexpr double ClusterLabels = 0.2;
constexpr double IdealVectors = 1.5;
constexpr double ElasticMapping = TessellationEdges + ClusterLabels + IdealVectors;
constexpr double InterfaceMesh = 2.5;
constexpr double BurgersCircuits = 6.0;
constexpr double CoreAtoms = 0.3;
constexpr double DefectMesh = 0.6;
constexpr double LineSmoothing = 0.1;

constexpr double DefectMeshFile = 1.0;
constexpr double DislocationsFile = 0.2;
constexpr double InterfaceMeshFile = 1.0;
constexpr double TessellationFile = 8.0;
constexpr double AtomsFile = 2.0;
constexpr double CoreAtomsFile = 0.2;
constexpr double CoherentRegionsFile = 3.0;
}

double costIf(bool runs, double cost){
    return runs ? cost : 0.0;
}

double outputsCost(const DislocationOutputs& outputs){
    return costIf(outputs.defectMesh, Cost::DefectMeshFile)
        + costIf(outputs.dislocations, Cost::DislocationsFile)
        + costIf(outputs.interfaceMesh, Cost::InterfaceMeshFile)
        + costIf(outputs.tessellation, Cost::TessellationFile)
        + costIf(outputs.atoms, Cost::AtomsFile)
        + costIf(outputs.coreAtoms, Cost::CoreAtomsFile)
        + costIf(outputs.coherentRegions, Cost::CoherentRegionsFile);
}

long long millisecondsSince(const Stopwatch& timer){
    return std::chrono::duration_cast<std::chrono::milliseconds>(timer.elapsedTime()).count();
}

long long lapMilliseconds(Stopwatch& timer){
    const long long milliseconds = millisecondsSince(timer);
    timer.restart();
    return milliseconds;
}

struct NetworkTopology{
    std::size_t closedLoops = 0;
    std::size_t infiniteLines = 0;
    std::size_t junctions = 0;
    std::size_t danglingEnds = 0;
};

NetworkTopology summarizeNetworkTopology(const DislocationNetwork& network){
    NetworkTopology topology;
    for(const DislocationSegment* segment : network.segments()){
        if(!segment) continue;
        if(segment->isClosedLoop()){
            if(!segment->line.empty() && segment->isInfiniteLine()) ++topology.infiniteLines;
            else ++topology.closedLoops;
            continue;
        }
        for(const DislocationNode* node : segment->nodes){
            if(node->isDangling()){
                ++topology.danglingEnds;
                continue;
            }
            // Every arm of a junction is visited; count the junction once, at its lowest-addressed arm.
            bool lowestArm = true;
            for(const DislocationNode* arm = node->junctionRing; arm != node; arm = arm->junctionRing){
                if(std::less<const DislocationNode*>{}(arm, node)){
                    lowestArm = false;
                    break;
                }
            }
            if(lowestArm) ++topology.junctions;
        }
    }
    return topology;
}

void reconstructStructure(
    const LammpsParser::Frame& frame,
    const std::string& neighborLatticePath,
    const ClusterGraphExportPaths& clusterPaths,
    StructureAnalysis& structure,
    ReconstructedStructureContext& context
){
    const Stopwatch timer;
    spdlog::info("Reconstructing structure from the upstream contract...");
    std::string error;
    if(!ReconstructedStructureLoader::load(frame, neighborLatticePath, clusterPaths, structure, context, &error)){
        throw std::runtime_error(error);
    }

    const ClusterGraph& clusterGraph = structure.clusterGraph();
    const auto clusterCount = std::ranges::count_if(
        clusterGraph.clusters(),
        [](const Cluster* cluster){ return cluster->id != 0; }
    );
    std::size_t clusteredAtoms = 0;
    if(context.atomClusters){
        const int* atomClusters = context.atomClusters->constDataInt();
        clusteredAtoms = static_cast<std::size_t>(std::count_if(
            atomClusters, atomClusters + context.atomClusters->size(),
            [](int clusterId){ return clusterId != 0; }
        ));
    }
    spdlog::info(
        "[{:>6}ms] Structure reconstruction: {} clusters, {} transitions, {} of {} atoms in clusters ({:.1f}%)",
        millisecondsSince(timer), clusterCount, clusterGraph.clusterTransitions().size(),
        clusteredAtoms, context.atomCount(), percentOf(clusteredAtoms, context.atomCount())
    );
}

void tessellate(
    DelaunayTessellation& tessellation,
    const StructureContext& context,
    double ghostLayerSize,
    bool coverDomainWithFiniteTets
){
    const Stopwatch timer;
    spdlog::info("Building Delaunay tessellation (ghost layer {:.2f})...", ghostLayerSize);
    tessellation.generateTessellation(
        context.simCell,
        context.positions->constDataPoint3(),
        context.atomCount(),
        ghostLayerSize,
        coverDomainWithFiniteTets,
        nullptr
    );
    spdlog::info(
        "[{:>6}ms] Delaunay tessellation: {} primary tetrahedra",
        millisecondsSince(timer), tessellation.numberOfPrimaryTetrahedra()
    );
}

void mapElastically(ElasticMapping& elasticMap, int crystalPathSteps){
    const Stopwatch total;
    Stopwatch lap;
    spdlog::info("Computing elastic mapping (crystal path search depth {})...", crystalPathSteps);
    ProgressSteps steps({Cost::TessellationEdges, Cost::ClusterLabels, Cost::IdealVectors});

    elasticMap.generateTessellationEdges();
    spdlog::info("[{:>6}ms]   ... tessellation edges: {}", lapMilliseconds(lap), elasticMap.edgeCount());

    steps.next("Assigning vertices to clusters");
    elasticMap.assignVerticesToClusters();
    spdlog::info("[{:>6}ms]   ... vertex-to-cluster labels", lapMilliseconds(lap));

    steps.next("Assigning ideal vectors to edges");
    elasticMap.assignIdealVectorsToEdges(false, crystalPathSteps);
    spdlog::info("[{:>6}ms]   ... ideal vectors on edges", lapMilliseconds(lap));

    elasticMap.shrinkVertexStorage();
    spdlog::info("[{:>6}ms] Elastic mapping", millisecondsSince(total));
}

void buildInterfaceMesh(InterfaceMesh& interfaceMesh, double alphaScale){
    const Stopwatch timer;
    const double maximumNeighborDistance = interfaceMesh.structureAnalysis().maximumNeighborDistance();
    spdlog::info("Building interface mesh (alpha {:.2f})...", alphaScale * maximumNeighborDistance);
    interfaceMesh.createMesh(maximumNeighborDistance, alphaScale);
    interfaceMesh.elasticMapping().releaseCaches();
    spdlog::info(
        "[{:>6}ms] Interface mesh: {} facets, {} vertices{}",
        millisecondsSince(timer), interfaceMesh.faceCount(), interfaceMesh.vertexCount(),
        interfaceMesh.isCompletelyGood() && !interfaceMesh.isCompletelyBad()
            ? " (no elastically incompatible tetrahedra)" : ""
    );
    if(interfaceMesh.isCompletelyBad()){
        spdlog::warn("No elastically compatible tetrahedra: there is no crystalline region to trace dislocations in");
    }
}

void traceDislocations(BurgersLoopBuilder& tracer, std::string_view referenceTopology){
    const Stopwatch timer;
    spdlog::info(
        "Tracing Burgers circuits (trial circuits up to {} edges, stretched up to {})...",
        tracer.maxTrialCircuitSize(), tracer.maxExtendedCircuitSize()
    );
    tracer.traceDislocationSegments();
    tracer.finishDislocationSegments(referenceTopology);
    const NetworkTopology topology = summarizeNetworkTopology(tracer.network());
    spdlog::info(
        "[{:>6}ms] Burgers loop tracing: {} segments ({} closed loops, {} infinite lines), {} junctions, {} dangling ends",
        millisecondsSince(timer), tracer.network().segments().size(), topology.closedLoops, topology.infiniteLines,
        topology.junctions, topology.danglingEnds
    );
}

std::vector<int> markCoreAtoms(const BurgersLoopBuilder& tracer, std::size_t atomCount){
    const Stopwatch timer;
    std::vector<int> dislocationIds = tracer.coreAtomDislocationIds(atomCount);
    const auto markedAtoms = std::ranges::count_if(dislocationIds, [](int id){ return id >= 0; });
    spdlog::info("[{:>6}ms] Core atom marking ({} atoms marked)", millisecondsSince(timer), markedAtoms);
    return dislocationIds;
}

void buildDefectMesh(const InterfaceMesh& interfaceMesh, const BurgersLoopBuilder& tracer, InterfaceMeshTopology& defectMesh){
    const Stopwatch timer;
    const bool watertight = interfaceMesh.generateDefectMesh(tracer, defectMesh);
    if(!watertight){
        spdlog::warn("Defect mesh is not watertight: at least one dislocation hole was left uncapped");
    }
    spdlog::info(
        "[{:>6}ms] Defect mesh ({} of {} interface facets kept, {} vertices{})",
        millisecondsSince(timer), defectMesh.faceCount(), interfaceMesh.faceCount(), defectMesh.vertexCount(),
        watertight ? ", watertight" : ""
    );
}

void smoothLines(DislocationNetwork& network, double smoothingLevel, double pointInterval, double cellVolume){
    const Stopwatch timer;
    network.smoothDislocationLines(smoothingLevel, pointInterval);
    double totalLineLength = 0.0;
    for(const DislocationSegment* segment : network.segments()){
        if(segment) totalLineLength += segment->calculateLength();
    }
    spdlog::info(
        "[{:>6}ms] Line smoothing: total line length {:.2f}, dislocation density {:.3e}",
        millisecondsSince(timer), totalLineLength, cellVolume > 0.0 ? totalLineLength / cellVolume : 0.0
    );
}

void writeCoreAtoms(
    const std::string& filePath,
    const LammpsParser::Frame& frame,
    const StructureAnalysis& structure,
    const std::vector<int>& dislocationIds
){
    std::vector<std::size_t> coreAtoms;
    for(std::size_t i = 0; i < dislocationIds.size(); ++i){
        if(dislocationIds[i] >= 0) coreAtoms.push_back(i);
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

    const StructureContext& context = structure.context();
    const int* structureTypes = context.structureTypes ? context.structureTypes->constDataInt() : nullptr;

    spdlog::info("Writing core atom data ({} atoms)", coreAtoms.size());
    streamAtomsToParquet(
        filePath,
        coreFrame,
        [&coreAtoms, structureTypes](std::size_t atomIndex){
            return std::string(structureTypeName(
                structureTypes
                    ? structureTypes[coreAtoms[atomIndex]]
                    : static_cast<int>(StructureType::OTHER)
            ));
        },
        [&coreAtoms, &dislocationIds](ColumnarAtomWriter& writer, std::size_t atomIndex){
            writer.field("dislocation_id", dislocationIds[coreAtoms[atomIndex]]);
        }
    );
}

struct AnalysisResults{
    const LammpsParser::Frame& frame;
    const StructureAnalysis& structure;
    const DelaunayTessellation& tessellation;
    const InterfaceMesh& interfaceMesh;
    InterfaceMeshTopology& defectMesh;
    const DislocationNetwork& network;
    const std::vector<int>& coreAtomDislocationIds;
};

void writeOutputs(const std::string& outputBase, const DislocationOutputs& outputs, const AnalysisResults& results){
    const Stopwatch timer;
    ProgressSteps steps({
        costIf(outputs.defectMesh, Cost::DefectMeshFile),
        costIf(outputs.dislocations, Cost::DislocationsFile),
        costIf(outputs.interfaceMesh, Cost::InterfaceMeshFile),
        costIf(outputs.tessellation, Cost::TessellationFile),
        costIf(outputs.atoms, Cost::AtomsFile),
        costIf(outputs.coreAtoms, Cost::CoreAtomsFile),
        costIf(outputs.coherentRegions, Cost::CoherentRegionsFile)
    });

    if(outputs.defectMesh){
        spdlog::info("Writing defect mesh data");
        DxaSerialization::streamMeshToFile(
            outputBase + "_defect_mesh.parquet", results.defectMesh, results.structure, true
        );
        results.defectMesh.clear();
    }
    steps.next();

    if(outputs.dislocations){
        spdlog::info("Writing dislocations data");
        DxaSerialization::streamDislocationsToFile(
            outputBase + "_dislocations.parquet",
            outputBase + "_dislocation_summary.parquet",
            results.network,
            results.frame.simulationCell,
            outputs.dislocationOptions,
            results.structure
        );
    }
    steps.next();

    if(outputs.interfaceMesh){
        spdlog::info("Writing interface mesh data");
        DxaSerialization::streamMeshToFile(
            outputBase + "_interface_mesh.parquet",
            results.interfaceMesh,
            results.structure,
            true,
            DxaSerialization::InterfaceMeshFlags{
                results.interfaceMesh.isCompletelyGood(),
                results.interfaceMesh.isCompletelyBad()
            }
        );
    }
    steps.next();

    if(outputs.tessellation){
        spdlog::info("Writing Delaunay tessellation data");
        DxaSerialization::streamDelaunayTessellationToFile(
            outputBase + "_delaunay_tessellation.parquet", results.tessellation
        );
    }
    steps.next();

    if(outputs.atoms){
        spdlog::info("Writing structure identification data");
        StructureIdentificationExport::streamStructureIdentificationToParquet(
            outputBase + "_atoms.parquet", results.frame, results.structure
        );
    }
    steps.next();

    if(outputs.coreAtoms){
        writeCoreAtoms(
            outputBase + "_core_atoms.parquet", results.frame, results.structure, results.coreAtomDislocationIds
        );
    }
    steps.next();

    if(outputs.coherentRegions){
        spdlog::info("Writing coherent crystalline region data");
        DxaSerialization::streamCoherentCrystallineRegionsToFile(
            outputBase + "_coherent_crystalline_regions.parquet", results.frame, results.structure
        );
    }

    spdlog::info("[{:>6}ms] Export", millisecondsSince(timer));
}

}

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

DislocationOutputs DislocationAnalysis::outputsFor(const std::string& outputBase) const{
    const bool writes = !outputBase.empty();
    return {
        .defectMesh = writes && _exportDefectMesh,
        .dislocations = writes && _exportDislocations,
        .interfaceMesh = writes && _exportInterfaceMesh,
        .tessellation = writes && _exportDelaunayTessellation,
        .atoms = writes && _exportStructureIdentification,
        .coreAtoms = writes && _markCoreAtoms,
        .coherentRegions = writes && _exportCoherentCrystallineRegions,
        .dislocationOptions = {
            _clipPbcSegments,
            _exportCircuitInformation,
            _exportDislocationNetworkStats,
            _exportJunctions
        },
    };
}

void DislocationAnalysis::compute(const LammpsParser::Frame& frame, const std::string& outputBase){
    if(_clustersTablePath.empty() || _clusterTransitionsPath.empty()){
        throw std::runtime_error("OpenDXA requires --clusters_table and --clusters_transitions for reconstruct Cluster Graph");
    }
    if(_neighborLatticePath.empty()){
        throw std::runtime_error("OpenDXA requires --neighbor_lattice (per-atom neighbor topology parquet) for reconstruct Cluster Graph");
    }

    spdlog::info("Processing frame {} with {} atoms", frame.timestep, frame.natoms);
    const Stopwatch timer;

    FrameAdapter::PreparedAnalysisInput input;
    std::string inputError;
    if(!FrameAdapter::prepareAnalysisInput(frame, input, &inputError)){
        throw std::runtime_error(inputError);
    }

    const DislocationOutputs outputs = outputsFor(outputBase);
    ProgressSteps steps({
        Cost::Reconstruction,
        Cost::Tessellation,
        Cost::ElasticMapping,
        Cost::InterfaceMesh,
        Cost::BurgersCircuits,
        costIf(_markCoreAtoms, Cost::CoreAtoms),
        costIf(outputs.defectMesh, Cost::DefectMesh),
        Cost::LineSmoothing,
        outputsCost(outputs)
    });

    JobProgress::instance().setText("Reconstructing structure from the upstream contract");
    ReconstructedStructureContext context(input.positions.get(), frame.simulationCell);
    StructureAnalysis structure(context);
    reconstructStructure(frame, _neighborLatticePath, {_clustersTablePath, _clusterTransitionsPath}, structure, context);

    steps.next("Building the Delaunay tessellation");
    DelaunayTessellation tessellation;
    tessellate(tessellation, context, _ghostLayerScale * structure.maximumNeighborDistance(), _coverDomainWithFiniteTets);

    steps.next("Computing the elastic mapping");
    ElasticMapping elasticMap(structure, tessellation);
    mapElastically(elasticMap, _crystalPathSteps);

    steps.next("Building the interface mesh");
    InterfaceMesh interfaceMesh(elasticMap);
    buildInterfaceMesh(interfaceMesh, _interfaceAlphaScale);

    steps.next("Tracing Burgers circuits");
    BurgersLoopBuilder tracer(
        interfaceMesh,
        static_cast<int>(_maxTrialCircuitSize),
        static_cast<int>(_circuitStretchability)
    );
    tracer.setMarkCoreAtoms(_markCoreAtoms);
    traceDislocations(tracer, _referenceTopologyName);

    steps.next("Marking core atoms");
    const std::vector<int> coreAtomDislocationIds = _markCoreAtoms
        ? markCoreAtoms(tracer, context.atomCount())
        : std::vector<int>{};

    steps.next("Building the defect mesh");
    InterfaceMeshTopology defectMesh;
    if(outputs.defectMesh){
        buildDefectMesh(interfaceMesh, tracer, defectMesh);
    }

    steps.next("Smoothing dislocation lines");
    smoothLines(tracer.network(), _lineSmoothingLevel, _linePointInterval, frame.simulationCell.volume3D());

    steps.next("Writing outputs");
    writeOutputs(outputBase, outputs, {
        frame, structure, tessellation, interfaceMesh, defectMesh, tracer.network(), coreAtomDislocationIds
    });

    spdlog::info("[{:>6}ms] TOTAL", millisecondsSince(timer));
}

}
