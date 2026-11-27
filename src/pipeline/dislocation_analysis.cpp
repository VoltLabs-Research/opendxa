#include <volt/pipeline/dislocation_analysis.h>

#include <volt/analysis/cluster_graph_io.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/progress.h>
#include <volt/core/reconstructed_structure.h>
#include <volt/pipeline/burgers_loop_builder.h>
#include <volt/pipeline/delaunay_tessellation.h>
#include <volt/helpers/dislocation_analysis_log.h>
#include <volt/pipeline/elastic_mapping.h>
#include <volt/pipeline/interface_mesh.h>

#include <stdexcept>
#include <string>
#include <vector>

namespace Volt{

void DislocationAnalysis::reconstructStructure(
    const LammpsParser::Frame& frame,
    const std::string& neighborLatticePath,
    const ClusterGraphExportPaths& clusterPaths,
    StructureAnalysis& structure,
    ReconstructedStructureContext& context
){
    std::string error;
    if(!ReconstructedStructureLoader::load(frame, neighborLatticePath, clusterPaths, structure, context, &error)){
        throw std::runtime_error(error);
    }
}

void DislocationAnalysis::mapElastically(ElasticMapping& elasticMap, int crystalPathSteps){
    ProgressSteps steps({
        {"Building tessellation edges", 1.0},
        {"Assigning vertices to clusters", 0.2},
        {"Assigning ideal vectors to edges", 1.5},
    });
    elasticMap.generateTessellationEdges();

    steps.next();
    elasticMap.assignVerticesToClusters();

    steps.next();
    elasticMap.assignIdealVectorsToEdges(crystalPathSteps);
}

void DislocationAnalysis::compute(const LammpsParser::Frame& frame, const std::string& outputBase){
    if(_clustersTablePath.empty() || _clusterTransitionsPath.empty()){
        throw std::runtime_error("OpenDXA requires --clusters_table and --clusters_transitions for reconstruct Cluster Graph");
    }
    if(_neighborLatticePath.empty()){
        throw std::runtime_error("OpenDXA requires --neighbor_lattice (per-atom neighbor topology parquet) for reconstruct Cluster Graph");
    }

    DislocationAnalysisLog log(frame);
    FrameAdapter::PreparedAnalysisInput input;
    std::string inputError;
    if(!FrameAdapter::prepareAnalysisInput(frame, input, &inputError)){
        throw std::runtime_error(inputError);
    }

    ProgressSteps steps({
        {"Reconstructing structure from the upstream contract", 2.0},
        {"Building the Delaunay tessellation", 4.0},
        {"Computing the elastic mapping", 2.7},
        {"Building the interface mesh", 2.5},
        {"Tracing Burgers circuits", 6.0},
        {"Marking core atoms", 0.3, _exports.coreAtoms},
        {"Building the defect mesh", 0.6, _exports.defectMesh},
        {"Smoothing dislocation lines", 0.1},
        {"Writing outputs", DislocationOutputs::weight(_exports)},
    });

    ReconstructedStructureContext context(input.positions.get(), frame.simulationCell);
    StructureAnalysis structure(context);
    reconstructStructure(frame, _neighborLatticePath, {_clustersTablePath, _clusterTransitionsPath}, structure, context);
    log.dump(structure);

    steps.next();
    DelaunayTessellation tessellation;
    tessellation.generateTessellation(
        context.simCell, context.positions->constDataPoint3(), context.atomCount(),
        _ghostLayerScale * structure.maximumNeighborDistance(), _coverDomainWithFiniteTets, nullptr
    );
    log.dump(tessellation);

    steps.next();
    ElasticMapping elasticMap(structure, tessellation);
    mapElastically(elasticMap, _crystalPathSteps);
    log.dump(elasticMap);

    steps.next();
    InterfaceMesh interfaceMesh(elasticMap);
    interfaceMesh.createMesh(structure.maximumNeighborDistance(), _interfaceAlphaScale);
    elasticMap.releaseCaches();
    log.dump(interfaceMesh);

    steps.next();
    BurgersLoopBuilder tracer(
        interfaceMesh,
        static_cast<int>(_maxTrialCircuitSize),
        static_cast<int>(_circuitStretchability)
    );
    tracer.setMarkCoreAtoms(_exports.coreAtoms);
    tracer.traceDislocationSegments();
    tracer.finishDislocationSegments(_referenceTopologyName);
    log.dump(tracer);

    steps.next();
    std::vector<int> coreAtomDislocationIds;
    if(_exports.coreAtoms){
        coreAtomDislocationIds = tracer.coreAtomDislocationIds(context.atomCount());
        log.dump(coreAtomDislocationIds);
    }

    steps.next();
    InterfaceMeshTopology defectMesh;
    if(_exports.defectMesh){
        interfaceMesh.generateDefectMesh(tracer, defectMesh);
        log.dump(defectMesh);
    }

    steps.next();
    tracer.network().smoothDislocationLines(_lineSmoothingLevel, _linePointInterval);
    log.dump(tracer.network());

    steps.next();
    DislocationOutputs{
        outputBase, _exports, frame, structure, tessellation, interfaceMesh, defectMesh, tracer.network(),
        coreAtomDislocationIds
    }.write();
    log.dump(_exports);
}

}
