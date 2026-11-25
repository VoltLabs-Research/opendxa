#include <volt/pipeline/dislocation_analysis_log.h>

#include <volt/analysis/structure_analysis.h>
#include <volt/helpers/dislocation_network.h>
#include <volt/pipeline/burgers_loop_builder.h>
#include <volt/pipeline/delaunay_tessellation.h>
#include <volt/pipeline/dislocation_outputs.h>
#include <volt/pipeline/elastic_mapping.h>

#include <spdlog/spdlog.h>

#include <algorithm>

namespace Volt{

DislocationAnalysisLog::NetworkTopology DislocationAnalysisLog::summarizeNetworkTopology(const DislocationNetwork& network){
    NetworkTopology topology;
    topology.junctions = network.junctionArmCounts().size();
    for(const DislocationSegment* segment : network.segments()){
        if(!segment) continue;
        if(segment->isClosedLoop()){
            if(!segment->line.empty() && segment->isInfiniteLine()) ++topology.infiniteLines;
            else ++topology.closedLoops;
            continue;
        }
        for(const DislocationNode* node : segment->nodes){
            if(node->isDangling()) ++topology.danglingEnds;
        }
    }
    return topology;
}

std::size_t DislocationAnalysisLog::clusteredAtomCount(const StructureContext& context){
    if(!context.atomClusters) return 0;
    const int* atomClusters = context.atomClusters->constDataInt();
    return static_cast<std::size_t>(std::count_if(
        atomClusters, atomClusters + context.atomClusters->size(),
        [](int clusterId){ return clusterId != 0; }
    ));
}

DislocationAnalysisLog::DislocationAnalysisLog(const LammpsParser::Frame& frame)
    : _cellVolume(frame.simulationCell.volume3D()){
    spdlog::info("Processing frame {} with {} atoms", frame.timestep, frame.natoms);
}

void DislocationAnalysisLog::dump(const StructureAnalysis& structure){
    const StructureContext& context = structure.context();
    const ClusterGraph& clusterGraph = structure.clusterGraph();
    const auto clusterCount = std::ranges::count_if(
        clusterGraph.clusters(),
        [](const Cluster* cluster){ return cluster->id != 0; }
    );
    const std::size_t clusteredAtoms = clusteredAtomCount(context);
    spdlog::info(
        "[{:>6}ms] Structure reconstruction: {} clusters, {} transitions, {} of {} atoms in clusters ({:.1f}%)",
        phaseMilliseconds(), clusterCount, clusterGraph.clusterTransitions().size(),
        clusteredAtoms, context.atomCount(), percentOf(clusteredAtoms, context.atomCount())
    );
}

void DislocationAnalysisLog::dump(const DelaunayTessellation& tessellation){
    spdlog::info(
        "[{:>6}ms] Delaunay tessellation: {} primary tetrahedra",
        phaseMilliseconds(), tessellation.numberOfPrimaryTetrahedra()
    );
}

void DislocationAnalysisLog::dump(const ElasticMapping& elasticMap){
    spdlog::info("[{:>6}ms] Elastic mapping: {} tessellation edges", phaseMilliseconds(), elasticMap.edgeCount());
}

void DislocationAnalysisLog::dump(const InterfaceMesh& interfaceMesh){
    _interfaceFacets = interfaceMesh.faceCount();
    spdlog::info(
        "[{:>6}ms] Interface mesh: {} facets, {} vertices{}",
        phaseMilliseconds(), interfaceMesh.faceCount(), interfaceMesh.vertexCount(),
        interfaceMesh.isCompletelyGood() && !interfaceMesh.isCompletelyBad()
            ? " (no elastically incompatible tetrahedra)" : ""
    );
    if(interfaceMesh.isCompletelyBad()){
        spdlog::warn("No elastically compatible tetrahedra: there is no crystalline region to trace dislocations in");
    }
}

void DislocationAnalysisLog::dump(const BurgersLoopBuilder& tracer){
    const DislocationNetwork& network = tracer.network();
    const NetworkTopology topology = summarizeNetworkTopology(network);
    spdlog::info(
        "[{:>6}ms] Burgers loop tracing: {} segments ({} closed loops, {} infinite lines), {} junctions, {} dangling ends",
        phaseMilliseconds(), network.segments().size(), topology.closedLoops, topology.infiniteLines,
        topology.junctions, topology.danglingEnds
    );
}

void DislocationAnalysisLog::dump(const std::vector<int>& coreAtomDislocationIds){
    const auto markedAtoms = std::ranges::count_if(coreAtomDislocationIds, [](int id){ return id >= 0; });
    spdlog::info("[{:>6}ms] Core atom marking ({} atoms marked)", phaseMilliseconds(), markedAtoms);
}

void DislocationAnalysisLog::dump(const InterfaceMeshTopology& defectMesh){
    const bool watertight = defectMesh.isClosed();
    if(!watertight){
        spdlog::warn("Defect mesh is not watertight: at least one dislocation hole was left uncapped");
    }
    spdlog::info(
        "[{:>6}ms] Defect mesh ({} of {} interface facets kept, {} vertices{})",
        phaseMilliseconds(), defectMesh.faceCount(), _interfaceFacets, defectMesh.vertexCount(),
        watertight ? ", watertight" : ""
    );
}

void DislocationAnalysisLog::dump(const DislocationNetwork& smoothedNetwork){
    double totalLineLength = 0.0;
    for(const DislocationSegment* segment : smoothedNetwork.segments()){
        if(segment) totalLineLength += segment->calculateLength();
    }
    spdlog::info(
        "[{:>6}ms] Line smoothing: total line length {:.2f}, dislocation density {:.3e}",
        phaseMilliseconds(), totalLineLength, _cellVolume > 0.0 ? totalLineLength / _cellVolume : 0.0
    );
}

void DislocationAnalysisLog::dump(const DislocationExports&){
    spdlog::info("[{:>6}ms] Export", phaseMilliseconds());
    spdlog::info("[{:>6}ms] TOTAL", runMilliseconds());
}

}
