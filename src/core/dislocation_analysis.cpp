#include <volt/core/dislocation_analysis.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/utilities/concurrence/parallel_system.h>
#include <volt/analysis/analysis_context.h>
#include <volt/analysis/cluster_connector.h>
#include <spdlog/spdlog.h>
#include <volt/utilities/json_utils.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <string_view>

namespace Volt{

using namespace Volt::Particles;

DislocationAnalysis::DislocationAnalysis()
    : _inputCrystalStructure(LATTICE_FCC),
      _maxTrialCircuitSize(14),
      _circuitStretchability(9),
      _lineSmoothingLevel(10),
      _linePointInterval(2.5),
      _defectMeshSmoothingLevel(8),
      _rmsd(0.12f),
      _identificationMode(StructureAnalysis::Mode::CNA),
      _markCoreAtoms(false),
      _structureIdentificationOnly(false),
      _onlyPerfectDislocations(false) {}

void DislocationAnalysis::setInputCrystalStructure(LatticeStructureType structure){
    _inputCrystalStructure = structure;
}

void DislocationAnalysis::setStructureIdentificationOnly(bool structureIdentificationOnly){
    _structureIdentificationOnly = structureIdentificationOnly;
}

void DislocationAnalysis::setMaxTrialCircuitSize(double size){
    _maxTrialCircuitSize= size;
}

void DislocationAnalysis::setCircuitStretchability(double stretch){
    _circuitStretchability = stretch;
}

void DislocationAnalysis::setMarkCoreAtoms(bool markCoreAtoms){
    _markCoreAtoms = markCoreAtoms;
}

void DislocationAnalysis::setOnlyPerfectDislocations(bool flag){
    _onlyPerfectDislocations = flag;
}

void DislocationAnalysis::setRmsd(float rmsd){
    _rmsd = rmsd;
}

void DislocationAnalysis::setLineSmoothingLevel(double lineSmoothingLevel){
    _lineSmoothingLevel = lineSmoothingLevel;
}

void DislocationAnalysis::setLinePointInterval(double linePointInterval){
    _linePointInterval = linePointInterval;
}

void DislocationAnalysis::setDefectMeshSmoothingLevel(double defectMeshSmoothingLevel){
    _defectMeshSmoothingLevel = defectMeshSmoothingLevel;
}

void DislocationAnalysis::setIdentificationMode(StructureAnalysis::Mode identificationMode){
    _identificationMode = identificationMode;
}

json DislocationAnalysis::compute(const LammpsParser::Frame &frame, const std::string& outputFile){
    auto start_time = std::chrono::high_resolution_clock::now();
    auto stage_start = start_time;
    spdlog::debug("Processing frame {} with {} atoms", frame.timestep, frame.natoms);
    
    json result;
    json stage_metrics = json::array();
    auto mark_stage = [&](std::string_view name){
        const auto now = std::chrono::high_resolution_clock::now();
        const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - stage_start).count();
        stage_metrics.push_back({
            {"stage", name},
            {"elapsed_ms", elapsed}
        });
        stage_start = now;
    };

    if(frame.natoms <= 0){
        return AnalysisResult::failure("Invalid number of atoms: " + std::to_string(frame.natoms));
    }

    if(frame.positions.empty()){
        return AnalysisResult::failure("No position data available");
    }

    std::shared_ptr<ParticleProperty> positions;
    {
        PROFILE("Create Position Property");
        positions = FrameAdapter::createPositionPropertyShared(frame);
        if(!positions){
            return AnalysisResult::failure("Failed to create position property");
        }
    }
    mark_stage("create_position_property");

    if(!FrameAdapter::validateSimulationCell(frame.simulationCell)){
        return AnalysisResult::failure("Invalid simulation cell");
    }

    // Default orientations (Identity)
    std::vector<Matrix3> preferredOrientations;
    preferredOrientations.push_back(Matrix3::Identity());

    auto structureTypes = std::make_unique<ParticleProperty>(frame.natoms, DataType::Int, 1, 0, true);
    AnalysisContext context(
        positions.get(),
        frame.simulationCell,
        _inputCrystalStructure,
        nullptr,
        structureTypes.get(),
        std::move(preferredOrientations)
    );

    std::unique_ptr<StructureAnalysis> structureAnalysis;
    {
        PROFILE("Structure Analysis Setup");
        structureAnalysis = std::make_unique<StructureAnalysis>(
            context,
            !_onlyPerfectDislocations,
            _identificationMode,
            _rmsd
        );
    }
    mark_stage("structure_analysis_setup");
    
    {
        PROFILE("Identify Structures");
        structureAnalysis->identifyStructures();
    }
    mark_stage("identify_structures");

    std::vector<int> extractedStructureTypes;
    if(outputFile.empty()){
        extractedStructureTypes.resize(frame.natoms);
        tbb::parallel_for(tbb::blocked_range<int>(0, frame.natoms),
            [&](const tbb::blocked_range<int>& r){
                for(int i = r.begin(); i < r.end(); ++i){
                    extractedStructureTypes[i] = structureAnalysis->context().structureTypes->getInt(i);
                }
            });
    }
    mark_stage("extract_structure_types");

    // If identification mode is PTM, export PTM data
    if(!outputFile.empty() && _identificationMode == StructureAnalysis::Mode::PTM){
         _jsonExporter.exportPTMData(
            structureAnalysis->context(),
            frame.ids,
            outputFile
        );
        mark_stage("export_ptm");
    }

    // If structure identification only is requested
    if(_structureIdentificationOnly && !outputFile.empty()){
        _jsonExporter.exportForStructureIdentification(frame, *structureAnalysis, outputFile);
        mark_stage("structure_identification_only_export");

        result["is_failed"] = false;
        result["stage_metrics"] = std::move(stage_metrics);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
        result["total_time"] = duration;
        
        return result;
    }

    // Standard Dislocation Analysis Pipeline
    ClusterConnector clusterConnector(*structureAnalysis, context);

    {
        PROFILE("Build Clusters");
        clusterConnector.buildClusters();
    }
    mark_stage("build_clusters");

    {
        PROFILE("Connect Clusters");
        clusterConnector.connectClusters();
    }
    mark_stage("connect_clusters");

    {
        PROFILE("Form Super Clusters");
        clusterConnector.formSuperClusters();
    }
    mark_stage("form_super_clusters");

    DelaunayTessellation tessellation;
    double ghostLayerSize;
    {
        PROFILE("Delaunay Tessellation");
        ghostLayerSize = 3.5f * structureAnalysis->maximumNeighborDistance();
        tessellation.generateTessellation(
            context.simCell,
            context.positions->constDataPoint3(),
            context.atomCount(),
            ghostLayerSize,
            false,
            nullptr
        );
    }
    mark_stage("delaunay_tessellation");

    ElasticMapping elasticMap(*structureAnalysis, tessellation);
    {
        PROFILE("Elastic Mapping - Generate Edges");
        elasticMap.generateTessellationEdges();
    }
    mark_stage("elastic_generate_edges");

    {
        PROFILE("Elastic Mapping - Assign Vertices");
        elasticMap.assignVerticesToClusters();
    }
    mark_stage("elastic_assign_vertices");

    {
        PROFILE("Elastic Mapping - Assign Ideal Vectors");
        elasticMap.assignIdealVectorsToEdges(false, 4);
    }
    mark_stage("elastic_assign_ideal_vectors");

    elasticMap.shrinkVertexStorage();
    
    structureAnalysis->freeNeighborLists();

    InterfaceMesh interfaceMesh(elasticMap);
    {
        PROFILE("InterfaceMesh - Create Mesh");
        interfaceMesh.createMesh(structureAnalysis->maximumNeighborDistance());
    }
    mark_stage("interface_mesh_create");

    elasticMap.releaseCaches();
    tessellation.releaseMemory();

    BurgersLoopBuilder tracer(
        interfaceMesh, 
        &structureAnalysis->clusterGraph(),
        static_cast<int>(_maxTrialCircuitSize), 
        static_cast<int>(_circuitStretchability),
        _markCoreAtoms
    );
    
    {
        PROFILE("Burgers Loop Builder - Trace Dislocation Segments");
        tracer.traceDislocationSegments();
    }
    mark_stage("trace_dislocation_segments");

    {
        PROFILE("Burgers Loop Builder - Finish Dislocation Segments");
        tracer.finishDislocationSegments(_inputCrystalStructure);
    }
    mark_stage("finish_dislocation_segments");

    if(!outputFile.empty()){
        PROFILE("Streaming Defect Mesh MsgPack");
        _jsonExporter.writeDefectMeshMsgpackToFile(
            interfaceMesh,
            tracer,
            interfaceMesh.structureAnalysis(),
            true,
            outputFile + "_defect_mesh.msgpack"
        );
        mark_stage("stream_defect_mesh_msgpack");
    }

    DislocationNetwork& network = tracer.network();
    spdlog::debug("Found {} dislocation segments", network.segments().size());

    HalfEdgeMesh<InterfaceMeshEdge, InterfaceMeshFace, InterfaceMeshVertex> defectMesh;

    {
        PROFILE("Post Processing - Smooth Vertices & Smooth Dislocation Lines");
        network.smoothDislocationLines(_lineSmoothingLevel, _linePointInterval);
    }
    mark_stage("smooth_dislocation_lines");

    double totalLineLength = 0.0;
    const auto& segments = network.segments();
    
    totalLineLength = tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, segments.size()),
        0.0,
        [&segments](const tbb::blocked_range<size_t>& r, double sum) -> double {
            for(size_t i = r.begin(); i < r.end(); ++i){
                DislocationSegment* segment = segments[i];
                if(segment && !segment->isDegenerate()){
                    sum += segment->calculateLength();
                }
            }
            return sum;
        },
        std::plus<double>()
    );

    spdlog::debug("Total line length: {} ", totalLineLength);

    result["is_failed"] = false;

    if(!outputFile.empty()){
        {
            PROFILE("Streaming Atoms MsgPack");
            _jsonExporter.exportForStructureIdentification(frame, interfaceMesh.structureAnalysis(), outputFile);
        }
        mark_stage("stream_atoms_msgpack");
        {
            PROFILE("Streaming Dislocations MsgPack");
            _jsonExporter.writeDislocationsMsgpackToFile(
                &network,
                frame.simulationCell,
                outputFile + "_dislocations.msgpack"
            );
        }
        mark_stage("stream_dislocations_msgpack");

        {
            PROFILE("Streaming Interface Mesh MsgPack");
            _jsonExporter.writeMeshMsgpackToFile(
                interfaceMesh,
                interfaceMesh.structureAnalysis(),
                true,
                &interfaceMesh,
                outputFile + "_interface_mesh.msgpack"
            );
        }
        mark_stage("stream_interface_mesh_msgpack");
        
        {
            PROFILE("Streaming Simulation Cell MsgPack");
            auto simCellInfo = _jsonExporter.getExtendedSimulationCellInfo(frame.simulationCell);
            JsonUtils::writeJsonMsgpackToFile(simCellInfo, outputFile + "_simulation_cell.msgpack", false);
        }
        mark_stage("stream_simulation_cell_msgpack");

        if(_markCoreAtoms){
            PROFILE("Streaming Core Atoms MsgPack");
            _jsonExporter.exportCoreAtoms(frame, tracer._coreAtomFlags, outputFile + "_core_atoms.msgpack");
            mark_stage("stream_core_atoms_msgpack");
        }
    }
    
    structureAnalysis.reset();
    structureTypes.reset();
    positions.reset();

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    result["total_time"] = duration;
    result["stage_metrics"] = std::move(stage_metrics);

    spdlog::debug("Total time {} ms ", duration);

    return result;
}

}
