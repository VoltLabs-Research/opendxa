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
#include <algorithm>

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
    spdlog::debug("Processing frame {} with {} atoms", frame.timestep, frame.natoms);
    
    json result;

    if(_identificationMode == StructureAnalysis::Mode::CNA && _inputCrystalStructure == LATTICE_SC){
        return AnalysisResult::failure("CNA does not support SC. Use PTM for simple cubic crystals.");
    }

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
    
    {
        PROFILE("Identify Structures");
        if(structureAnalysis->usingPTM()){
            structureAnalysis->determineLocalStructuresWithPTM();
            structureAnalysis->computeMaximumNeighborDistanceFromPTM();
        }else{
            structureAnalysis->identifyStructuresCNA();
        }
    }

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

    // If identification mode is PTM, export PTM data
    if(!outputFile.empty() && _identificationMode == StructureAnalysis::Mode::PTM){
         _jsonExporter.exportPTMData(
            structureAnalysis->context(),
            frame.ids,
            outputFile
        );
    }

    // If structure identification only is requested
    if(_structureIdentificationOnly && !outputFile.empty()){
        _jsonExporter.exportForStructureIdentification(frame, *structureAnalysis, outputFile);

        result["is_failed"] = false;
        
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

    {
        PROFILE("Connect Clusters");
        clusterConnector.connectClusters();
    }

    {
        PROFILE("Form Super Clusters");
        clusterConnector.formSuperClusters();
    }

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

    ElasticMapping elasticMap(*structureAnalysis, tessellation);
    {
        PROFILE("Elastic Mapping - Generate Edges");
        elasticMap.generateTessellationEdges();
    }

    {
        PROFILE("Elastic Mapping - Assign Vertices");
        elasticMap.assignVerticesToClusters();
    }

    {
        PROFILE("Elastic Mapping - Assign Ideal Vectors");
        elasticMap.assignIdealVectorsToEdges(false, 4);
    }

    elasticMap.shrinkVertexStorage();
    
    structureAnalysis->freeNeighborLists();

    InterfaceMesh interfaceMesh(elasticMap);
    {
        PROFILE("InterfaceMesh - Create Mesh");
        interfaceMesh.createMesh(structureAnalysis->maximumNeighborDistance());
    }

    elasticMap.releaseCaches();
    tessellation.releaseMemory();

    BurgersLoopBuilder tracer(
        interfaceMesh, 
        &structureAnalysis->clusterGraph(),
        static_cast<int>(_maxTrialCircuitSize), 
        static_cast<int>(_circuitStretchability)
    );
    
    {
        PROFILE("Burgers Loop Builder - Trace Dislocation Segments");
        tracer.traceDislocationSegments();
    }

    {
        PROFILE("Burgers Loop Builder - Finish Dislocation Segments");
        tracer.finishDislocationSegments(_inputCrystalStructure);
    }

    if(!outputFile.empty()){
        PROFILE("Streaming Defect Mesh MsgPack");
        _jsonExporter.writeDefectMeshMsgpackToFile(
            interfaceMesh,
            tracer,
            interfaceMesh.structureAnalysis(),
            true,
            outputFile + "_defect_mesh.msgpack"
        );
    }

    DislocationNetwork& network = tracer.network();
    spdlog::debug("Found {} dislocation segments", network.segments().size());

    HalfEdgeMesh<InterfaceMeshEdge, InterfaceMeshFace, InterfaceMeshVertex> defectMesh;

    {
        PROFILE("Post Processing - Smooth Vertices & Smooth Dislocation Lines");
        network.smoothDislocationLines(_lineSmoothingLevel, _linePointInterval);
    }

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

        // atoms.msgpack: atoms grouped by structure type for AtomisticExporter
        {
            PROFILE("Streaming Atoms By Structure Type MsgPack");
            const auto& sa = interfaceMesh.structureAnalysis();
            constexpr int K = static_cast<int>(StructureType::NUM_STRUCTURE_TYPES);

            std::vector<std::string> names(K);
            for(int st = 0; st < K; st++){
                names[st] = sa.getStructureTypeName(st);
            }

            std::vector<std::vector<size_t>> structureAtomIndices(K);
            for(size_t i = 0; i < static_cast<size_t>(frame.natoms); ++i){
                const int raw = sa.context().structureTypes->getInt(static_cast<int>(i));
                const int st = (0 <= raw && raw < K) ? raw : 0;
                structureAtomIndices[static_cast<size_t>(st)].push_back(i);
            }

            std::vector<int> structureOrder;
            structureOrder.reserve(K);
            for(int st = 0; st < K; st++){
                if(!structureAtomIndices[static_cast<size_t>(st)].empty()){
                    structureOrder.push_back(st);
                }
            }
            std::sort(structureOrder.begin(), structureOrder.end(), [&](int a, int b){
                return names[a] < names[b];
            });

            json atomsByStructure;
            for(int st : structureOrder){
                json atomsArray = json::array();
                for(size_t atomIndex : structureAtomIndices[static_cast<size_t>(st)]){
                    const auto& pos = frame.positions[atomIndex];
                    atomsArray.push_back({
                        {"id", frame.ids[atomIndex]},
                        {"pos", { pos.x(), pos.y(), pos.z() }}
                    });
                }
                atomsByStructure[names[st]] = atomsArray;
            }

            json exportWrapper;
            exportWrapper["export"] = json::object();
            exportWrapper["export"]["AtomisticExporter"] = atomsByStructure;

            const std::string atomsPath = outputFile + "_atoms.msgpack";
            JsonUtils::writeJsonMsgpackToFile(exportWrapper, atomsPath, false);
        }
        {
            PROFILE("Streaming Atoms Around Dislocations MsgPack");
            const double surroundingAtomRadius = 1.5 * interfaceMesh.structureAnalysis().maximumNeighborDistance();
            _jsonExporter.exportAtomsAroundDislocations(
                frame,
                network,
                interfaceMesh.structureAnalysis(),
                surroundingAtomRadius,
                outputFile + "_dislocation_surrounding_atoms.msgpack"
            );
        }
        {
            PROFILE("Streaming Dislocations MsgPack");
            _jsonExporter.writeDislocationsMsgpackToFile(
                &network,
                frame.simulationCell,
                outputFile + "_dislocations.msgpack"
            );
        }

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
        
        {
            PROFILE("Streaming Simulation Cell MsgPack");
            auto simCellInfo = _jsonExporter.getExtendedSimulationCellInfo(frame.simulationCell);
            JsonUtils::writeJsonMsgpackToFile(simCellInfo, outputFile + "_simulation_cell.msgpack", false);
        }
    }
    
    structureAnalysis.reset();
    structureTypes.reset();
    positions.reset();

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    result["total_time"] = duration;

    spdlog::debug("Total time {} ms ", duration);

    return result;
}

}
