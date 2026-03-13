#include <volt/utilities/json_exporter.h>
#include <volt/analysis/burgers_circuit.h>
#include <volt/analysis/burgers_loop_builder.h>
#include <iomanip>
#include <numeric>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <array>
#include <cstdint>
#include <volt/utilities/json_utils.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

namespace Volt {

void clipDislocationLine(
    const std::vector<Point3>& line,
    const SimulationCell& simulationCell,
    const std::function<void(const Point3&, const Point3&, bool)>& segmentCallback
){
    if(line.size() < 2) return;
    bool isInitialSegment = true;

    // initialize the first point and the shift vector
    auto v1Iter = line.cbegin();
    Point3 rp1 = simulationCell.absoluteToReduced(*v1Iter);
    Vector3 shiftVector = Vector3::Zero();
    for(size_t dimension = 0; dimension < 3; dimension++){
        if(simulationCell.pbcFlags()[dimension]){
            // move the start point to the main box [0,1) and record the offset
            double shift = -std::floor(rp1[dimension]);
            rp1[dimension] += shift;
            shiftVector[dimension] += shift;
        }
    }

    // iterate over the original line segments
    for(auto v2Iter = v1Iter + 1; v2Iter != line.cend(); v1Iter = v2Iter, ++v2Iter){
        Point3 rp2 = simulationCell.absoluteToReduced(*v2Iter) + shiftVector;
        // ugly hack
        int maxIterations = 10;
        int iterationCount = 0;
        do{
            iterationCount++;
            if(iterationCount > maxIterations){
                segmentCallback(
                    simulationCell.reducedToAbsolute(rp1),
                    simulationCell.reducedToAbsolute(rp2),
                    isInitialSegment
                );
                break;
            }

            size_t crossDim = -1;
            double crossDir = 0;
            double smallestT = std::numeric_limits<double>::max();
            for(size_t dimension = 0; dimension < 3; dimension++){
                if(simulationCell.pbcFlags()[dimension]){
                    // crossing detection
                    int d = (int) std::floor(rp2[dimension]) - (int) std::floor(rp1[dimension]);
                    if(d == 0) continue;

                    double dr = rp2[dimension] - rp1[dimension];
                    if(std::abs(dr) < 1e-9) continue;

                    double t = (d > 0) ? (std::ceil(rp1[dimension]) - rp1[dimension]) / dr
                                       : (std::floor(rp1[dimension]) - rp1[dimension]) / dr;
                    if(t > 1e-9 && t < smallestT){
                        smallestT = t;
                        crossDim = dimension;
                        crossDir = (d > 0) ? 1.0 : -1.0;
                    }
                }
            }

            // tolerance to avoid very small intersections
            if(smallestT < (1.0 - 1e-9)){
                Point3 intersection = rp1 + smallestT * (rp2 - rp1);
                intersection[crossDim] = std::round(intersection[crossDim]);
                segmentCallback(simulationCell.reducedToAbsolute(rp1), simulationCell.reducedToAbsolute(intersection), isInitialSegment);
                shiftVector[crossDim] -= crossDir;
                rp1 = intersection;
                rp1[crossDim] -= crossDir;
                rp2[crossDim] -= crossDir;
                isInitialSegment = true;
            }else{
                // no more intersections for this segment
                segmentCallback(simulationCell.reducedToAbsolute(rp1), simulationCell.reducedToAbsolute(rp2), isInitialSegment);
                isInitialSegment = false;
                break;
            }
        }while(true);
        rp1 = rp2;
    }
}

json DXAJsonExporter::exportDislocationsToJson(
    const DislocationNetwork* network, 
    const SimulationCell* simulationCell
){
    json dislocations;
    const auto& segments = network->segments();
    std::vector<const DislocationSegment*> validSegments;
    validSegments.reserve(segments.size());
    for(const auto* segment : segments){
        if(segment && !segment->isDegenerate()){
            validSegments.push_back(segment);
        }
    }

    json dataArray = json::array();
    double totalLength = 0.0;
    int totalPoints = 0;
    double maxLength = 0.0;
    double minLength = std::numeric_limits<double>::max();
    int globalChunkId = 0;

    auto saveChunk = [&](const std::vector<Point3>& chunk, const DislocationSegment* originalSegment){
        json segmentJson;
        json points = json::array();

        segmentJson["segment_id"] = globalChunkId++;
        double chunkLength = 0.0;
        for(size_t pointIdx = 0; pointIdx < chunk.size(); ++pointIdx){
            points.push_back({ chunk[pointIdx].x(), chunk[pointIdx].y(), chunk[pointIdx].z() });
            if(pointIdx > 0){
			chunkLength += (chunk[pointIdx] - chunk[pointIdx - 1]).length();
            }
        }

        segmentJson["points"] = points;
        segmentJson["length"] = chunkLength;
        segmentJson["num_points"] = chunk.size();

        Vector3 burgers = originalSegment->burgersVector.localVec();
        segmentJson["burgers_vector"] = { burgers.x(), burgers.y(), burgers.z() };
        segmentJson["magnitude"] = burgers.length();
        segmentJson["fractional"] = getBurgersVectorString(burgers);

        dataArray.push_back(segmentJson);

        totalLength += chunkLength;
        totalPoints += chunk.size();
        maxLength = std::max(maxLength, chunkLength);
        minLength = std::min(minLength, chunkLength);
    };

    for(size_t segmentId = 0; segmentId < validSegments.size(); ++segmentId){
        const auto* segment = validSegments[segmentId];
        if(simulationCell){
            std::vector<Point3> currentChunk;
            clipDislocationLine(segment->line, *simulationCell, 
                [&](const Point3& p1, const Point3& p2, bool isInitialSegment){
                    if(isInitialSegment && !currentChunk.empty()){
                        saveChunk(currentChunk, segment);
                        currentChunk.clear();
                    }

                    if(currentChunk.empty()){
                        currentChunk.push_back(p1);
                    }

                    currentChunk.push_back(p2);
            });

            if(!currentChunk.empty()){
                saveChunk(currentChunk, segment);
            }
        }else{
            std::vector<Point3> rawChunk(segment->line.begin(), segment->line.end());
            if(!rawChunk.empty()){
                saveChunk(rawChunk, segment);
            }
        }
    }

    if(dataArray.empty()){
        minLength = 0.0;
    }

    dislocations["main_listing"] = {
        { "dislocations", static_cast<int>(dataArray.size()) },
        { "total_points", totalPoints },
        { "average_segment_length", dataArray.empty() ? 0.0 : totalLength / dataArray.size() },
        { "max_segment_length", maxLength },
        { "min_segment_length", minLength },
        { "total_length", totalLength }
    };
    dislocations["sub_listings"] = { { "dislocation_segments", dataArray } };
    dislocations["export"]["DislocationExporter"]["segments"] = dataArray;

    dislocations["sub_listings"]["junction_information"] = getJunctionInformation(network);
    dislocations["sub_listings"]["circuit_information"] = getCircuitInformation(network);
    if(simulationCell){
        dislocations["sub_listings"]["network_statistics"] = getNetworkStatistics(network, simulationCell->volume3D());
    }

    return dislocations;
}


void DXAJsonExporter::writeDislocationsMsgpackToFile(
    const DislocationNetwork* network,
    const SimulationCell& simulationCell,
    const std::string& filePath
){
    if(network == nullptr) return;
    const json data = exportDislocationsToJson(network, &simulationCell);
    JsonUtils::writeJsonMsgpackToFile(data, filePath, false);
}


template <typename MeshType>
json DXAJsonExporter::getMeshData(
    const MeshType& mesh,
    const StructureAnalysis& structureAnalysis,
    bool includeTopologyInfo,
    const InterfaceMesh* interfaceMeshForTopology
){
    json meshData;
    const auto& originalVertices = mesh.vertices(); 
    const auto& originalFaces = mesh.faces();
    const auto& cell = structureAnalysis.context().simCell;

    std::vector<Point3> exportPoints;
    exportPoints.reserve(originalVertices.size());
    std::vector<int> originalToExportVertexMap(originalVertices.size());
    for(size_t i = 0; i < originalVertices.size(); ++i){
        exportPoints.push_back(originalVertices[i]->pos());
        originalToExportVertexMap[i] = i;
    }

    std::vector<std::vector<int>> exportFaces;
    exportFaces.reserve(originalFaces.size());
    for(const auto* face : originalFaces){
        if (!face || !face->edges()) continue;
        std::vector<int> faceVertexIndices;
        std::vector<Point3> faceVertexPositions;
        auto* startEdge = face->edges();
        auto* currentEdge = startEdge;
        do{
            faceVertexIndices.push_back(currentEdge->vertex1()->index());
            faceVertexPositions.push_back(currentEdge->vertex1()->pos());
            currentEdge = currentEdge->nextFaceEdge();
        }while(currentEdge != startEdge);

        cell.unwrapPositions(faceVertexPositions.data(), faceVertexPositions.size());

        std::vector<int> newFaceIndices;
        for(size_t i = 0; i < faceVertexIndices.size(); ++i){
            int originalIndex = faceVertexIndices[i];
            const Point3& originalPos = originalVertices[originalIndex]->pos();
            const Point3& unwrappedPos = faceVertexPositions[i];
			if(!originalPos.equals(unwrappedPos, 1e-6)){
                newFaceIndices.push_back(exportPoints.size());
                exportPoints.push_back(unwrappedPos);
            }else{
                newFaceIndices.push_back(originalToExportVertexMap[originalIndex]);
            }
        }
        exportFaces.push_back(newFaceIndices);
    }
    
    meshData["main_listing"] = {
        {"total_nodes", static_cast<int>(exportPoints.size())},
        {"total_facets", static_cast<int>(exportFaces.size())}
    };

    json points = json::array();
    for(size_t i = 0; i < exportPoints.size(); ++i){
        const auto& pos = exportPoints[i];
        points.push_back({
            {"index", static_cast<int>(i)},
            {"position", {pos.x(), pos.y(), pos.z()}}
        });
    }

    json facets = json::array();
    for(const auto& faceIndices : exportFaces){
        assert(faceIndices.size() == 3 && "The mesh does not contain any triangular faces.");
        facets.push_back({
            {"vertices", faceIndices}
        });
    }

    meshData["sub_listings"] = {
        {"points", points},
        {"facets", facets}
    };
    meshData["export"]["MeshExporter"]["vertices"] = points;
    meshData["export"]["MeshExporter"]["facets"] = facets;
    
    if(includeTopologyInfo && interfaceMeshForTopology != nullptr){
        std::unordered_set<uint64_t> originalEdgeSet;
        originalEdgeSet.reserve(originalFaces.size() * 3);
        for(const auto* face : originalFaces){
            if(!face || !face->edges()) continue;
            auto* edge = face->edges();
            do{
                int v1 = edge->vertex1()->index();
                int v2 = edge->vertex2()->index();
                if(v1 > v2) std::swap(v1, v2);
                uint64_t key = (static_cast<uint64_t>(static_cast<uint32_t>(v1)) << 32)
                    | static_cast<uint32_t>(v2);
                originalEdgeSet.insert(key);
                edge = edge->nextFaceEdge();
            }while(edge != face->edges());
        }

        meshData["topology"] = {
            {"euler_characteristic", static_cast<int>(originalVertices.size()) - static_cast<int>(originalEdgeSet.size()) + static_cast<int>(originalFaces.size())},
            {"is_completely_good", interfaceMeshForTopology->isCompletelyGood()},
            {"is_completely_bad", interfaceMeshForTopology->isCompletelyBad()}
        };
    }
    
    return meshData;
}

template <typename MeshType>
void DXAJsonExporter::writeMeshMsgpackToFile(
    const MeshType& mesh,
    const StructureAnalysis& structureAnalysis,
    bool includeTopologyInfo,
    const InterfaceMesh* interfaceMeshForTopology,
    const std::string& filePath
){
    const json data = getMeshData(mesh, structureAnalysis, includeTopologyInfo, interfaceMeshForTopology);
    JsonUtils::writeJsonMsgpackToFile(data, filePath, false);
}


void DXAJsonExporter::writeDefectMeshMsgpackToFile(
    const InterfaceMesh& interfaceMesh,
    const BurgersLoopBuilder& tracer,
    const StructureAnalysis& structureAnalysis,
    bool includeTopologyInfo,
    const std::string& filePath
){
    const json data = getMeshData(
        interfaceMesh, structureAnalysis, includeTopologyInfo, &interfaceMesh
    );
    JsonUtils::writeJsonMsgpackToFile(data, filePath, false);
}

void DXAJsonExporter::exportPTMData(
    const AnalysisContext& context,
    const std::vector<int>& ids,
    const std::string& outputFilename
){
    auto ptmProp = context.ptmOrientation;
    auto corrProp = context.correspondencesCode;
    if(!ptmProp || !corrProp) return;

    const bool includeStructureType = context.structureTypes && context.structureTypes->size() >= ids.size();
    json perAtom = json::array();
    for(size_t i = 0; i < ids.size(); ++i){
        json atom;
        atom["id"] = ids[i];
        atom["correspondences"] = static_cast<uint64_t>(corrProp->getInt64(i));
        if(includeStructureType){
            atom["structure_type"] = context.structureTypes->getInt(i);
        }
        json orient = json::array();
        for(int c = 0; c < 4; ++c) orient.push_back(ptmProp->getDoubleComponent(i, c));
        atom["orientation"] = orient;
        perAtom.push_back(atom);
    }

    json data;
    data["main_listing"] = {
        {"total_atoms", static_cast<int>(ids.size())},
        {"include_structure_type", includeStructureType}
    };
    data["per-atom-properties"] = perAtom;

    const std::string ptmPath = outputFilename + "_ptm_data.msgpack";
    if(JsonUtils::writeJsonMsgpackToFile(data, ptmPath, false)){
        spdlog::info("PTM data written to {}", ptmPath);
    }
}


void DXAJsonExporter::exportCoreAtoms(
    const LammpsParser::Frame& frame,
    const std::vector<uint8_t>& coreAtomFlags,
    const std::string& outputFilename
){
    if(coreAtomFlags.empty()) return;

    const size_t maxAtoms = std::min(static_cast<size_t>(frame.natoms), coreAtomFlags.size());

    json coreAtomsArray = json::array();
    for(size_t atomIdx = 0; atomIdx < maxAtoms; ++atomIdx){
        if(coreAtomFlags[atomIdx] == 0) continue;
        json atom;
        atom["id"] = frame.ids[atomIdx];
        if(atomIdx < frame.positions.size()){
            const auto& pos = frame.positions[atomIdx];
            atom["pos"] = { pos.x(), pos.y(), pos.z() };
        }
        coreAtomsArray.push_back(atom);
    }

    if(coreAtomsArray.empty()) return;

    json data;
    data["export"] = json::object();
    data["export"]["AtomisticExporter"] = json::object();
    data["export"]["AtomisticExporter"]["Core"] = coreAtomsArray;

    if(JsonUtils::writeJsonMsgpackToFile(data, outputFilename, false)){
        spdlog::info("Core atoms data written to {} ({} atoms)", outputFilename, coreAtomsArray.size());
    }
}

json DXAJsonExporter::vectorToJson(const Vector3& vector){
    json vectorJson;
    vectorJson["x"] = vector.x();
    vectorJson["y"] = vector.y();
    vectorJson["z"] = vector.z();
    return vectorJson;
}

json DXAJsonExporter::affineTransformationToJson(const AffineTransformation& transform){
    json transformJson = json::array();
    for(size_t i = 0; i < 3; ++i){
        json row = json::array();
        for(size_t j = 0; j < 3; ++j){
            row.push_back(transform(i, j));
        }
        transformJson.push_back(row);
    }
    return transformJson;
}

json DXAJsonExporter::simulationCellToJson(const SimulationCell& cell){
    json cellJson;
    
    cellJson["matrix"] = affineTransformationToJson(cell.matrix());
    cellJson["volume"] = cell.volume3D();
    cellJson["is_2d"] = cell.is2D();

    Vector3 a = cell.matrix().column(0);
    Vector3 b = cell.matrix().column(1);
    Vector3 c = cell.matrix().column(2);
    
    cellJson["lattice_vectors"] = {
        {"a", vectorToJson(a)},
        {"b", vectorToJson(b)},
        {"c", vectorToJson(c)}
    };
    
    cellJson["lattice_parameters"] = {
		{"a_length", a.length()},
		{"b_length", b.length()},
		{"c_length", c.length()}
    };
    
    return cellJson;
}

json DXAJsonExporter::getExtendedSimulationCellInfo(const SimulationCell& cell){
    json cellJson = simulationCellToJson(cell);
    
    const auto& pbcFlags = cell.pbcFlags();
    cellJson["periodic_boundary_conditions"] = {
        {"x", pbcFlags[0]},
        {"y", pbcFlags[1]}, 
        {"z", pbcFlags[2]}
    };
    
    Vector3 a = cell.matrix().column(0);
    Vector3 b = cell.matrix().column(1);
    Vector3 c = cell.matrix().column(2);
    
    cellJson["angles"] = {
        {"alpha", calculateAngle(b, c)},
        {"beta", calculateAngle(a, c)},
        {"gamma", calculateAngle(a, b)}
    };
    
    cellJson["reciprocal_lattice"] = {
        {"matrix", affineTransformationToJson(cell.inverseMatrix())},
        {"volume", 1.0 / cell.volume3D()}
    };
    
    cellJson["dimensionality"] = {
        {"is_2d", cell.is2D()},
        {"effective_dimensions", cell.is2D() ? 2 : 3}
    };
    
    return cellJson;
}

std::string DXAJsonExporter::getBurgersVectorString(const Vector3& burgers){
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3);
    oss << "[" << burgers.x() << " " << burgers.y() << " " << burgers.z() << "]";
    return oss.str();
}

void DXAJsonExporter::exportForStructureIdentification(
    const LammpsParser::Frame &frame,
    const StructureAnalysis& structureAnalysis,
    const std::string& outputFilename
){
    const size_t N = frame.natoms;

    std::vector<int> atomStructureTypes(N);
    for(size_t i = 0; i < N; ++i){
        atomStructureTypes[i] = static_cast<int>(
            structureAnalysis.context().structureTypes->getInt(static_cast<int>(i))
        );
    }

    json perAtomProperties = json::array();
    if(structureAnalysis.usingPTM()){
        const auto& ctx = structureAnalysis.context();
        const auto ptmOrientation = ctx.ptmOrientation;
        const auto correspondences = ctx.correspondencesCode;

        for(size_t i = 0; i < N; ++i){
            int atomId = i < frame.ids.size()
                ? frame.ids[i]
                : static_cast<int>(i);

            json atom;
            atom["id"] = atomId;
            atom["structure_type"] = atomStructureTypes[i];

            if(i < static_cast<size_t>(frame.positions.size())){
                const auto& pos = frame.positions[i];
                atom["pos"] = {pos.x(), pos.y(), pos.z()};
            }else{
                atom["pos"] = {0.0, 0.0, 0.0};
            }

            if(correspondences){
                atom["correspondence"] = static_cast<uint64_t>(correspondences->getInt64(i));
            }

            if(ptmOrientation){
                atom["orientation"] = {
                    ptmOrientation->getDoubleComponent(i, 0),
                    ptmOrientation->getDoubleComponent(i, 1),
                    ptmOrientation->getDoubleComponent(i, 2),
                    ptmOrientation->getDoubleComponent(i, 3)
                };
            }

            perAtomProperties.push_back(atom);
        }
    }else{
        for(size_t i = 0; i < N; ++i){
            int atomId = i < frame.ids.size()
                ? frame.ids[i]
                : static_cast<int>(i);

            json atom;
            atom["id"] = atomId;
            atom["structure_type"] = atomStructureTypes[i];
            atom["structure_name"] = structureAnalysis.getStructureTypeName(atomStructureTypes[i]);

            if(i < static_cast<size_t>(frame.positions.size())){
                const auto& pos = frame.positions[i];
                atom["pos"] = {pos.x(), pos.y(), pos.z()};
            }else{
                atom["pos"] = {0.0, 0.0, 0.0};
            }

            perAtomProperties.push_back(atom);
        }
    }

    json structureStats;
    structureStats["main_listing"] = structureAnalysis.buildMainListing();
    structureStats["per-atom-properties"] = perAtomProperties;
    JsonUtils::writeJsonMsgpackToFile(structureStats, outputFilename + "_structure_analysis_stats.msgpack", false);
}

json DXAJsonExporter::getNetworkStatistics(const DislocationNetwork* network, double cellVolume){
    json stats;
    const auto& segments = network->segments();
    
    double totalLength = 0.0;
    int validSegments = 0;
    
    struct SegmentStats {
        double length = 0.0;
        int count = 0;
    };
    
    auto stats_result = tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, segments.size()),
        SegmentStats{},
        [&segments](const tbb::blocked_range<size_t>& r, SegmentStats val) -> SegmentStats {
            for(size_t i = r.begin(); i < r.end(); ++i){
                const auto* segment = segments[i];
                if(segment && !segment->isDegenerate()){
                    val.length += segment->calculateLength();
                    val.count++;
                }
            }
            return val;
        },
        [](const SegmentStats& a, const SegmentStats& b) -> SegmentStats {
            return {a.length + b.length, a.count + b.count};
        }
    );
    totalLength = stats_result.length;
    validSegments = stats_result.count;
    
    stats = {
        {"total_network_length", totalLength},
        {"segment_count", validSegments},
        {"junction_count", countJunctions(network)},
        {"dangling_segments", countDanglingSegments(network)},
        {"average_segment_length", validSegments > 0 ? totalLength / validSegments : 0.0},
        {"density", cellVolume > 0 ? totalLength / cellVolume : 0.0},
        {"total_segments_including_degenerate", static_cast<int>(segments.size())}
    };
    
    return stats;
}

json DXAJsonExporter::getJunctionInformation(const DislocationNetwork* network){
    json junctionInfo;
    const auto& segments = network->segments();
    
    std::map<int, int> junctionArmDistribution;
    int totalJunctions = 0;
    
    for(const auto* segment : segments){
        if(segment){
            int forwardArms = segment->forwardNode().countJunctionArms();
            int backwardArms = segment->backwardNode().countJunctionArms();
            
            if(forwardArms > 1){
                junctionArmDistribution[forwardArms]++;
                totalJunctions++;
            }
            if(backwardArms > 1){
                junctionArmDistribution[backwardArms]++;
                totalJunctions++;
            }
        }
    }
    
    junctionInfo = {
        {"total_junctions", totalJunctions},
        {"junction_arm_distribution", junctionArmDistribution}
    };
    
    return junctionInfo;
}

json DXAJsonExporter::getCircuitInformation(const DislocationNetwork* network){
    json circuitInfo;
    const auto& segments = network->segments();
    
    std::vector<int> edgeCounts;
    int totalCircuits = 0;
    int danglingCircuits = 0;
    int blockedCircuits = 0;
    
    for(const auto* segment : segments){
        if(segment){
            if(segment->forwardNode().circuit){
                auto* circuit = segment->forwardNode().circuit;
                edgeCounts.push_back(circuit->edgeCount);
                totalCircuits++;
                if(circuit->isDangling) danglingCircuits++;
                if(circuit->isCompletelyBlocked) blockedCircuits++;
            }
            
            if(segment->backwardNode().circuit){
                auto* circuit = segment->backwardNode().circuit;
                edgeCounts.push_back(circuit->edgeCount);
                totalCircuits++;
                if(circuit->isDangling) danglingCircuits++;
                if(circuit->isCompletelyBlocked) blockedCircuits++;
            }
        }
    }
    
    double averageEdgeCount = 0.0;
    if(!edgeCounts.empty()){
        averageEdgeCount = std::accumulate(edgeCounts.begin(), edgeCounts.end(), 0.0) / edgeCounts.size();
    }
    
    circuitInfo = {
        {"total_circuits", totalCircuits},
        {"dangling_circuits", danglingCircuits},
        {"blocked_circuits", blockedCircuits},
        {"average_edge_count", averageEdgeCount},
        {"edge_count_range", {
            {"min", edgeCounts.empty() ? 0 : *std::min_element(edgeCounts.begin(), edgeCounts.end())},
            {"max", edgeCounts.empty() ? 0 : *std::max_element(edgeCounts.begin(), edgeCounts.end())}
        }}
    };
    
    return circuitInfo;
}



int DXAJsonExporter::countJunctions(const DislocationNetwork* network){
    int junctions = 0;
    const auto& segments = network->segments();
    
    for(const auto* segment : segments){
        if(segment){
            if(!segment->forwardNode().isDangling()) junctions++;
            if(!segment->backwardNode().isDangling()) junctions++;
        }
    }
    
    return junctions / 2;
}

int DXAJsonExporter::countDanglingSegments(const DislocationNetwork* network){
    int dangling = 0;
    const auto& segments = network->segments();
    
    for(const auto* segment : segments){
        if(segment && (segment->forwardNode().isDangling() || segment->backwardNode().isDangling())){
            dangling++;
        }
    }
    
    return dangling;
}

double DXAJsonExporter::calculateAngle(const Vector3& a, const Vector3& b){
    double dot = a.dot(b);
	double magnitudes = a.length() * b.length();
    
    if(magnitudes == 0.0) return 0.0;
    
    double cosAngle = dot / magnitudes;
    cosAngle = std::max(-1.0, std::min(1.0, cosAngle));
    
    return std::acos(cosAngle) * 180.0 / PI;
}

}

template void Volt::DXAJsonExporter::writeMeshMsgpackToFile<Volt::InterfaceMesh>(
    const Volt::InterfaceMesh& mesh,
    const Volt::StructureAnalysis& structureAnalysis,
    bool includeTopologyInfo,
    const Volt::InterfaceMesh* interfaceMeshForTopology,
    const std::string& filePath
);
