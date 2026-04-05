#include <volt/helpers/dxa_serialization.h>
#include <volt/helpers/burgers_circuit.h>
#include <cassert>
#include <functional>
#include <iomanip>
#include <numeric>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <array>
#include <cstdint>
#include <map>
#include <volt/utilities/json_utils.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

namespace Volt::DxaSerialization {

namespace Detail {

Vector3 getGlobalBurgersVector(const ClusterVector& burgersVector){
    if(burgersVector.cluster() == nullptr){
        return burgersVector.localVec();
    }
    return burgersVector.toSpatialVector();
}

}  // namespace Detail

using namespace Detail;

std::string getBurgersVectorString(const Vector3& burgers);
json getNetworkStatistics(const DislocationNetwork* network, double cellVolume);
json getJunctionInformation(const DislocationNetwork* network);
json getCircuitInformation(const DislocationNetwork* network);
int countJunctions(const DislocationNetwork* network);
int countDanglingSegments(const DislocationNetwork* network);

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

json buildDislocationsJson(
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

        const Vector3 burgersLocal = originalSegment->burgersVector.localVec();
        const Vector3 burgersGlobal = getGlobalBurgersVector(originalSegment->burgersVector);
        segmentJson["burgers_vector"] = { burgersLocal.x(), burgersLocal.y(), burgersLocal.z() };
        segmentJson["burgers_vector_local"] = { burgersLocal.x(), burgersLocal.y(), burgersLocal.z() };
        segmentJson["burgers_vector_global"] = { burgersGlobal.x(), burgersGlobal.y(), burgersGlobal.z() };
        segmentJson["burgers"] = {
            {"vector", { burgersLocal.x(), burgersLocal.y(), burgersLocal.z() }},
            {"vector_local", { burgersLocal.x(), burgersLocal.y(), burgersLocal.z() }},
            {"vector_global", { burgersGlobal.x(), burgersGlobal.y(), burgersGlobal.z() }},
            {"magnitude", burgersLocal.length()}
        };
        segmentJson["magnitude"] = burgersLocal.length();
        segmentJson["fractional"] = getBurgersVectorString(burgersLocal);

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


json buildMeshJson(
    const InterfaceMesh& mesh,
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

json buildDefectMeshJson(
    const InterfaceMesh& interfaceMesh,
    const StructureAnalysis& structureAnalysis,
    bool includeTopologyInfo
){
    return buildMeshJson(interfaceMesh, structureAnalysis, includeTopologyInfo, &interfaceMesh);
}

std::string getBurgersVectorString(const Vector3& burgers){
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3);
    oss << "[" << burgers.x() << " " << burgers.y() << " " << burgers.z() << "]";
    return oss.str();
}

json getNetworkStatistics(const DislocationNetwork* network, double cellVolume){
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

json getJunctionInformation(const DislocationNetwork* network){
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

json getCircuitInformation(const DislocationNetwork* network){
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



int countJunctions(const DislocationNetwork* network){
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

int countDanglingSegments(const DislocationNetwork* network){
    int dangling = 0;
    const auto& segments = network->segments();
    
    for(const auto* segment : segments){
        if(segment && (segment->forwardNode().isDangling() || segment->backwardNode().isDangling())){
            dangling++;
        }
    }
    
    return dangling;
}

}
