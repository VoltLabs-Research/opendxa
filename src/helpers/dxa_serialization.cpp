#include <volt/helpers/dxa_serialization.h>
#include <volt/helpers/burgers_circuit.h>
#include <volt/structures/crystal_structure_types.h>
#include <cassert>
#include <functional>
#include <numeric>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <array>
#include <cstdint>
#include <map>
#include <set>
#include <string>
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

std::string burgersVectorLabel(const Vector3& vector){
    return "[" + std::to_string(vector.x()) + ", "
        + std::to_string(vector.y()) + ", "
        + std::to_string(vector.z()) + "]";
}

std::string structureTypeNameForExport(int structureType){
    return structureTypeName(structureType);
}

int atomIdForExport(const LammpsParser::Frame& frame, std::size_t atomIndex){
    return atomIndex < frame.ids.size()
        ? frame.ids[atomIndex]
        : static_cast<int>(atomIndex);
}

Point3 atomPositionForExport(const LammpsParser::Frame& frame, std::size_t atomIndex){
    if(atomIndex < frame.positions.size()){
        return frame.positions[atomIndex];
    }
    return Point3::Origin();
}

std::string topologyNameForAtomExport(const StructureAnalysis& structureAnalysis, std::size_t atomIndex, int structureType){
    if(const Cluster* cluster = structureAnalysis.atomCluster(static_cast<int>(atomIndex));
       cluster && !cluster->topologyName.empty()){
        return cluster->topologyName;
    }
    return {};
}

json buildAtomExportRecord(
    const LammpsParser::Frame& frame,
    const StructureAnalysis& structureAnalysis,
    std::size_t atomIndex
){
    const StructureContext& context = structureAnalysis.context();
    const int structureType = context.structureTypes
        ? context.structureTypes->getInt(atomIndex)
        : static_cast<int>(StructureType::OTHER);
    const int clusterId = context.atomClusters
        ? context.atomClusters->getInt(atomIndex)
        : 0;
    const Point3 position = atomPositionForExport(frame, atomIndex);
    json atom = {
        {"id", atomIdForExport(frame, atomIndex)},
        {"pos", {position.x(), position.y(), position.z()}},
        {"structure_id", structureType},
        {"structure_type", structureType},
        {"structure_name", structureTypeNameForExport(structureType)},
        {"cluster_id", clusterId}
    };

    const std::string topologyName = topologyNameForAtomExport(structureAnalysis, atomIndex, structureType);
    if(!topologyName.empty()){
        atom["topology_name"] = topologyName;
    }

    return atom;
}

std::array<std::uint64_t, 3> canonicalFaceKey(
    DelaunayTessellation::VertexHandle v0,
    DelaunayTessellation::VertexHandle v1,
    DelaunayTessellation::VertexHandle v2
){
    std::array<std::uint64_t, 3> key{
        static_cast<std::uint64_t>(v0),
        static_cast<std::uint64_t>(v1),
        static_cast<std::uint64_t>(v2)
    };
    std::sort(key.begin(), key.end());
    return key;
}

}  // namespace Detail

using namespace Detail;

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
    const SimulationCell* simulationCell,
    const DislocationsExportOptions& options
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
    std::map<std::string, std::pair<int, double>> burgersSummary;

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
        segmentJson["magnitude"] = burgersLocal.length();

        dataArray.push_back(segmentJson);

        totalLength += chunkLength;
        totalPoints += chunk.size();
        maxLength = std::max(maxLength, chunkLength);
        minLength = std::min(minLength, chunkLength);
        auto& summary = burgersSummary[burgersVectorLabel(burgersLocal)];
        summary.first += 1;
        summary.second += chunkLength;
    };

    for(size_t segmentId = 0; segmentId < validSegments.size(); ++segmentId){
        const auto* segment = validSegments[segmentId];
        if(simulationCell && options.clipPbcSegments){
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

    json burgersLabels = json::array();
    json burgersCounts = json::array();
    json burgersLengths = json::array();
    for(const auto& [label, summary] : burgersSummary){
        burgersLabels.push_back(label);
        burgersCounts.push_back(summary.first);
        burgersLengths.push_back(summary.second);
    }
    dislocations["export"]["ChartExporter"]["burgers_counts"] = {
        {"burgers_vector", burgersLabels},
        {"segment_count", burgersCounts}
    };
    dislocations["export"]["ChartExporter"]["burgers_lengths"] = {
        {"burgers_vector", burgersLabels},
        {"total_length", burgersLengths}
    };

    if(options.exportJunctions){
        dislocations["sub_listings"]["junction_information"] = getJunctionInformation(network);
    }
    if(options.exportCircuitInformation){
        dislocations["sub_listings"]["circuit_information"] = getCircuitInformation(network);
    }
    if(options.exportDislocationNetworkStats && simulationCell){
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

json buildDelaunayTessellationJson(const DelaunayTessellation& tessellation){
    std::unordered_map<std::uint64_t, int> exportVertexIndices;
    std::vector<Point3> exportPoints;
    std::set<std::array<std::uint64_t, 3>> exportedFaces;

    exportPoints.reserve(tessellation.numberOfPrimaryTetrahedra() * 2);

    json facets = json::array();
    int primaryTetrahedra = 0;
    int boundaryFacets = 0;
    int internalFacets = 0;

    for(DelaunayTessellation::CellHandle cell : tessellation.cells()){
        if(!tessellation.isValidCell(cell) || tessellation.isGhostCell(cell)){
            continue;
        }

        ++primaryTetrahedra;

        for(int face = 0; face < 4; ++face){
            const auto facet = tessellation.mirrorFacet(cell, face);
            const auto v0 = tessellation.cellVertex(cell, DelaunayTessellation::cellFacetVertexIndex(face, 0));
            const auto v1 = tessellation.cellVertex(cell, DelaunayTessellation::cellFacetVertexIndex(face, 1));
            const auto v2 = tessellation.cellVertex(cell, DelaunayTessellation::cellFacetVertexIndex(face, 2));

            if(!exportedFaces.insert(canonicalFaceKey(v0, v1, v2)).second){
                continue;
            }

            if(!tessellation.isValidCell(facet.first) || tessellation.isGhostCell(facet.first)){
                ++boundaryFacets;
            }else{
                ++internalFacets;
            }

            std::array<int, 3> faceIndices{};
            const std::array<DelaunayTessellation::VertexHandle, 3> handles{v0, v1, v2};
            for(std::size_t vertexIndex = 0; vertexIndex < handles.size(); ++vertexIndex){
                const std::uint64_t handleKey = static_cast<std::uint64_t>(handles[vertexIndex]);
                auto [entry, inserted] = exportVertexIndices.emplace(handleKey, static_cast<int>(exportPoints.size()));
                if(inserted){
                    exportPoints.push_back(tessellation.vertexPosition(handles[vertexIndex]));
                }
                faceIndices[vertexIndex] = entry->second;
            }

            facets.push_back({
                {"vertices", {faceIndices[0], faceIndices[1], faceIndices[2]}}
            });
        }
    }

    json points = json::array();
    for(std::size_t index = 0; index < exportPoints.size(); ++index){
        const Point3& point = exportPoints[index];
        points.push_back({
            {"index", static_cast<int>(index)},
            {"position", {point.x(), point.y(), point.z()}}
        });
    }

    json meshData;
    meshData["main_listing"] = {
        {"total_primary_tetrahedra", primaryTetrahedra},
        {"total_nodes", static_cast<int>(exportPoints.size())},
        {"total_facets", static_cast<int>(facets.size())},
        {"boundary_facets", boundaryFacets},
        {"internal_facets", internalFacets}
    };
    meshData["sub_listings"] = {
        {"points", points},
        {"facets", facets}
    };
    meshData["export"]["MeshExporter"]["vertices"] = points;
    meshData["export"]["MeshExporter"]["facets"] = facets;
    return meshData;
}

json buildStructureIdentificationJson(
    const LammpsParser::Frame& frame,
    const StructureAnalysis& structureAnalysis
){
    const StructureContext& context = structureAnalysis.context();

    std::map<int, std::vector<std::size_t>> atomIndicesByStructure;
    int clusteredAtoms = 0;
    for(std::size_t atomIndex = 0; atomIndex < context.atomCount(); ++atomIndex){
        const int structureType = context.structureTypes
            ? context.structureTypes->getInt(atomIndex)
            : static_cast<int>(StructureType::OTHER);
        atomIndicesByStructure[structureType].push_back(atomIndex);
        if(context.atomClusters && context.atomClusters->getInt(atomIndex) != 0){
            ++clusteredAtoms;
        }
    }

    json structures = json::array();
    json atomsByStructure = json::object();
    for(const auto& [structureType, atomIndices] : atomIndicesByStructure){
        if(atomIndices.empty()){
            continue;
        }

        const std::string structureName = structureTypeNameForExport(structureType);
        json atoms = json::array();
        for(std::size_t atomIndex : atomIndices){
            atoms.push_back(buildAtomExportRecord(frame, structureAnalysis, atomIndex));
        }

        atomsByStructure[structureName] = std::move(atoms);
        structures.push_back({
            {"structure_id", structureType},
            {"structure_name", structureName},
            {"atom_count", static_cast<int>(atomIndices.size())}
        });
    }

    json result;
    result["main_listing"] = {
        {"total_atoms", static_cast<int>(context.atomCount())},
        {"structure_count", static_cast<int>(structures.size())},
        {"clustered_atoms", clusteredAtoms},
        {"unclustered_atoms", static_cast<int>(context.atomCount()) - clusteredAtoms}
    };
    result["sub_listings"] = {
        {"structures", structures}
    };
    result["export"]["AtomisticExporter"] = std::move(atomsByStructure);
    return result;
}

json buildCoherentCrystallineRegionsJson(
    const LammpsParser::Frame& frame,
    const StructureAnalysis& structureAnalysis
){
    const StructureContext& context = structureAnalysis.context();

    std::map<int, std::vector<std::size_t>> atomIndicesByCluster;
    int unassignedAtoms = 0;
    for(std::size_t atomIndex = 0; atomIndex < context.atomCount(); ++atomIndex){
        const int clusterId = context.atomClusters
            ? context.atomClusters->getInt(atomIndex)
            : 0;
        if(clusterId == 0){
            ++unassignedAtoms;
            continue;
        }

        atomIndicesByCluster[clusterId].push_back(atomIndex);
    }

    json coherentRegions = json::array();
    json atomsByCluster = json::object();
    int largestRegionSize = 0;

    for(const auto& [clusterId, atomIndices] : atomIndicesByCluster){
        if(atomIndices.empty()){
            continue;
        }

        const std::size_t representativeAtomIndex = atomIndices.front();
        const int structureType = context.structureTypes
            ? context.structureTypes->getInt(representativeAtomIndex)
            : static_cast<int>(StructureType::OTHER);
        const std::string structureName = structureTypeNameForExport(structureType);
        const Cluster* cluster = structureAnalysis.clusterGraph().findCluster(clusterId);

        json atoms = json::array();
        for(std::size_t atomIndex : atomIndices){
            atoms.push_back(buildAtomExportRecord(frame, structureAnalysis, atomIndex));
        }

        const std::string clusterName = "Cluster " + std::to_string(clusterId);
        atomsByCluster[clusterName] = std::move(atoms);
        coherentRegions.push_back({
            {"cluster_id", clusterId},
            {"cluster_name", clusterName},
            {"atom_count", static_cast<int>(atomIndices.size())},
            {"structure_id", structureType},
            {"structure_name", structureName},
            {"topology_name", cluster && !cluster->topologyName.empty() ? cluster->topologyName : topologyNameForAtomExport(structureAnalysis, representativeAtomIndex, structureType)}
        });
        largestRegionSize = std::max(largestRegionSize, static_cast<int>(atomIndices.size()));
    }

    const int assignedAtoms = static_cast<int>(context.atomCount()) - unassignedAtoms;
    json result;
    result["main_listing"] = {
        {"total_atoms", static_cast<int>(context.atomCount())},
        {"coherent_region_count", static_cast<int>(coherentRegions.size())},
        {"assigned_atoms", assignedAtoms},
        {"unassigned_atoms", unassignedAtoms},
        {"largest_region_size", largestRegionSize}
    };
    result["sub_listings"] = {
        {"coherent_crystalline_regions", coherentRegions}
    };
    result["export"]["AtomisticExporter"] = std::move(atomsByCluster);
    return result;
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
    
    stats = json::array({
        json{
        {"total_network_length", totalLength},
        {"segment_count", validSegments},
        {"junction_count", countJunctions(network)},
        {"dangling_segments", countDanglingSegments(network)},
        {"average_segment_length", validSegments > 0 ? totalLength / validSegments : 0.0},
        {"density", cellVolume > 0 ? totalLength / cellVolume : 0.0},
        {"total_segments_including_degenerate", static_cast<int>(segments.size())}
        }
    });
    
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
    
    junctionInfo = json::array();
    if(junctionArmDistribution.empty()){
        junctionInfo.push_back({
            {"junction_arms", 0},
            {"junction_count", 0},
            {"total_junctions", totalJunctions}
        });
    }else{
        for(const auto& [junctionArms, junctionCount] : junctionArmDistribution){
            junctionInfo.push_back({
                {"junction_arms", junctionArms},
                {"junction_count", junctionCount},
                {"total_junctions", totalJunctions}
            });
        }
    }
    
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
    
    circuitInfo = json::array({
        json{
        {"total_circuits", totalCircuits},
        {"dangling_circuits", danglingCircuits},
        {"blocked_circuits", blockedCircuits},
        {"average_edge_count", averageEdgeCount},
        {"edge_count_min", edgeCounts.empty() ? 0 : *std::min_element(edgeCounts.begin(), edgeCounts.end())},
        {"edge_count_max", edgeCounts.empty() ? 0 : *std::max_element(edgeCounts.begin(), edgeCounts.end())}
        }
    });
    
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
