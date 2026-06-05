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
#include <fstream>

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

void streamDelaunayTessellationToFile(
    const std::string& filePath,
    const DelaunayTessellation& tessellation
){
    // Pass 1: collect unique vertices and faces
    std::unordered_map<std::uint64_t, int> vertexMap;
    std::vector<Point3> exportPoints;
    std::set<std::array<std::uint64_t, 3>> exportedFaces;
    struct FaceEntry{ std::array<int,3> verts; bool isBoundary; };
    std::vector<FaceEntry> exportFaces;
    exportPoints.reserve(tessellation.numberOfPrimaryTetrahedra() * 2);

    int primaryTetrahedra = 0, boundaryFacets = 0, internalFacets = 0;

    for(DelaunayTessellation::CellHandle cell : tessellation.cells()){
        if(!tessellation.isValidCell(cell) || tessellation.isGhostCell(cell)) continue;
        ++primaryTetrahedra;
        for(int face = 0; face < 4; ++face){
            const auto facet = tessellation.mirrorFacet(cell, face);
            const auto v0 = tessellation.cellVertex(cell, DelaunayTessellation::cellFacetVertexIndex(face, 0));
            const auto v1 = tessellation.cellVertex(cell, DelaunayTessellation::cellFacetVertexIndex(face, 1));
            const auto v2 = tessellation.cellVertex(cell, DelaunayTessellation::cellFacetVertexIndex(face, 2));
            if(!exportedFaces.insert(canonicalFaceKey(v0, v1, v2)).second) continue;

            const bool isBoundary = !tessellation.isValidCell(facet.first) || tessellation.isGhostCell(facet.first);
            isBoundary ? ++boundaryFacets : ++internalFacets;

            std::array<int,3> faceIdx{};
            for(int vi = 0; vi < 3; ++vi){
                const auto h = (vi == 0 ? v0 : vi == 1 ? v1 : v2);
                auto [it, ins] = vertexMap.emplace(static_cast<uint64_t>(h), static_cast<int>(exportPoints.size()));
                if(ins) exportPoints.push_back(tessellation.vertexPosition(h));
                faceIdx[vi] = it->second;
            }
            exportFaces.push_back({faceIdx, isBoundary});
        }
    }

    // Pass 2: stream to file
    std::ofstream of(filePath, std::ios::binary);
    MsgpackWriter w(of);

    w.write_map_header(3);

    w.write_key("main_listing");
    w.write_map_header(5);
    w.write_key("total_primary_tetrahedra"); w.write_int(primaryTetrahedra);
    w.write_key("total_nodes"); w.write_int(static_cast<int64_t>(exportPoints.size()));
    w.write_key("total_facets"); w.write_int(static_cast<int64_t>(exportFaces.size()));
    w.write_key("boundary_facets"); w.write_int(boundaryFacets);
    w.write_key("internal_facets"); w.write_int(internalFacets);

    w.write_key("sub_listings");
    w.write_map_header(2);
    w.write_key("points");
    w.write_array_header(static_cast<uint32_t>(exportPoints.size()));
    for(size_t i = 0; i < exportPoints.size(); ++i){
        w.write_map_header(2);
        w.write_key("index"); w.write_int(static_cast<int64_t>(i));
        w.write_key("position"); w.write_array_header(3);
        w.write_double(exportPoints[i].x()); w.write_double(exportPoints[i].y()); w.write_double(exportPoints[i].z());
    }
    w.write_key("facets");
    w.write_array_header(static_cast<uint32_t>(exportFaces.size()));
    for(const auto& f : exportFaces){
        w.write_map_header(1);
        w.write_key("vertices"); w.write_array_header(3);
        w.write_int(f.verts[0]); w.write_int(f.verts[1]); w.write_int(f.verts[2]);
    }

    w.write_key("export");
    w.write_map_header(1);
    w.write_key("MeshExporter");
    w.write_map_header(2);
    w.write_key("vertices");
    w.write_array_header(static_cast<uint32_t>(exportPoints.size()));
    for(size_t i = 0; i < exportPoints.size(); ++i){
        w.write_map_header(2);
        w.write_key("index"); w.write_int(static_cast<int64_t>(i));
        w.write_key("position"); w.write_array_header(3);
        w.write_double(exportPoints[i].x()); w.write_double(exportPoints[i].y()); w.write_double(exportPoints[i].z());
    }
    w.write_key("facets");
    w.write_array_header(static_cast<uint32_t>(exportFaces.size()));
    for(const auto& f : exportFaces){
        w.write_map_header(1);
        w.write_key("vertices"); w.write_array_header(3);
        w.write_int(f.verts[0]); w.write_int(f.verts[1]); w.write_int(f.verts[2]);
    }
    of.flush();
}

void streamCoherentCrystallineRegionsToFile(
    const std::string& filePath,
    const LammpsParser::Frame& frame,
    const StructureAnalysis& structureAnalysis
){
    const StructureContext& context = structureAnalysis.context();

    // Pass 1: group atoms by cluster
    std::map<int, std::vector<std::size_t>> atomIndicesByCluster;
    int unassignedAtoms = 0;
    for(std::size_t i = 0; i < context.atomCount(); ++i){
        const int clusterId = context.atomClusters ? context.atomClusters->getInt(i) : 0;
        if(clusterId == 0){ ++unassignedAtoms; continue; }
        atomIndicesByCluster[clusterId].push_back(i);
    }

    int largestRegionSize = 0;
    for(const auto& [_, indices] : atomIndicesByCluster)
        largestRegionSize = std::max(largestRegionSize, static_cast<int>(indices.size()));

    const int assignedAtoms = static_cast<int>(context.atomCount()) - unassignedAtoms;
    const int baseAtomFields = 6; // id, pos, structure_id, structure_name, cluster_id, topology_name (optional)

    // Pass 2: stream
    std::ofstream of(filePath, std::ios::binary);
    MsgpackWriter w(of);

    w.write_map_header(3);

    w.write_key("main_listing");
    w.write_map_header(5);
    w.write_key("total_atoms"); w.write_int(static_cast<int64_t>(context.atomCount()));
    w.write_key("coherent_region_count"); w.write_int(static_cast<int64_t>(atomIndicesByCluster.size()));
    w.write_key("assigned_atoms"); w.write_int(assignedAtoms);
    w.write_key("unassigned_atoms"); w.write_int(unassignedAtoms);
    w.write_key("largest_region_size"); w.write_int(largestRegionSize);

    w.write_key("sub_listings");
    w.write_map_header(1);
    w.write_key("coherent_crystalline_regions");
    w.write_array_header(static_cast<uint32_t>(atomIndicesByCluster.size()));
    for(const auto& [clusterId, atomIndices] : atomIndicesByCluster){
        const std::size_t rep = atomIndices.front();
        const int stype = context.structureTypes ? context.structureTypes->getInt(rep) : static_cast<int>(StructureType::OTHER);
        const Cluster* cluster = structureAnalysis.clusterGraph().findCluster(clusterId);
        const std::string topo = cluster && !cluster->topologyName.empty()
            ? cluster->topologyName : topologyNameForAtomExport(structureAnalysis, rep, stype);
        const std::string clusterName = "Cluster " + std::to_string(clusterId);
        w.write_map_header(6);
        w.write_key("cluster_id"); w.write_int(clusterId);
        w.write_key("cluster_name"); w.write_str(clusterName);
        w.write_key("atom_count"); w.write_int(static_cast<int64_t>(atomIndices.size()));
        w.write_key("structure_id"); w.write_int(stype);
        w.write_key("structure_name"); w.write_str(structureTypeNameForExport(stype));
        w.write_key("topology_name"); w.write_str(topo);
    }

    w.write_key("export");
    w.write_map_header(1);
    w.write_key("AtomisticExporter");
    w.write_map_header(static_cast<uint32_t>(atomIndicesByCluster.size()));
    for(const auto& [clusterId, atomIndices] : atomIndicesByCluster){
        const std::string clusterName = "Cluster " + std::to_string(clusterId);
        w.write_key(clusterName);
        w.write_array_header(static_cast<uint32_t>(atomIndices.size()));
        for(std::size_t atomIndex : atomIndices){
            const int stype = context.structureTypes ? context.structureTypes->getInt(atomIndex) : static_cast<int>(StructureType::OTHER);
            const int clustId = context.atomClusters ? context.atomClusters->getInt(atomIndex) : 0;
            const std::string topo = topologyNameForAtomExport(structureAnalysis, atomIndex, stype);
            const int nfields = topo.empty() ? (baseAtomFields - 1) : baseAtomFields;
            w.write_map_header(static_cast<uint32_t>(nfields));
            w.write_key("id"); w.write_int(atomIndex < frame.ids.size() ? frame.ids[atomIndex] : static_cast<int>(atomIndex));
            w.write_key("pos"); w.write_array_header(3);
            const auto& pos = atomIndex < frame.positions.size() ? frame.positions[atomIndex] : Point3::Origin();
            w.write_double(pos.x()); w.write_double(pos.y()); w.write_double(pos.z());
            w.write_key("structure_id"); w.write_int(stype);
            w.write_key("structure_type"); w.write_int(stype);
            w.write_key("structure_name"); w.write_str(structureTypeNameForExport(stype));
            w.write_key("cluster_id"); w.write_int(clustId);
            if(!topo.empty()){ w.write_key("topology_name"); w.write_str(topo); }
        }
    }
    of.flush();
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

// ============ STREAMING EXPORT ============

void streamDislocationsToFile(
    const std::string& filePath,
    const DislocationNetwork* network,
    const SimulationCell* simulationCell,
    const DislocationsExportOptions& options
){
    std::ofstream of(filePath, std::ios::binary);
    MsgpackWriter w(of);

    const auto& segments = network->segments();
    std::vector<const DislocationSegment*> validSegments;
    validSegments.reserve(segments.size());
    for(const auto* seg : segments){
        if(seg && !seg->isDegenerate()) validSegments.push_back(seg);
    }

    // Pre-compute all chunks to know counts for headers
    struct Chunk{
        std::vector<Point3> points;
        double length;
        Vector3 burgersLocal;
        Vector3 burgersGlobal;
        double magnitude;
    };
    std::vector<Chunk> chunks;
    chunks.reserve(validSegments.size() * 2);

    double totalLength = 0;
    int totalPoints = 0;
    double maxLength = 0, minLength = std::numeric_limits<double>::max();
    std::map<std::string, std::pair<int, double>> burgersSummary;

    for(const auto* segment : validSegments){
        auto emitChunk = [&](const std::vector<Point3>& pts){
            double len = 0;
            for(size_t i = 1; i < pts.size(); ++i)
                len += (pts[i] - pts[i-1]).length();

            Chunk c;
            c.points = pts;
            c.length = len;
            c.burgersLocal = segment->burgersVector.localVec();
            c.burgersGlobal = getGlobalBurgersVector(segment->burgersVector);
            c.magnitude = c.burgersLocal.length();
            chunks.push_back(std::move(c));

            totalLength += len;
            totalPoints += pts.size();
            maxLength = std::max(maxLength, len);
            minLength = std::min(minLength, len);
            auto& s = burgersSummary[burgersVectorLabel(c.burgersLocal)];
            s.first++; s.second += len;
        };

        if(simulationCell && options.clipPbcSegments){
            std::vector<Point3> currentChunk;
            clipDislocationLine(segment->line, *simulationCell,
                [&](const Point3& p1, const Point3& p2, bool isInitial){
                    if(isInitial && !currentChunk.empty()){
                        emitChunk(currentChunk);
                        currentChunk.clear();
                    }
                    if(currentChunk.empty()) currentChunk.push_back(p1);
                    currentChunk.push_back(p2);
                });
            if(!currentChunk.empty()) emitChunk(currentChunk);
        }else{
            if(!segment->line.empty())
                emitChunk(std::vector<Point3>(segment->line.begin(), segment->line.end()));
        }
    }
    if(chunks.empty()) minLength = 0;

    // Write msgpack structure directly
    // Root: map with 3 keys: "export", "main_listing", "sub_listings"
    w.write_map_header(3);

    // "export"
    w.write_key("export");
    w.write_map_header(2);
    {
        // "DislocationExporter"
        w.write_key("DislocationExporter");
        w.write_map_header(1);
        w.write_key("segments");
        w.write_array_header(static_cast<uint32_t>(chunks.size()));
        for(size_t i = 0; i < chunks.size(); ++i){
            const auto& c = chunks[i];
            w.write_map_header(8);
            w.write_key("segment_id"); w.write_int(static_cast<int64_t>(i));
            w.write_key("points");
            w.write_array_header(static_cast<uint32_t>(c.points.size()));
            for(const auto& p : c.points){
                w.write_array_header(3);
                w.write_double(p.x()); w.write_double(p.y()); w.write_double(p.z());
            }
            w.write_key("length"); w.write_double(c.length);
            w.write_key("num_points"); w.write_int(static_cast<int64_t>(c.points.size()));
            w.write_key("burgers_vector");
            w.write_array_header(3);
            w.write_double(c.burgersLocal.x()); w.write_double(c.burgersLocal.y()); w.write_double(c.burgersLocal.z());
            w.write_key("burgers_vector_local");
            w.write_array_header(3);
            w.write_double(c.burgersLocal.x()); w.write_double(c.burgersLocal.y()); w.write_double(c.burgersLocal.z());
            w.write_key("burgers_vector_global");
            w.write_array_header(3);
            w.write_double(c.burgersGlobal.x()); w.write_double(c.burgersGlobal.y()); w.write_double(c.burgersGlobal.z());
            w.write_key("magnitude"); w.write_double(c.magnitude);
        }

        // "ChartExporter"
        w.write_key("ChartExporter");
        w.write_map_header(2);
        w.write_key("burgers_counts");
        w.write_map_header(2);
        w.write_key("burgers_vector");
        w.write_array_header(static_cast<uint32_t>(burgersSummary.size()));
        for(const auto& [label, _] : burgersSummary) w.write_str(label);
        w.write_key("segment_count");
        w.write_array_header(static_cast<uint32_t>(burgersSummary.size()));
        for(const auto& [_, s] : burgersSummary) w.write_int(s.first);
        w.write_key("burgers_lengths");
        w.write_map_header(2);
        w.write_key("burgers_vector");
        w.write_array_header(static_cast<uint32_t>(burgersSummary.size()));
        for(const auto& [label, _] : burgersSummary) w.write_str(label);
        w.write_key("total_length");
        w.write_array_header(static_cast<uint32_t>(burgersSummary.size()));
        for(const auto& [_, s] : burgersSummary) w.write_double(s.second);
    }

    // "main_listing"
    w.write_key("main_listing");
    w.write_map_header(6);
    w.write_key("dislocations"); w.write_int(static_cast<int64_t>(chunks.size()));
    w.write_key("total_points"); w.write_int(static_cast<int64_t>(totalPoints));
    w.write_key("average_segment_length"); w.write_double(chunks.empty() ? 0.0 : totalLength / chunks.size());
    w.write_key("max_segment_length"); w.write_double(maxLength);
    w.write_key("min_segment_length"); w.write_double(minLength);
    w.write_key("total_length"); w.write_double(totalLength);

    // "sub_listings" - empty for streaming (data is in export)
    w.write_key("sub_listings");
    w.write_nil();

    of.flush();
}

void streamDefectMeshToFile(
    const std::string& filePath,
    const InterfaceMesh& interfaceMesh,
    const StructureAnalysis& structureAnalysis,
    bool includeTopologyInfo
){
    std::ofstream of(filePath, std::ios::binary);
    MsgpackWriter w(of);

    const auto& originalVertices = interfaceMesh.vertices();
    const auto& originalFaces = interfaceMesh.faces();
    const auto& cell = structureAnalysis.context().simCell;

    // Pre-compute export data (positions + face indices with PBC unwrapping)
    std::vector<Point3> exportPoints;
    exportPoints.reserve(originalVertices.size());
    std::vector<int> originalToExportMap(originalVertices.size());
    for(size_t i = 0; i < originalVertices.size(); ++i){
        exportPoints.push_back(originalVertices[i]->pos());
        originalToExportMap[i] = static_cast<int>(i);
    }

    std::vector<std::array<int, 3>> exportFaces;
    exportFaces.reserve(originalFaces.size());
    for(const auto* face : originalFaces){
        if(!face || !face->edges()) continue;
        std::array<int, 3> faceVerts{};
        std::array<Point3, 3> facePos{};
        auto* edge = face->edges();
        for(int i = 0; i < 3; ++i, edge = edge->nextFaceEdge()){
            faceVerts[i] = edge->vertex1()->index();
            facePos[i] = edge->vertex1()->pos();
        }

        cell.unwrapPositions(facePos.data(), 3);

        std::array<int, 3> newFaceVerts{};
        for(int i = 0; i < 3; ++i){
            const Point3& orig = originalVertices[faceVerts[i]]->pos();
            if(!orig.equals(facePos[i], 1e-6)){
                newFaceVerts[i] = static_cast<int>(exportPoints.size());
                exportPoints.push_back(facePos[i]);
            }else{
                newFaceVerts[i] = originalToExportMap[faceVerts[i]];
            }
        }
        exportFaces.push_back(newFaceVerts);
    }

    // Write msgpack
    int numKeys = includeTopologyInfo ? 4 : 3;
    w.write_map_header(numKeys);

    // "main_listing"
    w.write_key("main_listing");
    w.write_map_header(2);
    w.write_key("total_nodes"); w.write_int(static_cast<int64_t>(exportPoints.size()));
    w.write_key("total_facets"); w.write_int(static_cast<int64_t>(exportFaces.size()));

    // "sub_listings"
    w.write_key("sub_listings");
    w.write_map_header(2);
    w.write_key("points");
    w.write_array_header(static_cast<uint32_t>(exportPoints.size()));
    for(size_t i = 0; i < exportPoints.size(); ++i){
        w.write_map_header(2);
        w.write_key("index"); w.write_int(static_cast<int64_t>(i));
        w.write_key("position");
        w.write_array_header(3);
        w.write_double(exportPoints[i].x());
        w.write_double(exportPoints[i].y());
        w.write_double(exportPoints[i].z());
    }
    w.write_key("facets");
    w.write_array_header(static_cast<uint32_t>(exportFaces.size()));
    for(const auto& f : exportFaces){
        w.write_map_header(1);
        w.write_key("vertices");
        w.write_array_header(3);
        w.write_int(f[0]); w.write_int(f[1]); w.write_int(f[2]);
    }

    // "export"
    w.write_key("export");
    w.write_map_header(1);
    w.write_key("MeshExporter");
    w.write_map_header(2);
    w.write_key("vertices");
    w.write_array_header(static_cast<uint32_t>(exportPoints.size()));
    for(size_t i = 0; i < exportPoints.size(); ++i){
        w.write_map_header(2);
        w.write_key("index"); w.write_int(static_cast<int64_t>(i));
        w.write_key("position");
        w.write_array_header(3);
        w.write_double(exportPoints[i].x());
        w.write_double(exportPoints[i].y());
        w.write_double(exportPoints[i].z());
    }
    w.write_key("facets");
    w.write_array_header(static_cast<uint32_t>(exportFaces.size()));
    for(const auto& f : exportFaces){
        w.write_map_header(1);
        w.write_key("vertices");
        w.write_array_header(3);
        w.write_int(f[0]); w.write_int(f[1]); w.write_int(f[2]);
    }

    if(includeTopologyInfo){
        w.write_key("topology");
        w.write_map_header(3);
        std::unordered_set<uint64_t> edgeSet;
        edgeSet.reserve(originalFaces.size() * 3);
        for(const auto* face : originalFaces){
            if(!face || !face->edges()) continue;
            auto* edge = face->edges();
            do{
                int v1 = edge->vertex1()->index(), v2 = edge->vertex2()->index();
                if(v1 > v2) std::swap(v1, v2);
                edgeSet.insert((static_cast<uint64_t>(static_cast<uint32_t>(v1)) << 32) | static_cast<uint32_t>(v2));
                edge = edge->nextFaceEdge();
            }while(edge != face->edges());
        }
        w.write_key("euler_characteristic");
        w.write_int(static_cast<int64_t>(originalVertices.size()) - static_cast<int64_t>(edgeSet.size()) + static_cast<int64_t>(originalFaces.size()));
        w.write_key("is_completely_good"); w.write_bool(interfaceMesh.isCompletelyGood());
        w.write_key("is_completely_bad"); w.write_bool(interfaceMesh.isCompletelyBad());
    }

    of.flush();
}

}
