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
#include <volt/utilities/parquet_line_writer.h>
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

struct BurgersFamilyMatch{
    const char* name;
    const char* label;
};

constexpr BurgersFamilyMatch OTHER_BURGERS_FAMILY{"Other", "Other"};

// Burgers vectors arrive in the cluster lattice frame defined by the reference
// lattice YAMLs (lattices/*.yml): cubic lattices use the conventional cubic cell
// with unit lattice constant, so prototypes are exact lattice-fraction vectors.
// Cubic matching is permutation/sign invariant (sorted absolute components).
struct CubicBurgersFamily{
    std::array<double, 3> sortedAbsComponents;
    BurgersFamilyMatch family;
};

const CubicBurgersFamily FCC_BURGERS_FAMILIES[] = {
    {{0.5, 0.5, 0.0}, {"1/2<110>", "1/2<110> (Perfect)"}},
    {{1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0}, {"1/6<112>", "1/6<112> (Shockley)"}},
    {{1.0 / 6.0, 1.0 / 6.0, 0.0}, {"1/6<110>", "1/6<110> (Stair-rod)"}},
    {{1.0 / 3.0, 0.0, 0.0}, {"1/3<100>", "1/3<100> (Hirth)"}},
    {{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}, {"1/3<111>", "1/3<111> (Frank)"}}
};

const CubicBurgersFamily BCC_BURGERS_FAMILIES[] = {
    {{0.5, 0.5, 0.5}, {"1/2<111>", "1/2<111> (Perfect)"}},
    {{1.0, 0.0, 0.0}, {"<100>", "<100>"}},
    {{1.0, 1.0, 0.0}, {"<110>", "<110>"}}
};

const CubicBurgersFamily SC_BURGERS_FAMILIES[] = {
    {{1.0, 0.0, 0.0}, {"<100>", "<100>"}},
    {{1.0, 1.0, 0.0}, {"<110>", "<110>"}},
    {{1.0, 1.0, 1.0}, {"<111>", "<111>"}}
};

// The hexagonal reference frame (lattices/hcp.yml, hex_diamond.yml) keeps the
// c axis along z with the FCC-equivalent length scale: in-plane nearest
// neighbor distance a = 1/sqrt(2), ideal c = 2*sqrt(2/3)*a. Families are
// rotation invariant about c, so they match on (basal, axial) magnitudes.
struct HexagonalBurgersFamily{
    double basalLength;
    double axialLength;
    BurgersFamilyMatch family;
};

const HexagonalBurgersFamily HCP_BURGERS_FAMILIES[] = {
    {0.7071067811865476, 0.0, {"1/3<1-210>", "1/3<1-210> (Perfect basal)"}},
    {0.4082482904638630, 0.0, {"1/3<1-100>", "1/3<1-100> (Shockley)"}},
    {0.0, 1.1547005383792515, {"<0001>", "<0001> (Perfect c)"}},
    {0.0, 0.5773502691896257, {"1/2<0001>", "1/2<0001> (Partial c)"}},
    {0.7071067811865476, 1.1547005383792515, {"1/3<1-213>", "1/3<1-213> (Perfect c+a)"}}
};

template<std::size_t N>
BurgersFamilyMatch matchCubicBurgersFamily(
    const Vector3& localBurgers,
    const CubicBurgersFamily (&families)[N],
    double tolerance
){
    std::array<double, 3> components{
        std::abs(localBurgers.x()),
        std::abs(localBurgers.y()),
        std::abs(localBurgers.z())
    };
    std::sort(components.begin(), components.end(), std::greater<double>());
    for(const auto& candidate : families){
        if(std::abs(components[0] - candidate.sortedAbsComponents[0]) < tolerance
            && std::abs(components[1] - candidate.sortedAbsComponents[1]) < tolerance
            && std::abs(components[2] - candidate.sortedAbsComponents[2]) < tolerance){
            return candidate.family;
        }
    }
    return OTHER_BURGERS_FAMILY;
}

BurgersFamilyMatch matchHexagonalBurgersFamily(const Vector3& localBurgers, double tolerance){
    const double basal = std::hypot(localBurgers.x(), localBurgers.y());
    const double axial = std::abs(localBurgers.z());
    for(const auto& candidate : HCP_BURGERS_FAMILIES){
        if(std::abs(basal - candidate.basalLength) < tolerance
            && std::abs(axial - candidate.axialLength) < tolerance){
            return candidate.family;
        }
    }
    return OTHER_BURGERS_FAMILY;
}

BurgersFamilyMatch classifyBurgersFamily(const Vector3& localBurgers, const std::string& crystalStructure){
    LatticeStructureType structure = LATTICE_OTHER;
    if(crystalStructure.empty() || !parseLatticeStructureType(crystalStructure, structure)){
        return OTHER_BURGERS_FAMILY;
    }

    constexpr double tolerance = 0.01;
    switch(structure){
        case LATTICE_FCC:
        case LATTICE_CUBIC_DIAMOND:
            return matchCubicBurgersFamily(localBurgers, FCC_BURGERS_FAMILIES, tolerance);
        case LATTICE_BCC:
            return matchCubicBurgersFamily(localBurgers, BCC_BURGERS_FAMILIES, tolerance);
        case LATTICE_SC:
            return matchCubicBurgersFamily(localBurgers, SC_BURGERS_FAMILIES, tolerance);
        case LATTICE_HCP:
        case LATTICE_HEX_DIAMOND:
            return matchHexagonalBurgersFamily(localBurgers, tolerance);
        default:
            return OTHER_BURGERS_FAMILY;
    }
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

    // Pass 2: build JSON document and persist as Parquet payload.
    json points = json::array();
    for(size_t i = 0; i < exportPoints.size(); ++i){
        points.push_back({
            {"index", static_cast<int64_t>(i)},
            {"position", {exportPoints[i].x(), exportPoints[i].y(), exportPoints[i].z()}}
        });
    }
    json facets = json::array();
    for(const auto& f : exportFaces){
        facets.push_back({{"vertices", {f.verts[0], f.verts[1], f.verts[2]}}});
    }

    json doc;
    doc["main_listing"] = {
        {"total_primary_tetrahedra", primaryTetrahedra},
        {"total_nodes", static_cast<int64_t>(exportPoints.size())},
        {"total_facets", static_cast<int64_t>(exportFaces.size())},
        {"boundary_facets", boundaryFacets},
        {"internal_facets", internalFacets}
    };
    doc["sub_listings"] = {{"points", points}, {"facets", facets}};
    doc["export"]["MeshExporter"] = {{"vertices", points}, {"facets", facets}};

    JsonUtils::writeJsonToParquet(doc, filePath);
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
    const int baseAtomFields = 7; // id, pos, structure_id, structure_type, structure_name, cluster_id, topology_name (optional)

    // Pass 2: build JSON document and persist as Parquet payload.
    (void)baseAtomFields;
    json regions = json::array();
    json atomisticExporter = json::object();
    for(const auto& [clusterId, atomIndices] : atomIndicesByCluster){
        const std::size_t rep = atomIndices.front();
        const int stype = context.structureTypes ? context.structureTypes->getInt(rep) : static_cast<int>(StructureType::OTHER);
        const Cluster* cluster = structureAnalysis.clusterGraph().findCluster(clusterId);
        const std::string topo = cluster && !cluster->topologyName.empty()
            ? cluster->topologyName : topologyNameForAtomExport(structureAnalysis, rep, stype);
        const std::string clusterName = "Cluster " + std::to_string(clusterId);
        regions.push_back({
            {"cluster_id", clusterId},
            {"cluster_name", clusterName},
            {"atom_count", static_cast<int64_t>(atomIndices.size())},
            {"structure_id", stype},
            {"structure_name", structureTypeNameForExport(stype)},
            {"topology_name", topo}
        });

        json atoms = json::array();
        for(std::size_t atomIndex : atomIndices){
            const int astype = context.structureTypes ? context.structureTypes->getInt(atomIndex) : static_cast<int>(StructureType::OTHER);
            const int clustId = context.atomClusters ? context.atomClusters->getInt(atomIndex) : 0;
            const std::string atopo = topologyNameForAtomExport(structureAnalysis, atomIndex, astype);
            const auto& pos = atomIndex < frame.positions.size() ? frame.positions[atomIndex] : Point3::Origin();
            json atom = {
                {"id", atomIndex < frame.ids.size() ? frame.ids[atomIndex] : static_cast<int>(atomIndex)},
                {"pos", {pos.x(), pos.y(), pos.z()}},
                {"structure_id", astype},
                {"structure_type", astype},
                {"structure_name", structureTypeNameForExport(astype)},
                {"cluster_id", clustId}
            };
            if(!atopo.empty()) atom["topology_name"] = atopo;
            atoms.push_back(std::move(atom));
        }
        atomisticExporter[clusterName] = std::move(atoms);
    }

    json doc;
    doc["main_listing"] = {
        {"total_atoms", static_cast<int64_t>(context.atomCount())},
        {"coherent_region_count", static_cast<int64_t>(atomIndicesByCluster.size())},
        {"assigned_atoms", assignedAtoms},
        {"unassigned_atoms", unassignedAtoms},
        {"largest_region_size", largestRegionSize}
    };
    doc["sub_listings"] = {{"coherent_crystalline_regions", regions}};
    doc["export"]["AtomisticExporter"] = std::move(atomisticExporter);

    JsonUtils::writeJsonToParquet(doc, filePath);
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

BurgersFamily classifyBurgersFamily(const Vector3& localBurgers, const std::string& crystalStructure){
    const BurgersFamilyMatch match = Detail::classifyBurgersFamily(localBurgers, crystalStructure);
    return BurgersFamily{match.name, match.label};
}

void streamDislocationsToFile(
    const std::string& linesFilePath,
    const std::string& summaryFilePath,
    const DislocationNetwork* network,
    const SimulationCell* simulationCell,
    const DislocationsExportOptions& options
){
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
        int clusterId = 0;
        std::string crystalStructure;
        std::string burgersFamily;
        std::string burgersFamilyLabel;
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
            const Cluster* chunkCluster = segment->burgersVector.cluster();
            c.clusterId = chunkCluster ? chunkCluster->id : 0;
            c.crystalStructure = chunkCluster
                ? (!chunkCluster->topologyName.empty()
                    ? chunkCluster->topologyName
                    : structureTypeNameForExport(chunkCluster->structure))
                : std::string();
            const BurgersFamilyMatch family = Detail::classifyBurgersFamily(c.burgersLocal, c.crystalStructure);
            c.burgersFamily = family.name;
            c.burgersFamilyLabel = family.label;
            // Unclassified vectors keep the numeric label so distinct "Other"
            // vectors stay distinguishable in the charts.
            const std::string summaryKey = c.burgersFamily != "Other"
                ? c.burgersFamilyLabel
                : burgersVectorLabel(c.burgersLocal);
            chunks.push_back(std::move(c));

            totalLength += len;
            totalPoints += pts.size();
            maxLength = std::max(maxLength, len);
            minLength = std::min(minLength, len);
            auto& s = burgersSummary[summaryKey];
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

    // One row per segment in the standard line entity table (id, points,
    // per-segment property columns). VOLT discovers, queries and styles these
    // properties generically — no dislocation knowledge outside this plugin.
    auto vecAsList = [](const Vector3& v){ return std::vector<double>{v.x(), v.y(), v.z()}; };
    auto pointAsList = [](const Point3& p){ return std::vector<double>{p.x(), p.y(), p.z()}; };
    streamLinesToParquet(
        linesFilePath,
        chunks.size(),
        [&](std::size_t i, std::vector<Point3>& outPoints){ outPoints = chunks[i].points; },
        [&](ColumnarLineWriter& writer, std::size_t i){
            const auto& c = chunks[i];
            writer.field("length", c.length);
            writer.field("num_points", static_cast<std::int64_t>(c.points.size()));
            writer.field("magnitude", c.magnitude);
            writer.field("burgers_vector_local", vecAsList(c.burgersLocal));
            writer.field("burgers_vector_global", vecAsList(c.burgersGlobal));
            writer.field("crystal_structure", c.crystalStructure);
            writer.field("burgers_family", c.burgersFamily);
            writer.field("burgers_family_label", c.burgersFamilyLabel);
            writer.field("cluster_id", static_cast<std::int64_t>(c.clusterId));
            writer.field("head_vertex", pointAsList(c.points.empty() ? Point3::Origin() : c.points.front()));
            writer.field("tail_vertex", pointAsList(c.points.empty() ? Point3::Origin() : c.points.back()));
        }
    );

    // Network-level statistics and chart data ride a separate JSON-payload
    // summary file; they are plugin-specific aggregates, not line entities.
    json burgersLabels = json::array(), segmentCounts = json::array(), burgersLengths = json::array();
    for(const auto& [label, s] : burgersSummary){
        burgersLabels.push_back(label);
        segmentCounts.push_back(s.first);
        burgersLengths.push_back(s.second);
    }

    json doc;
    doc["export"]["ChartExporter"] = {
        {"burgers_counts", {{"burgers_vector", burgersLabels}, {"segment_count", segmentCounts}}},
        {"burgers_lengths", {{"burgers_vector", burgersLabels}, {"total_length", burgersLengths}}}
    };
    doc["main_listing"] = {
        {"dislocations", static_cast<int64_t>(chunks.size())},
        {"total_points", static_cast<int64_t>(totalPoints)},
        {"average_segment_length", chunks.empty() ? 0.0 : totalLength / chunks.size()},
        {"max_segment_length", maxLength},
        {"min_segment_length", minLength},
        {"total_length", totalLength}
    };

    json subListings;
    const double cellVolume = simulationCell ? simulationCell->volume3D() : 0.0;
    if(options.exportDislocationNetworkStats)
        subListings["network_statistics"] = getNetworkStatistics(network, cellVolume);
    if(options.exportCircuitInformation)
        subListings["circuit_information"] = getCircuitInformation(network);
    if(options.exportJunctions)
        subListings["junction_information"] = getJunctionInformation(network);
    doc["sub_listings"] = std::move(subListings);

    JsonUtils::writeJsonToParquet(doc, summaryFilePath);
}

void streamDefectMeshToFile(
    const std::string& filePath,
    const InterfaceMesh& interfaceMesh,
    const StructureAnalysis& structureAnalysis,
    bool includeTopologyInfo
){
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

    // Build JSON document and persist as Parquet payload.
    json points = json::array();
    for(size_t i = 0; i < exportPoints.size(); ++i){
        points.push_back({
            {"index", static_cast<int64_t>(i)},
            {"position", {exportPoints[i].x(), exportPoints[i].y(), exportPoints[i].z()}}
        });
    }
    json facets = json::array();
    for(const auto& f : exportFaces){
        facets.push_back({{"vertices", {f[0], f[1], f[2]}}});
    }

    json doc;
    doc["main_listing"] = {
        {"total_nodes", static_cast<int64_t>(exportPoints.size())},
        {"total_facets", static_cast<int64_t>(exportFaces.size())}
    };
    doc["sub_listings"] = {{"points", points}, {"facets", facets}};
    doc["export"]["MeshExporter"] = {{"vertices", points}, {"facets", facets}};

    if(includeTopologyInfo){
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
        doc["topology"] = {
            {"euler_characteristic", static_cast<int64_t>(originalVertices.size()) - static_cast<int64_t>(edgeSet.size()) + static_cast<int64_t>(originalFaces.size())},
            {"is_completely_good", interfaceMesh.isCompletelyGood()},
            {"is_completely_bad", interfaceMesh.isCompletelyBad()}
        };
    }

    JsonUtils::writeJsonToParquet(doc, filePath);
}

}
