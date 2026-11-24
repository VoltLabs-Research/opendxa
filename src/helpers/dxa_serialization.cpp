#include <utility>
#include <volt/helpers/dxa_serialization.h>
#include <volt/helpers/burgers_circuit.h>
#include <volt/helpers/edge_key.h>
#include <volt/structures/crystal_structure_types.h>
#include <cassert>
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
#include <volt/utilities/parquet_mesh_writer.h>
#include <spdlog/spdlog.h>
#include <fmt/format.h>
#include <filesystem>
#include <fstream>

namespace Volt::DxaSerialization {

namespace {

constexpr const char* kParquetSourceKey = "__parquet_source__";

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

const BurgersFamily OTHER_BURGERS_FAMILY{"Other", "Other"};

struct CubicBurgersFamily{
    std::array<double, 3> sortedAbsComponents;
    BurgersFamily family;
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

struct HexagonalBurgersFamily{
    double basalLength;
    double axialLength;
    BurgersFamily family;
};

// HCP lengths in the units of the local Burgers vectors, where the nearest-neighbour distance is a = 1/√2 and
// c/a takes its ideal value √(8/3).
constexpr double kHcpA = 0.7071067811865476;         // 1/√2
constexpr double kHcpShockley = 0.4082482904638630;  // a/√3
constexpr double kHcpC = 1.1547005383792515;         // a·√(8/3)
constexpr double kHcpHalfC = 0.5773502691896257;     // c/2

const HexagonalBurgersFamily HCP_BURGERS_FAMILIES[] = {
    {kHcpA, 0.0, {"1/3<1-210>", "1/3<1-210> (Perfect basal)"}},
    {kHcpShockley, 0.0, {"1/3<1-100>", "1/3<1-100> (Shockley)"}},
    {0.0, kHcpC, {"<0001>", "<0001> (Perfect c)"}},
    {0.0, kHcpHalfC, {"1/2<0001>", "1/2<0001> (Partial c)"}},
    {kHcpA, kHcpC, {"1/3<1-213>", "1/3<1-213> (Perfect c+a)"}}
};

template<std::size_t N>
const BurgersFamily& matchCubicBurgersFamily(
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

const BurgersFamily& matchHexagonalBurgersFamily(const Vector3& localBurgers, double tolerance){
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

int structureTypeOf(const StructureContext& context, std::size_t atomIndex){
    return context.structureTypes
        ? context.structureTypes->getInt(atomIndex)
        : static_cast<int>(StructureType::OTHER);
}

std::array<std::uint64_t, 3> canonicalFaceKey(const std::array<DelaunayTessellation::VertexHandle, 3>& corners){
    std::array<std::uint64_t, 3> key{
        static_cast<std::uint64_t>(corners[0]),
        static_cast<std::uint64_t>(corners[1]),
        static_cast<std::uint64_t>(corners[2])
    };
    std::sort(key.begin(), key.end());
    return key;
}

// Past this many periodic crossings within one segment, the rest of it is kept as it is.
constexpr int kMaxCrossingsPerSegment = 10;
constexpr double kCrossingEpsilon = 1e-9;

}

std::vector<std::vector<Point3>> splitAtPeriodicBoundaries(const std::vector<Point3>& line, const SimulationCell& cell){
    std::vector<std::vector<Point3>> pieces;
    if(line.size() < 2) return pieces;

    // Works in reduced coordinates, with the line shifted so that it starts inside the cell. Each crossing of a
    // periodic face ends the current piece and wraps the rest of the line back by one cell.
    bool startsPiece = true;
    auto addSegment = [&](const Point3& reducedFrom, const Point3& reducedTo){
        if(startsPiece) pieces.push_back({cell.reducedToAbsolute(reducedFrom)});
        pieces.back().push_back(cell.reducedToAbsolute(reducedTo));
    };

    Point3 rp1 = cell.absoluteToReduced(line.front());
    Vector3 shiftVector = Vector3::Zero();
    for(size_t dimension = 0; dimension < 3; dimension++){
        if(cell.pbcFlags()[dimension]){
            double shift = -std::floor(rp1[dimension]);
            rp1[dimension] += shift;
            shiftVector[dimension] += shift;
        }
    }

    for(size_t i = 1; i < line.size(); ++i){
        Point3 rp2 = cell.absoluteToReduced(line[i]) + shiftVector;
        for(int crossings = 0; ; ++crossings){
            if(crossings == kMaxCrossingsPerSegment){
                addSegment(rp1, rp2);
                break;
            }

            size_t crossDim = -1;
            double crossDir = 0;
            double smallestT = std::numeric_limits<double>::max();
            for(size_t dimension = 0; dimension < 3; dimension++){
                if(!cell.pbcFlags()[dimension]) continue;
                int d = (int) std::floor(rp2[dimension]) - (int) std::floor(rp1[dimension]);
                if(d == 0) continue;

                double dr = rp2[dimension] - rp1[dimension];
                if(std::abs(dr) < kCrossingEpsilon) continue;

                double t = (d > 0) ? (std::ceil(rp1[dimension]) - rp1[dimension]) / dr
                                   : (std::floor(rp1[dimension]) - rp1[dimension]) / dr;
                if(t > kCrossingEpsilon && t < smallestT){
                    smallestT = t;
                    crossDim = dimension;
                    crossDir = (d > 0) ? 1.0 : -1.0;
                }
            }

            if(smallestT < (1.0 - kCrossingEpsilon)){
                Point3 intersection = rp1 + smallestT * (rp2 - rp1);
                intersection[crossDim] = std::round(intersection[crossDim]);
                addSegment(rp1, intersection);
                startsPiece = true;
                shiftVector[crossDim] -= crossDir;
                rp1 = intersection;
                rp1[crossDim] -= crossDir;
                rp2[crossDim] -= crossDir;
            }else{
                addSegment(rp1, rp2);
                startsPiece = false;
                break;
            }
        }
        rp1 = rp2;
    }
    return pieces;
}

void streamDelaunayTessellationToFile(
    const std::string& filePath,
    const DelaunayTessellation& tessellation
){
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
            std::array<DelaunayTessellation::VertexHandle, 3> corners;
            for(int corner = 0; corner < 3; ++corner){
                corners[corner] = tessellation.cellVertex(cell, DelaunayTessellation::cellFacetVertexIndex(face, corner));
            }
            if(!exportedFaces.insert(canonicalFaceKey(corners)).second) continue;

            const DelaunayTessellation::CellHandle neighbor = tessellation.mirrorFacet(cell, face).first;
            const bool isBoundary = !tessellation.isValidCell(neighbor) || tessellation.isGhostCell(neighbor);
            isBoundary ? ++boundaryFacets : ++internalFacets;

            std::array<int, 3> faceIdx{};
            for(int corner = 0; corner < 3; ++corner){
                const auto handle = corners[corner];
                auto [it, inserted] = vertexMap.emplace(static_cast<uint64_t>(handle), static_cast<int>(exportPoints.size()));
                if(inserted) exportPoints.push_back(tessellation.vertexPosition(handle));
                faceIdx[corner] = it->second;
            }
            exportFaces.push_back({faceIdx, isBoundary});
        }
    }

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
    json regions = json::array();
    json atomisticExporter = json::object();
    for(const auto& [clusterId, atomIndices] : atomIndicesByCluster){
        const int stype = structureTypeOf(context, atomIndices.front());
        const Cluster* cluster = structureAnalysis.clusterGraph().findCluster(clusterId);
        const std::string clusterName = "Cluster " + std::to_string(clusterId);
        regions.push_back({
            {"cluster_id", clusterId},
            {"cluster_name", clusterName},
            {"atom_count", static_cast<int64_t>(atomIndices.size())},
            {"structure_id", stype},
            {"structure_name", structureTypeName(stype)},
            {"topology_name", cluster ? cluster->topologyName : std::string()}
        });

        json atoms = json::array();
        for(std::size_t atomIndex : atomIndices){
            const int astype = structureTypeOf(context, atomIndex);
            const auto& pos = atomIndex < frame.positions.size() ? frame.positions[atomIndex] : Point3::Origin();
            atoms.push_back({
                {"id", atomIndex < frame.ids.size() ? frame.ids[atomIndex] : static_cast<int>(atomIndex)},
                {"pos", {pos.x(), pos.y(), pos.z()}},
                {"structure_id", astype},
                {"structure_name", structureTypeName(astype)},
                {"cluster_id", clusterId}
            });
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

namespace {

std::array<const DislocationNode*, 2> endNodes(const DislocationSegment& segment){
    return {&segment.forwardNode(), &segment.backwardNode()};
}

int countDanglingSegments(const DislocationNetwork& network){
    return static_cast<int>(std::ranges::count_if(network.segments(), [](const DislocationSegment* segment){
        return segment && (segment->forwardNode().isDangling() || segment->backwardNode().isDangling());
    }));
}

json getNetworkStatistics(const DislocationNetwork& network, double cellVolume){
    const auto& segments = network.segments();

    // Summed in segment order, so the totals are the same with any thread count.
    double totalLength = 0.0;
    int segmentCount = 0;
    for(const DislocationSegment* segment : segments){
        if(!segment || segment->isDegenerate()) continue;
        totalLength += segment->calculateLength();
        ++segmentCount;
    }

    return json::array({json{
        {"total_network_length", totalLength},
        {"segment_count", segmentCount},
        {"junction_count", static_cast<int>(network.junctionArmCounts().size())},
        {"dangling_segments", countDanglingSegments(network)},
        {"average_segment_length", segmentCount > 0 ? totalLength / segmentCount : 0.0},
        {"density", cellVolume > 0 ? totalLength / cellVolume : 0.0},
        {"total_segments_including_degenerate", static_cast<int>(segments.size())}
    }});
}

json getJunctionInformation(const DislocationNetwork& network){
    const std::vector<int> junctionArms = network.junctionArmCounts();
    std::map<int, int> junctionArmDistribution;
    for(const int arms : junctionArms){
        junctionArmDistribution[arms]++;
    }
    const int totalJunctions = static_cast<int>(junctionArms.size());

    json junctionInfo = json::array();
    if(junctionArmDistribution.empty()){
        junctionInfo.push_back({
            {"junction_arms", 0},
            {"junction_count", 0},
            {"total_junctions", totalJunctions}
        });
        return junctionInfo;
    }
    for(const auto& [junctionArms, junctionCount] : junctionArmDistribution){
        junctionInfo.push_back({
            {"junction_arms", junctionArms},
            {"junction_count", junctionCount},
            {"total_junctions", totalJunctions}
        });
    }
    return junctionInfo;
}

json getCircuitInformation(const DislocationNetwork& network){
    std::vector<int> edgeCounts;
    int danglingCircuits = 0;
    int blockedCircuits = 0;
    for(const DislocationSegment* segment : network.segments()){
        if(!segment) continue;
        for(const DislocationNode* node : endNodes(*segment)){
            const BurgersCircuit* circuit = node->circuit;
            if(!circuit) continue;
            edgeCounts.push_back(circuit->edgeCount);
            if(circuit->isDangling) danglingCircuits++;
            if(circuit->isCompletelyBlocked) blockedCircuits++;
        }
    }

    const double averageEdgeCount = edgeCounts.empty()
        ? 0.0
        : std::accumulate(edgeCounts.begin(), edgeCounts.end(), 0.0) / edgeCounts.size();

    return json::array({json{
        {"total_circuits", static_cast<int>(edgeCounts.size())},
        {"dangling_circuits", danglingCircuits},
        {"blocked_circuits", blockedCircuits},
        {"average_edge_count", averageEdgeCount},
        {"edge_count_min", edgeCounts.empty() ? 0 : *std::ranges::min_element(edgeCounts)},
        {"edge_count_max", edgeCounts.empty() ? 0 : *std::ranges::max_element(edgeCounts)}
    }});
}

}

BurgersFamily classifyBurgersFamily(const Vector3& localBurgers, const std::string& crystalStructure){
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

namespace {

// What every piece of a dislocation segment shares, computed once per segment.
struct SegmentProperties{
    Vector3 burgersLocal;
    Vector3 burgersGlobal;
    double magnitude;
    int clusterId;
    int segmentId;
    std::string crystalStructure;
    BurgersFamily family;
};

struct DislocationLine{
    std::vector<Point3> points;
    double length;
    SegmentProperties segment;
};

struct FamilyTally{
    std::string label;
    int segmentCount = 0;
    double totalLength = 0.0;
};

struct LineSummary{
    double totalLength = 0;
    int totalPoints = 0;
    double minLength = 0;
    double meanLength = 0;
    double maxLength = 0;
    std::map<std::string, FamilyTally> byBurgersVector;
    std::map<std::string, FamilyTally> byFamily;
};

std::vector<const DislocationSegment*> nonDegenerateSegments(const DislocationNetwork& network){
    std::vector<const DislocationSegment*> segments;
    segments.reserve(network.segments().size());
    for(const auto* segment : network.segments()){
        if(segment && !segment->isDegenerate()) segments.push_back(segment);
    }
    return segments;
}

// The name the Burgers family is looked up by: the cluster's topology, or its structure type when the topology
// has no name.
std::string crystalStructureOf(const Cluster* cluster){
    if(!cluster) return {};
    return !cluster->topologyName.empty() ? cluster->topologyName : std::string(structureTypeName(cluster->structure));
}

SegmentProperties segmentProperties(const DislocationSegment& segment){
    const Cluster* cluster = segment.burgersVector.cluster();
    SegmentProperties properties;
    properties.burgersLocal = segment.burgersVector.localVec();
    properties.burgersGlobal = getGlobalBurgersVector(segment.burgersVector);
    properties.magnitude = properties.burgersLocal.length();
    properties.clusterId = cluster ? cluster->id : 0;
    properties.segmentId = segment.id;
    properties.crystalStructure = crystalStructureOf(cluster);
    properties.family = classifyBurgersFamily(properties.burgersLocal, properties.crystalStructure);
    return properties;
}

double polylineLength(const std::vector<Point3>& points){
    double length = 0;
    for(size_t i = 1; i < points.size(); ++i){
        length += (points[i] - points[i - 1]).length();
    }
    return length;
}

std::vector<DislocationLine> collectLines(
    const std::vector<const DislocationSegment*>& segments,
    const SimulationCell& cell,
    bool clipPbcSegments
){
    std::vector<DislocationLine> lines;
    lines.reserve(segments.size() * 2);
    for(const DislocationSegment* segment : segments){
        const SegmentProperties properties = segmentProperties(*segment);
        auto addLine = [&](std::vector<Point3> points){
            const double length = polylineLength(points);
            lines.push_back({std::move(points), length, properties});
        };
        if(clipPbcSegments){
            for(auto& piece : splitAtPeriodicBoundaries(segment->line, cell)){
                addLine(std::move(piece));
            }
        }else if(!segment->line.empty()){
            addLine(segment->line);
        }
    }
    return lines;
}

// Lines of a named family are charted under its label, the others under their own Burgers vector.
std::string burgersChartKey(const SegmentProperties& segment){
    return segment.family.name != OTHER_BURGERS_FAMILY.name
        ? segment.family.label
        : burgersVectorLabel(segment.burgersLocal);
}

void addToTally(FamilyTally& tally, double length){
    tally.segmentCount++;
    tally.totalLength += length;
}

LineSummary summarizeLines(const std::vector<DislocationLine>& lines){
    LineSummary summary;
    double minLength = std::numeric_limits<double>::max();
    for(const DislocationLine& line : lines){
        summary.totalLength += line.length;
        summary.totalPoints += line.points.size();
        summary.maxLength = std::max(summary.maxLength, line.length);
        minLength = std::min(minLength, line.length);
        addToTally(summary.byBurgersVector[burgersChartKey(line.segment)], line.length);
        FamilyTally& family = summary.byFamily[line.segment.family.name];
        family.label = line.segment.family.label;
        addToTally(family, line.length);
    }
    summary.minLength = lines.empty() ? 0 : minLength;
    summary.meanLength = lines.empty() ? 0.0 : summary.totalLength / lines.size();
    return summary;
}

void logLineSummary(std::size_t skippedSegments, std::size_t lineCount, const LineSummary& summary, bool clipped){
    if(skippedSegments > 0){
        spdlog::info("  Skipped {} degenerate segments (fewer than two points)", skippedSegments);
    }
    spdlog::info(
        "  {} dislocation lines{} ({} points), length min {:.2f} / mean {:.2f} / max {:.2f}",
        lineCount, clipped ? " after PBC clipping" : "", summary.totalPoints,
        summary.minLength, summary.meanLength, summary.maxLength
    );
    if(summary.byFamily.empty()) return;

    std::string families;
    for(const auto& [_, tally] : summary.byFamily){
        if(!families.empty()) families += ", ";
        families += fmt::format("{} {} ({:.1f}%)", tally.label, tally.segmentCount,
            summary.totalLength > 0 ? 100.0 * tally.totalLength / summary.totalLength : 0.0);
    }
    spdlog::info("  Burgers vector families (lines, % of length): {}", families);
}

void writeLinesTable(const std::string& filePath, const std::vector<DislocationLine>& lines){
    auto vecAsList = [](const Vector3& v){ return std::vector<double>{v.x(), v.y(), v.z()}; };
    auto pointAsList = [](const Point3& p){ return std::vector<double>{p.x(), p.y(), p.z()}; };
    streamLinesToParquet(
        filePath,
        lines.size(),
        [&](std::size_t i, std::vector<Point3>& outPoints){ outPoints = lines[i].points; },
        [&](ColumnarLineWriter& writer, std::size_t i){
            const DislocationLine& line = lines[i];
            const SegmentProperties& segment = line.segment;
            writer.field("length", line.length);
            writer.field("num_points", static_cast<std::int64_t>(line.points.size()));
            writer.field("magnitude", segment.magnitude);
            writer.field("burgers_vector_local", vecAsList(segment.burgersLocal));
            writer.field("burgers_vector_global", vecAsList(segment.burgersGlobal));
            writer.field("crystal_structure", segment.crystalStructure);
            writer.field("burgers_family", segment.family.name);
            writer.field("burgers_family_label", segment.family.label);
            writer.field("cluster_id", static_cast<std::int64_t>(segment.clusterId));
            writer.field("segment_id", static_cast<std::int64_t>(segment.segmentId));
            writer.field("head_vertex", pointAsList(line.points.empty() ? Point3::Origin() : line.points.front()));
            writer.field("tail_vertex", pointAsList(line.points.empty() ? Point3::Origin() : line.points.back()));
        }
    );
}

json burgersVectorChart(const LineSummary& summary){
    json burgersLabels = json::array(), segmentCounts = json::array(), burgersLengths = json::array();
    for(const auto& [label, tally] : summary.byBurgersVector){
        burgersLabels.push_back(label);
        segmentCounts.push_back(tally.segmentCount);
        burgersLengths.push_back(tally.totalLength);
    }
    return {
        {"burgers_counts", {{"burgers_vector", burgersLabels}, {"segment_count", segmentCounts}}},
        {"burgers_lengths", {{"burgers_vector", burgersLabels}, {"total_length", burgersLengths}}}
    };
}

json lineListing(std::size_t lineCount, const LineSummary& summary){
    return {
        {"dislocations", static_cast<int64_t>(lineCount)},
        {"total_points", static_cast<int64_t>(summary.totalPoints)},
        {"average_segment_length", summary.meanLength},
        {"max_segment_length", summary.maxLength},
        {"min_segment_length", summary.minLength},
        {"total_length", summary.totalLength}
    };
}

json burgersFamilyRows(const LineSummary& summary){
    json rows = json::array();
    for(const auto& [familyName, tally] : summary.byFamily){
        rows.push_back({
            {"burgers_family", familyName},
            {"burgers_family_label", tally.label},
            {"segment_count", static_cast<int64_t>(tally.segmentCount)},
            {"total_length", tally.totalLength},
            {"length_fraction", summary.totalLength > 0 ? tally.totalLength / summary.totalLength : 0.0}
        });
    }
    return rows;
}

json structureCountRows(const StructureContext& context){
    const std::size_t atomCount = context.atomCount();
    std::map<int, std::int64_t> structureCounts;
    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        structureCounts[structureTypeOf(context, atomIndex)]++;
    }

    const double atomsForFraction = atomCount > 0 ? static_cast<double>(atomCount) : 1.0;
    json rows = json::array();
    for(const auto& [structureType, count] : structureCounts){
        rows.push_back({
            {"structure_id", structureType},
            {"structure_name", structureTypeName(structureType)},
            {"count", count},
            {"fraction", static_cast<double>(count) / atomsForFraction}
        });
    }
    return rows;
}

json networkListings(
    const DislocationNetwork& network,
    double cellVolume,
    const DislocationsExportOptions& options,
    const LineSummary& summary,
    const StructureContext& structureContext
){
    json listings;
    if(options.exportDislocationNetworkStats)
        listings["network_statistics"] = getNetworkStatistics(network, cellVolume);
    if(options.exportCircuitInformation)
        listings["circuit_information"] = getCircuitInformation(network);
    if(options.exportJunctions)
        listings["junction_information"] = getJunctionInformation(network);
    listings["burgers_families"] = burgersFamilyRows(summary);
    listings["structure_counts"] = structureCountRows(structureContext);
    return listings;
}

}

void streamDislocationsToFile(
    const std::string& linesFilePath,
    const std::string& summaryFilePath,
    const DislocationNetwork& network,
    const SimulationCell& simulationCell,
    const DislocationsExportOptions& options,
    const StructureAnalysis& structureAnalysis
){
    const std::vector<const DislocationSegment*> segments = nonDegenerateSegments(network);
    const std::vector<DislocationLine> lines = collectLines(segments, simulationCell, options.clipPbcSegments);
    const LineSummary summary = summarizeLines(lines);
    logLineSummary(network.segments().size() - segments.size(), lines.size(), summary, options.clipPbcSegments);
    writeLinesTable(linesFilePath, lines);

    json doc;
    doc["export"]["ChartExporter"] = burgersVectorChart(summary);
    doc["main_listing"] = lineListing(lines.size(), summary);
    doc["sub_listings"] = networkListings(
        network, simulationCell.volume3D(), options, summary, structureAnalysis.context()
    );
    JsonUtils::writeJsonToParquet(doc, summaryFilePath);
}

namespace {

// How far a corner may move when its face is unwrapped and still count as the same vertex.
constexpr double kUnwrapTolerance = 1e-6;

struct UnwrappedMesh{
    std::vector<Point3> points;
    std::vector<std::array<int, 3>> faces;
};

// A face that straddles a periodic boundary is unwrapped, and each corner that moves gets a vertex of its own, so
// that every face can be drawn as one triangle. Vertex i of the mesh stays point i.
UnwrappedMesh unwrapFaces(const InterfaceMeshTopology& mesh, const SimulationCell& cell){
    const auto& vertices = mesh.vertices();
    UnwrappedMesh unwrapped;
    unwrapped.points.reserve(vertices.size());
    for(const auto* vertex : vertices){
        unwrapped.points.push_back(vertex->pos());
    }

    unwrapped.faces.reserve(mesh.faces().size());
    for(const auto* face : mesh.faces()){
        if(!face || !face->edges()) continue;
        std::array<int, 3> corners{};
        std::array<Point3, 3> positions{};
        auto* edge = face->edges();
        for(int i = 0; i < 3; ++i, edge = edge->nextFaceEdge()){
            corners[i] = edge->vertex1()->index();
            positions[i] = edge->vertex1()->pos();
        }

        cell.unwrapPositions(positions.data(), 3);

        for(int i = 0; i < 3; ++i){
            if(vertices[corners[i]]->pos().equals(positions[i], kUnwrapTolerance)) continue;
            corners[i] = static_cast<int>(unwrapped.points.size());
            unwrapped.points.push_back(positions[i]);
        }
        unwrapped.faces.push_back(corners);
    }
    return unwrapped;
}

std::string withoutParquetExtension(const std::string& filePath){
    constexpr std::string_view extension = ".parquet";
    return filePath.ends_with(extension) ? filePath.substr(0, filePath.size() - extension.size()) : filePath;
}

json columnJson(const AffineTransformation& matrix, int column){
    return {matrix.column(column).x(), matrix.column(column).y(), matrix.column(column).z()};
}

json cellJson(const SimulationCell& cell){
    const auto& matrix = cell.matrix();
    const auto& pbc = cell.pbcFlags();
    return {
        {"matrix", {columnJson(matrix, 0), columnJson(matrix, 1), columnJson(matrix, 2)}},
        {"origin", columnJson(matrix, 3)},
        {"pbc", {pbc[0], pbc[1], pbc[2]}}
    };
}

// V - E + F of the mesh as it is, counting each undirected edge once.
std::int64_t eulerCharacteristic(const InterfaceMeshTopology& mesh){
    std::unordered_set<uint64_t> edges;
    edges.reserve(mesh.faces().size() * 3);
    for(const auto* face : mesh.faces()){
        if(!face || !face->edges()) continue;
        auto* edge = face->edges();
        do{
            edges.insert(packEdgeKey(edge->vertex1()->index(), edge->vertex2()->index()));
            edge = edge->nextFaceEdge();
        }while(edge != face->edges());
    }
    return static_cast<int64_t>(mesh.vertices().size())
        - static_cast<int64_t>(edges.size())
        + static_cast<int64_t>(mesh.faces().size());
}

}

void streamMeshToFile(
    const std::string& filePath,
    const InterfaceMeshTopology& mesh,
    const StructureAnalysis& structureAnalysis,
    const std::optional<InterfaceMeshFlags>& interfaceFlags
){
    const SimulationCell& cell = structureAnalysis.context().simCell;
    const UnwrappedMesh unwrapped = unwrapFaces(mesh, cell);

    const std::string basePath = withoutParquetExtension(filePath);
    const std::string verticesPath = basePath + ".vertices.parquet";
    const std::string facetsPath = basePath + ".facets.parquet";
    streamMeshTablesToParquet(verticesPath, facetsPath, unwrapped.points, unwrapped.faces);

    json doc;
    doc["main_listing"] = {
        {"total_nodes", static_cast<int64_t>(unwrapped.points.size())},
        {"total_facets", static_cast<int64_t>(unwrapped.faces.size())}
    };
    doc["export"]["MeshExporter"][kParquetSourceKey] = {
        {"vertices", std::filesystem::path(verticesPath).filename().string()},
        {"facets", std::filesystem::path(facetsPath).filename().string()},
        {"cell", cellJson(cell)}
    };
    doc["topology"] = {
        {"euler_characteristic", eulerCharacteristic(mesh)}
    };
    if(interfaceFlags){
        doc["topology"]["is_completely_good"] = interfaceFlags->isCompletelyGood;
        doc["topology"]["is_completely_bad"] = interfaceFlags->isCompletelyBad;
    }

    JsonUtils::writeJsonToParquet(doc, filePath);
}

}
