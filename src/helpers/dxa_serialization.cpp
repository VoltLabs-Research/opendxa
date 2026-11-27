#include <volt/helpers/dxa_serialization.h>
#include <volt/helpers/dxa_serialization_detail.h>
#include <volt/helpers/dxa_serialization_log.h>
#include <volt/helpers/burgers_circuit.h>
#include <volt/helpers/edge_key.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/utilities/json_utils.h>
#include <volt/utilities/parquet_line_writer.h>
#include <volt/utilities/parquet_mesh_writer.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <numeric>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace Volt::DxaSerialization {

using namespace Detail;

BurgersFamily classifyBurgersFamily(const Vector3& localBurgers, const std::string& crystalStructure){
    LatticeStructureType structure = LATTICE_OTHER;
    if(!parseLatticeStructureType(crystalStructure, structure)) return OTHER_BURGERS_FAMILY;

    const std::array<double, 3> signature = burgersSignature(localBurgers, structure);
    const auto close = [](double a, double b){ return std::abs(a - b) < kBurgersFamilyTolerance; };
    for(const BurgersFamilyMatch& candidate : burgersFamilies(structure)){
        if(std::ranges::equal(signature, candidate.signature, close)) return candidate.family;
    }
    return OTHER_BURGERS_FAMILY;
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

void streamDislocationsToFile(
    const std::string& linesFilePath,
    const std::string& summaryFilePath,
    const DislocationNetwork& network,
    const StructureAnalysis& structureAnalysis,
    const DislocationsExportOptions& options
){
    DxaSerializationLog log;
    const StructureContext& context = structureAnalysis.context();
    const std::vector<const DislocationSegment*> segments = nonDegenerateSegments(network);
    const std::vector<DislocationLine> lines = collectLines(segments, context.simCell, options.clipPbcSegments);
    const LineSummary summary = summarizeLines(
        lines, network.segments().size() - segments.size(), options.clipPbcSegments
    );
    log.dump(summary);
    writeLinesTable(linesFilePath, lines);

    json doc;
    doc["main_listing"] = lineListing(summary);
    doc["export"]["ChartExporter"] = burgersVectorChart(summary);
    json& listings = doc["sub_listings"];
    listings["burgers_families"] = burgersFamilyRows(summary);
    listings["structure_counts"] = structureCountRows(context);
    if(options.exportDislocationNetworkStats)
        listings["network_statistics"] = networkStatisticsRows(network, segments, context.simCell.volume3D());
    if(options.exportCircuitInformation)
        listings["circuit_information"] = circuitRows(network);
    if(options.exportJunctions)
        listings["junction_information"] = junctionRows(network);
    JsonUtils::writeJsonToParquet(doc, summaryFilePath);
}

void streamMeshToFile(const std::string& filePath, const InterfaceMeshTopology& mesh, const SimulationCell& cell){
    JsonUtils::writeJsonToParquet(writeMeshTables(filePath, mesh, cell), filePath);
}

void streamInterfaceMeshToFile(const std::string& filePath, const InterfaceMesh& interfaceMesh){
    json doc = writeMeshTables(filePath, interfaceMesh, interfaceMesh.structureAnalysis().context().simCell);
    doc["topology"]["is_completely_good"] = interfaceMesh.isCompletelyGood();
    doc["topology"]["is_completely_bad"] = interfaceMesh.isCompletelyBad();
    JsonUtils::writeJsonToParquet(doc, filePath);
}

void streamDelaunayTessellationToFile(const std::string& filePath, const DelaunayTessellation& tessellation){
    using CellHandle = DelaunayTessellation::CellHandle;
    const auto isPrimary = [&](CellHandle cell){
        return tessellation.isValidCell(cell) && !tessellation.isGhostCell(cell);
    };

    std::unordered_map<DelaunayTessellation::VertexHandle, int> pointOfVertex;
    json points = json::array();
    json facets = json::array();
    int primaryTetrahedra = 0, boundaryFacets = 0, internalFacets = 0;
    for(CellHandle cell : tessellation.cells()){
        if(!isPrimary(cell)) continue;
        ++primaryTetrahedra;
        for(int face = 0; face < 4; ++face){
            const CellHandle neighbor = tessellation.mirrorFacet(cell, face).first;
            const bool isBoundary = !isPrimary(neighbor);
            // A facet between two primary cells is written once, by the first of them.
            if(!isBoundary && neighbor < cell) continue;
            isBoundary ? ++boundaryFacets : ++internalFacets;

            json corners = json::array();
            for(int corner = 0; corner < 3; ++corner){
                const auto vertex = tessellation.cellVertex(cell, DelaunayTessellation::cellFacetVertexIndex(face, corner));
                const auto [it, inserted] = pointOfVertex.emplace(vertex, static_cast<int>(points.size()));
                if(inserted){
                    const Point3 position = tessellation.vertexPosition(vertex);
                    points.push_back({
                        {"index", static_cast<std::int64_t>(it->second)},
                        {"position", {position.x(), position.y(), position.z()}}
                    });
                }
                corners.push_back(it->second);
            }
            facets.push_back({{"vertices", std::move(corners)}});
        }
    }

    json doc;
    doc["main_listing"] = {
        {"total_primary_tetrahedra", primaryTetrahedra},
        {"total_nodes", static_cast<std::int64_t>(points.size())},
        {"total_facets", static_cast<std::int64_t>(facets.size())},
        {"boundary_facets", boundaryFacets},
        {"internal_facets", internalFacets}
    };
    doc["sub_listings"] = {{"points", points}, {"facets", facets}};
    doc["export"]["MeshExporter"] = {{"vertices", std::move(points)}, {"facets", std::move(facets)}};
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
    json regions = json::array();
    json atomisticExporter = json::object();
    for(const auto& [clusterId, atomIndices] : atomIndicesByCluster){
        largestRegionSize = std::max(largestRegionSize, static_cast<int>(atomIndices.size()));
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
        {"assigned_atoms", static_cast<int>(context.atomCount()) - unassignedAtoms},
        {"unassigned_atoms", unassignedAtoms},
        {"largest_region_size", largestRegionSize}
    };
    doc["sub_listings"] = {{"coherent_crystalline_regions", regions}};
    doc["export"]["AtomisticExporter"] = std::move(atomisticExporter);
    JsonUtils::writeJsonToParquet(doc, filePath);
}

namespace Detail {

std::array<double, 3> burgersSignature(const Vector3& localBurgers, LatticeStructureType structure){
    if(structure == LATTICE_HCP || structure == LATTICE_HEX_DIAMOND){
        return {std::hypot(localBurgers.x(), localBurgers.y()), std::abs(localBurgers.z()), 0.0};
    }
    std::array<double, 3> components{std::abs(localBurgers.x()), std::abs(localBurgers.y()), std::abs(localBurgers.z())};
    std::ranges::sort(components, std::greater<>());
    return components;
}

std::span<const BurgersFamilyMatch> burgersFamilies(LatticeStructureType structure){
    static const BurgersFamilyMatch fcc[] = {
        {{0.5, 0.5, 0.0}, {"1/2<110>", "1/2<110> (Perfect)"}},
        {{1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0}, {"1/6<112>", "1/6<112> (Shockley)"}},
        {{1.0 / 6.0, 1.0 / 6.0, 0.0}, {"1/6<110>", "1/6<110> (Stair-rod)"}},
        {{1.0 / 3.0, 0.0, 0.0}, {"1/3<100>", "1/3<100> (Hirth)"}},
        {{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}, {"1/3<111>", "1/3<111> (Frank)"}}
    };
    static const BurgersFamilyMatch bcc[] = {
        {{0.5, 0.5, 0.5}, {"1/2<111>", "1/2<111> (Perfect)"}},
        {{1.0, 0.0, 0.0}, {"<100>", "<100>"}},
        {{1.0, 1.0, 0.0}, {"<110>", "<110>"}}
    };
    static const BurgersFamilyMatch sc[] = {
        {{1.0, 0.0, 0.0}, {"<100>", "<100>"}},
        {{1.0, 1.0, 0.0}, {"<110>", "<110>"}},
        {{1.0, 1.0, 1.0}, {"<111>", "<111>"}}
    };
    constexpr double a = 0.7071067811865476;
    constexpr double shockley = 0.4082482904638630; 
    constexpr double c = 1.1547005383792515;
    static const BurgersFamilyMatch hcp[] = {
        {{a, 0.0, 0.0}, {"1/3<1-210>", "1/3<1-210> (Perfect basal)"}},
        {{shockley, 0.0, 0.0}, {"1/3<1-100>", "1/3<1-100> (Shockley)"}},
        {{0.0, c, 0.0}, {"<0001>", "<0001> (Perfect c)"}},
        {{0.0, c / 2, 0.0}, {"1/2<0001>", "1/2<0001> (Partial c)"}},
        {{a, c, 0.0}, {"1/3<1-213>", "1/3<1-213> (Perfect c+a)"}}
    };
    switch(structure){
        case LATTICE_FCC:
        case LATTICE_CUBIC_DIAMOND:
            return fcc;
        case LATTICE_BCC:
            return bcc;
        case LATTICE_SC:
            return sc;
        case LATTICE_HCP:
        case LATTICE_HEX_DIAMOND:
            return hcp;
        default:
            return {};
    }
}

Vector3 getGlobalBurgersVector(const ClusterVector& burgersVector){
    if(burgersVector.cluster() == nullptr){
        return burgersVector.localVec();
    }
    return burgersVector.toSpatialVector();
}

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
    if(segment.family.name != OTHER_BURGERS_FAMILY.name) return segment.family.label;
    const Vector3& b = segment.burgersLocal;
    return "[" + std::to_string(b.x()) + ", " + std::to_string(b.y()) + ", " + std::to_string(b.z()) + "]";
}

LineSummary summarizeLines(const std::vector<DislocationLine>& lines, std::size_t skippedSegments, bool clipped){
    LineSummary summary{.skippedSegments = skippedSegments, .lineCount = lines.size(), .clipped = clipped};
    double minLength = std::numeric_limits<double>::max();
    for(const DislocationLine& line : lines){
        summary.totalLength += line.length;
        summary.totalPoints += line.points.size();
        summary.maxLength = std::max(summary.maxLength, line.length);
        minLength = std::min(minLength, line.length);
        summary.byBurgersVector[burgersChartKey(line.segment)].add(line.length);
        FamilyTally& family = summary.byFamily[line.segment.family.name];
        family.label = line.segment.family.label;
        family.add(line.length);
    }
    summary.minLength = lines.empty() ? 0 : minLength;
    summary.meanLength = lines.empty() ? 0.0 : summary.totalLength / lines.size();
    return summary;
}

void writeLinesTable(const std::string& filePath, const std::vector<DislocationLine>& lines){
    const auto xyz = [](const auto& v){ return std::vector<double>{v.x(), v.y(), v.z()}; };
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
            writer.field("burgers_vector_local", xyz(segment.burgersLocal));
            writer.field("burgers_vector_global", xyz(segment.burgersGlobal));
            writer.field("crystal_structure", segment.crystalStructure);
            writer.field("burgers_family", segment.family.name);
            writer.field("burgers_family_label", segment.family.label);
            writer.field("cluster_id", static_cast<std::int64_t>(segment.clusterId));
            writer.field("segment_id", static_cast<std::int64_t>(segment.segmentId));
            writer.field("head_vertex", xyz(line.points.front()));
            writer.field("tail_vertex", xyz(line.points.back()));
        }
    );
}

json lineListing(const LineSummary& summary){
    return {
        {"dislocations", static_cast<int64_t>(summary.lineCount)},
        {"total_points", static_cast<int64_t>(summary.totalPoints)},
        {"average_segment_length", summary.meanLength},
        {"max_segment_length", summary.maxLength},
        {"min_segment_length", summary.minLength},
        {"total_length", summary.totalLength}
    };
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

json networkStatisticsRows(
    const DislocationNetwork& network,
    const std::vector<const DislocationSegment*>& segments,
    double cellVolume
){
    // Summed in segment order, so the total is the same with any thread count.
    double totalLength = 0.0;
    for(const DislocationSegment* segment : segments){
        totalLength += segment->calculateLength();
    }
    const int segmentCount = static_cast<int>(segments.size());
    const auto danglingSegments = std::ranges::count_if(network.segments(), [](const DislocationSegment* segment){
        return segment && (segment->forwardNode().isDangling() || segment->backwardNode().isDangling());
    });

    return json::array({json{
        {"total_network_length", totalLength},
        {"segment_count", segmentCount},
        {"junction_count", static_cast<int>(network.junctionArmCounts().size())},
        {"dangling_segments", static_cast<int>(danglingSegments)},
        {"average_segment_length", segmentCount > 0 ? totalLength / segmentCount : 0.0},
        {"density", cellVolume > 0 ? totalLength / cellVolume : 0.0},
        {"total_segments_including_degenerate", static_cast<int>(network.segments().size())}
    }});
}

json circuitRows(const DislocationNetwork& network){
    std::vector<int> edgeCounts;
    int danglingCircuits = 0;
    int blockedCircuits = 0;
    for(const DislocationSegment* segment : network.segments()){
        if(!segment) continue;
        for(const DislocationNode* node : {&segment->forwardNode(), &segment->backwardNode()}){
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

json junctionRows(const DislocationNetwork& network){
    const std::vector<int> junctionArms = network.junctionArmCounts();
    std::map<int, int> junctionsByArms;
    for(const int arms : junctionArms){
        junctionsByArms[arms]++;
    }
    // A network without junctions gets one row, of junctions with no arms.
    if(junctionsByArms.empty()) junctionsByArms[0] = 0;

    json rows = json::array();
    for(const auto& [arms, count] : junctionsByArms){
        rows.push_back({
            {"junction_arms", arms},
            {"junction_count", count},
            {"total_junctions", static_cast<int>(junctionArms.size())}
        });
    }
    return rows;
}

int structureTypeOf(const StructureContext& context, std::size_t atomIndex){
    return context.structureTypes
        ? context.structureTypes->getInt(atomIndex)
        : static_cast<int>(StructureType::OTHER);
}

json cellJson(const SimulationCell& cell){
    const auto column = [&](int i){
        const Vector3 v = cell.matrix().column(i);
        return json{v.x(), v.y(), v.z()};
    };
    const auto& pbc = cell.pbcFlags();
    return {
        {"matrix", {column(0), column(1), column(2)}},
        {"origin", column(3)},
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

json writeMeshTables(const std::string& filePath, const InterfaceMeshTopology& mesh, const SimulationCell& cell){
    const UnwrappedMesh unwrapped = unwrapFaces(mesh, cell);

    constexpr std::string_view extension = ".parquet";
    const std::string basePath = filePath.ends_with(extension)
        ? filePath.substr(0, filePath.size() - extension.size())
        : filePath;
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
    return doc;
}

}

}
