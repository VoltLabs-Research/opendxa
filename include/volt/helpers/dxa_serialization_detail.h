#pragma once

#include <volt/helpers/dxa_serialization.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

// What the serializer builds its files from. LineReconstructionDXA writes its lines and its interface mesh with some
// of it too: polylineLength, getGlobalBurgersVector and unwrapFaces.
namespace Volt::DxaSerialization::Detail {

inline constexpr const char* kParquetSourceKey = "__parquet_source__";

Vector3 getGlobalBurgersVector(const ClusterVector& burgersVector);
std::string burgersVectorLabel(const Vector3& vector);

struct CubicBurgersFamily{
    std::array<double, 3> sortedAbsComponents;
    BurgersFamily family;
};

struct HexagonalBurgersFamily{
    double basalLength;
    double axialLength;
    BurgersFamily family;
};

// HCP lengths in the units of the local Burgers vectors, where the nearest-neighbour distance is a = 1/√2 and
// c/a takes its ideal value √(8/3).
inline constexpr double kHcpA = 0.7071067811865476;         // 1/√2
inline constexpr double kHcpShockley = 0.4082482904638630;  // a/√3
inline constexpr double kHcpC = 1.1547005383792515;         // a·√(8/3)
inline constexpr double kHcpHalfC = 0.5773502691896257;     // c/2

extern const BurgersFamily OTHER_BURGERS_FAMILY;
extern const CubicBurgersFamily FCC_BURGERS_FAMILIES[];
extern const CubicBurgersFamily BCC_BURGERS_FAMILIES[];
extern const CubicBurgersFamily SC_BURGERS_FAMILIES[];
extern const HexagonalBurgersFamily HCP_BURGERS_FAMILIES[];

template<std::size_t N>
const BurgersFamily& matchCubicBurgersFamily(
    const Vector3& localBurgers,
    const CubicBurgersFamily (&families)[N],
    double tolerance
);
const BurgersFamily& matchHexagonalBurgersFamily(const Vector3& localBurgers, double tolerance);

int structureTypeOf(const StructureContext& context, std::size_t atomIndex);
std::array<std::uint64_t, 3> canonicalFaceKey(const std::array<DelaunayTessellation::VertexHandle, 3>& corners);

// Past this many periodic crossings within one segment, the rest of it is kept as it is.
inline constexpr int kMaxCrossingsPerSegment = 10;
inline constexpr double kCrossingEpsilon = 1e-9;

std::array<const DislocationNode*, 2> endNodes(const DislocationSegment& segment);
int countDanglingSegments(const DislocationNetwork& network);
json getNetworkStatistics(const DislocationNetwork& network, double cellVolume);
json getJunctionInformation(const DislocationNetwork& network);
json getCircuitInformation(const DislocationNetwork& network);

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
    // The degenerate segments, which have no line.
    std::size_t skippedSegments = 0;
    std::size_t lineCount = 0;
    // Whether the lines were split where they cross a periodic face.
    bool clipped = false;
    double totalLength = 0;
    int totalPoints = 0;
    double minLength = 0;
    double meanLength = 0;
    double maxLength = 0;
    std::map<std::string, FamilyTally> byBurgersVector;
    std::map<std::string, FamilyTally> byFamily;
};

std::vector<const DislocationSegment*> nonDegenerateSegments(const DislocationNetwork& network);
std::string crystalStructureOf(const Cluster* cluster);
SegmentProperties segmentProperties(const DislocationSegment& segment);
// Inline so that collectLines inlines it: an out-of-line copy rounds some lengths one ulp differently.
inline double polylineLength(const std::vector<Point3>& points){
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
);
std::string burgersChartKey(const SegmentProperties& segment);
void addToTally(FamilyTally& tally, double length);
LineSummary summarizeLines(const std::vector<DislocationLine>& lines, std::size_t skippedSegments, bool clipped);
void writeLinesTable(const std::string& filePath, const std::vector<DislocationLine>& lines);
json burgersVectorChart(const LineSummary& summary);
json lineListing(std::size_t lineCount, const LineSummary& summary);
json burgersFamilyRows(const LineSummary& summary);
json structureCountRows(const StructureContext& context);
json networkListings(
    const DislocationNetwork& network,
    double cellVolume,
    const DislocationsExportOptions& options,
    const LineSummary& summary,
    const StructureContext& structureContext
);

// How far a corner may move when its face is unwrapped and still count as the same vertex.
inline constexpr double kUnwrapTolerance = 1e-6;

struct UnwrappedMesh{
    std::vector<Point3> points;
    std::vector<std::array<int, 3>> faces;
};

// A face that straddles a periodic boundary is unwrapped, and each corner that moves gets a vertex of its own, so
// that every face can be drawn as one triangle. Vertex i of the mesh stays point i. Any triangle mesh of
// HalfEdgeMesh will do.
template<typename Mesh>
UnwrappedMesh unwrapFaces(const Mesh& mesh, const SimulationCell& cell){
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

std::string withoutParquetExtension(const std::string& filePath);
json columnJson(const AffineTransformation& matrix, int column);
json cellJson(const SimulationCell& cell);
std::int64_t eulerCharacteristic(const InterfaceMeshTopology& mesh);
json writeMeshTables(const std::string& filePath, const InterfaceMeshTopology& mesh, const SimulationCell& cell);

}
