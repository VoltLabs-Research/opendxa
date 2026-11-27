#pragma once

#include <volt/helpers/dxa_serialization.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/utilities/columnar_row_writer.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <span>
#include <string>
#include <vector>

// What the serializer builds its files from. LineReconstructionDXA writes its lines and its interface mesh with some
// of it too: polylineLength, getGlobalBurgersVector and unwrapFaces.
namespace Volt::DxaSerialization::Detail {

// Past this many periodic crossings within one segment, the rest of it is kept as it is.
inline constexpr int kMaxCrossingsPerSegment = 10;
inline constexpr double kCrossingEpsilon = 1e-9;

// How far a corner may move when its face is unwrapped and still count as the same vertex.
inline constexpr double kUnwrapTolerance = 1e-6;

// How far each number of a Burgers vector's signature may be from a family's, in the units of the local vectors.
inline constexpr double kBurgersFamilyTolerance = 0.01;

inline const BurgersFamily OTHER_BURGERS_FAMILY{"Other", "Other"};

// A family and the signature of its vectors (see burgersSignature).
struct BurgersFamilyMatch{
    std::array<double, 3> signature;
    BurgersFamily family;
};

// The sorted absolute components on a cubic lattice, the basal and axial lengths on a hexagonal one.
std::array<double, 3> burgersSignature(const Vector3& localBurgers, LatticeStructureType structure);
std::span<const BurgersFamilyMatch> burgersFamilies(LatticeStructureType structure);

Vector3 getGlobalBurgersVector(const ClusterVector& burgersVector);

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

    void add(double length){
        segmentCount++;
        totalLength += length;
    }
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
SegmentProperties segmentProperties(
    const ClusterVector& burgersVector,
    std::string crystalStructure,
    int clusterId,
    int segmentId
);

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
LineSummary summarizeLines(const std::vector<DislocationLine>& lines, std::size_t skippedSegments, bool clipped);
// The columns every dislocation line table starts with, LineReconstructionDXA's too: the line's length and point
// count, then its Burgers vector, crystal, family and cluster.
void writeBurgersColumns(ColumnarRowWriter& writer, const DislocationLine& line);
void writeLinesTable(const std::string& filePath, const std::vector<DislocationLine>& lines);

json lineListing(const LineSummary& summary);
json burgersVectorChart(const LineSummary& summary);
json burgersFamilyRows(const LineSummary& summary);
std::map<int, int> structureCounts(const StructureContext& context);
json networkStatisticsRows(
    const DislocationNetwork& network,
    const std::vector<const DislocationSegment*>& segments,
    double cellVolume
);
json circuitRows(const DislocationNetwork& network);
json junctionRows(const DislocationNetwork& network);

int structureTypeOf(const StructureContext& context, std::size_t atomIndex);

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

std::int64_t eulerCharacteristic(const InterfaceMeshTopology& mesh);
// Writes the unwrapped mesh's vertex and facet tables next to filePath, and returns the document that points to them.
json meshDocument(const std::string& filePath, const InterfaceMeshTopology& mesh, const SimulationCell& cell);

}
