#pragma once

#include <volt/core/volt.h>
#include <volt/utilities/memory_pool.h>
#include <volt/helpers/dislocation_network.h>
#include <volt/helpers/defect_cell_index.h>
#include <volt/pipeline/interface_mesh.h>
#include <tbb/enumerable_thread_specific.h>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <string_view>
#include <optional>
#include <vector>

namespace Volt{

// A vertex the breadth-first search for Burgers circuits reached, with the lattice coordinate and the frame rotation
// accumulated along the path from the search's start vertex.
struct BurgersLoopSearchNode{
	InterfaceMesh::Vertex* vertex;
	int vertexIndex;
	Point3 coord;
	Matrix3 tm;
	int depth;
	InterfaceMesh::Edge* viaEdge;
	BurgersLoopSearchNode* parent;
};

class BurgersLoopBuilder{
public:
	using SearchNode = BurgersLoopSearchNode;

	struct CoreCellClaim{
		DislocationNode* node = nullptr;
		bool claimedWhileCapStored = false;
	};

	BurgersLoopBuilder(InterfaceMesh& mesh, int maxTrialCircuitSize, int maxCircuitElongation) :
		_mesh(mesh),
		_maxTrialCircuitSize(maxTrialCircuitSize),
		_maxExtendedCircuitSize(maxTrialCircuitSize + maxCircuitElongation){}

	const InterfaceMesh& mesh() const{
		return _mesh;
	}

	DislocationNetwork& network(){
		return _network;
	}

	const DislocationNetwork& network() const{
		return _network;
	}

	const SimulationCell& cell() const{
		return mesh().structureAnalysis().context().simCell;
	}

	void traceDislocationSegments();
	void finishDislocationSegments(std::string_view referenceTopologyName);

	const std::vector<DislocationNode*>& danglingNodes() const{
		return _danglingNodes;
	}

	void setMarkCoreAtoms(bool markCoreAtoms){
		_markCoreAtoms = markCoreAtoms;
	}

	[[nodiscard]] int maxTrialCircuitSize() const noexcept{
		return _maxTrialCircuitSize;
	}

	[[nodiscard]] int maxExtendedCircuitSize() const noexcept{
		return _maxExtendedCircuitSize;
	}

	[[nodiscard]] std::vector<int> coreAtomDislocationIds(std::size_t particleCount) const;

private:
	// Bounds the walk around a vertex, should its fan of facets not close.
	static constexpr int kMaxFanSteps = 1000;

	// An edge the search may walk, with what the search reads of it.
	struct FreeEdge{
		Vector3 clusterVector;
		ClusterTransition* transition;
		InterfaceMesh::Edge* edge;
		InterfaceMesh::Vertex* target;
		int targetIndex;
	};
	struct LiveFreeEdges;
	class FrozenFreeEdges;
	struct CandidateLinks;
	class ChangedRegion;
	class RingSets;

	// Breadth-first search over free edges, at most maxDepth edges away from its start vertex. One search object is
	// reused for many searches and allocates nothing once warm: its nodes live in blocks kept between searches, and the
	// vertices it reached sit in a small open-addressing table whose slots carry the generation of the search that
	// filled them, so starting a search clears nothing. Each thread keeps its own.
	class CircuitSearch{
	public:
		// Calls onClosure(edge, from, reached) for each edge that closes a circuit with a non-zero Burgers vector, and
		// stops when it returns true.
		// Edges is LiveFreeEdges or FrozenFreeEdges: where the search reads the free edges leaving a vertex.
		template<typename Edges, typename OnClosure>
		void run(const Edges& freeEdges, InterfaceMesh::Vertex* start, int maxDepth, OnClosure&& onClosure);

		// Scratch space for the circuit a closure builds.
		std::vector<InterfaceMesh::Edge*>& candidate() noexcept{
			return _candidate;
		}

	private:
		static constexpr std::size_t kBlockSize = 1024;
		static constexpr std::size_t kInitialSlots = 1024;

		struct Slot{
			std::uint32_t generation = 0;
			int vertex = -1;
			SearchNode* node = nullptr;
		};

		SearchNode* newNode();
		void startSearch();
		SearchNode* reached(int vertex) const;
		void enqueue(SearchNode* node);
		void insertSlot(int vertex, SearchNode* node);
		void growSlots();

		std::size_t slotOf(int vertex) const{
			return (static_cast<std::uint32_t>(vertex) * 0x9E3779B1u) & (_slots.size() - 1);
		}

		std::vector<std::unique_ptr<SearchNode[]>> _blocks;
		std::size_t _usedNodes = 0;
		std::vector<SearchNode*> _queue;
		std::vector<Slot> _slots = std::vector<Slot>(kInitialSlots);
		std::uint32_t _generation = 0;
		std::vector<InterfaceMesh::Edge*> _candidate;
	};

	// What a search does with a circuit it closed.
	enum class PrimaryVerdict : std::uint8_t{
		KeepSearching,
		Stops,
		StartsSegment
	};

	void sweepFacet(InterfaceMesh::Face* facet, BurgersCircuit* circuit, bool markSwept);
	static void storeCaps(const std::vector<DislocationNode*>& nodes);
	static void trimPreliminaryPoints(DislocationSegment& segment);
	static void expressInReferenceStructure(DislocationSegment& segment, std::string_view referenceTopologyName);
	static void orientCanonically(DislocationSegment& segment);
	static bool isFreeEdge(const InterfaceMesh::Edge* edge);
	static bool hasFreeEdge(const InterfaceMesh::Vertex* vertex);
	static void initRoot(SearchNode* root, InterfaceMesh::Vertex* vertex);
	static void initChild(SearchNode* child, SearchNode* parent, const FreeEdge& edge, const Point3& coord);
	static bool closesNonZeroCircuit(const SearchNode* from, const FreeEdge& edge, const Point3& coord, const SearchNode* reached);
	std::vector<InterfaceMesh::Edge*> findClosingEdges(const FrozenFreeEdges& freeEdges, int searchDepth);
	CircuitSearch& threadSearch();
	static void collectPrimaryCircuit(InterfaceMesh::Edge* edge, const SearchNode* from, const SearchNode* reached, std::vector<InterfaceMesh::Edge*>& circuit);
	static bool leavesOtherCircuit(InterfaceMesh::Edge* aIn, InterfaceMesh::Edge* aOut, InterfaceMesh::Edge* bIn, InterfaceMesh::Edge* bOut);
	static bool crossedAtVertex(InterfaceMesh::Edge* edge1, InterfaceMesh::Edge* edge2, const CandidateLinks& links);
	static bool intersectsOtherCircuits(std::span<InterfaceMesh::Edge* const> candidate);
	static PrimaryVerdict judgePrimaryCircuit(
		InterfaceMesh::Edge* edge,
		const SearchNode* from,
		const SearchNode* reached,
		std::vector<InterfaceMesh::Edge*>& circuit,
		Vector3* burgersVector = nullptr
	);

	void indexDefectCells();
	BurgersCircuit* buildReverseCircuit(BurgersCircuit* forwardCircuit);

	void createAndTraceSegment(const ClusterVector& burgersVector, BurgersCircuit* forwardCircuit, int maxCircuitLength);
	void traceSegment(DislocationNode& node, int maxCircuitLength, bool markSwept);
	InterfaceMesh::Edge* nextStartEdge(const BurgersCircuit& circuit);
	void sweepCircuit(DislocationNode& node, bool markSwept);
	bool lengthenCircuit(DislocationNode& node, bool markSwept);
	void appendLinePoint(DislocationNode& node);
	void markCoreCells(DislocationNode& node, const Point3& center);
	[[nodiscard]] std::optional<std::size_t> claimSlot(DelaunayTessellation::CellHandle cell) const;
	void createSecondarySegment(InterfaceMesh::Edge* firstEdge, BurgersCircuit* outerCircuit, int maxCircuitLength);

	void findPrimarySegments(int maxCircuitLength);
	bool createPrimaryCircuit(InterfaceMesh::Edge* edge, SearchNode* from, SearchNode* reached, int maxCircuitLength, std::vector<InterfaceMesh::Edge*>& circuit);
	bool tryRemoveTwoCircuitEdges(InterfaceMesh::Edge*& edge0, InterfaceMesh::Edge*& edge1, InterfaceMesh::Edge*& edge2);
	bool tryRemoveThreeCircuitEdges(InterfaceMesh::Edge*& edge0, InterfaceMesh::Edge*& edge1, InterfaceMesh::Edge*& edge2, bool markSwept);
	bool tryRemoveOneCircuitEdge(InterfaceMesh::Edge*& edge0, InterfaceMesh::Edge*& edge1, InterfaceMesh::Edge*& edge2, bool markSwept);
	bool trySweepTwoFacets(InterfaceMesh::Edge*& edge0, InterfaceMesh::Edge*& edge1, InterfaceMesh::Edge*& edge2, bool markSwept);
	bool tryInsertOneCircuitEdge(InterfaceMesh::Edge*& edge0, InterfaceMesh::Edge*& edge1, bool markSwept);

	void joinSegments(int maxCircuitLength);
	void createSecondarySegments(int maxCircuitLength);
	void linkAdjacentCircuits();
	void joinArms(DislocationNode* node);
	void formJunction(DislocationNode* node, const Point3& centerOfMass);
	void closeLoop(DislocationNode* node1, DislocationNode* node2);
	void mergeSegments(DislocationNode* node1, DislocationNode* node2);

	Vector3 calculateShiftVector(const Point3& a, const Point3& b) const{
		Vector3 d = cell().absoluteToReduced(b - a);
		d.x() = cell().pbcFlags()[0] ? floor(d.x() + double(0.5)) : double(0);
		d.y() = cell().pbcFlags()[1] ? floor(d.y() + double(0.5)) : double(0);
		d.z() = cell().pbcFlags()[2] ? floor(d.z() + double(0.5)) : double(0);
		return cell().reducedToAbsolute(d);
	}

	InterfaceMesh& _mesh;
	DislocationNetwork _network;

	int _maxTrialCircuitSize;
	int _maxExtendedCircuitSize;

	MemoryPool<BurgersCircuit> _circuitPool;
	std::vector<DislocationNode*> _danglingNodes;
	size_t _edgeStartIndex = 0;

	// While the serial primary pass creates a segment, the vertices whose edges or facets it changes.
	std::vector<InterfaceMesh::Vertex*>* _changedVertices = nullptr;

	bool _markCoreAtoms = false;
	std::optional<DefectCellIndex> _defectCells;
	std::vector<CoreCellClaim> _coreCellClaims;

	// One circuit search per thread, made on first use and kept for the whole trace.
	tbb::enumerable_thread_specific<std::unique_ptr<CircuitSearch>> _searches;
};

}
