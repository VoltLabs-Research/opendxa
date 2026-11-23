#pragma once

#include <volt/core/volt.h>
#include <volt/utilities/memory_pool.h>
#include <volt/helpers/dislocation_network.h>
#include <volt/pipeline/defect_cell_index.h>
#include <volt/pipeline/interface_mesh.h>
#include <string_view>
#include <optional>

namespace Volt{

// A vertex the breadth-first search for Burgers circuits reached, with the lattice coordinate and the frame rotation
// accumulated along the path from the search's start vertex.
struct BurgersLoopSearchNode{
	InterfaceMesh::Vertex* vertex;
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
		_unusedCircuit(nullptr),
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
	void indexDefectCells();
	BurgersCircuit* allocateCircuit();
	BurgersCircuit* buildReverseCircuit(BurgersCircuit* forwardCircuit);

	void discardCircuit(BurgersCircuit* circuit);
	void createAndTraceSegment(const ClusterVector& burgersVector, BurgersCircuit* forwardCircuit, int maxCircuitLength);
	void traceSegment(DislocationNode& node, int maxCircuitLength, bool markSwept);
	InterfaceMesh::Edge* nextStartEdge(const BurgersCircuit& circuit);
	void sweepCircuit(DislocationNode& node, bool markSwept);
	bool lengthenCircuit(DislocationNode& node, bool markSwept);
	void appendLinePoint(DislocationNode& node);
	void markCoreCells(DislocationNode& node, const Point3& center);
	void createSecondarySegment(InterfaceMesh::Edge* firstEdge, BurgersCircuit* outerCircuit, int maxCircuitLength);

	void findPrimarySegments(int maxCircuitLength);
	bool createPrimaryCircuit(InterfaceMesh::Edge* edge, SearchNode* from, SearchNode* reached, int maxCircuitLength);
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
	BurgersCircuit* _unusedCircuit;
	size_t _edgeStartIndex = 0;

	bool _markCoreAtoms = false;
	std::optional<DefectCellIndex> _defectCells;
	std::vector<CoreCellClaim> _coreCellClaims;
};

}
