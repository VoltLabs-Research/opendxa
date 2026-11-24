#pragma once

#include <volt/core/volt.h>
#include <volt/core/simulation_cell.h>
#include <volt/pipeline/interface_mesh.h>
#include <volt/helpers/burgers_circuit.h>
#include <volt/helpers/cluster_vector.h>
#include <vector>
#include <algorithm>
#include <ranges>
#include <cmath>

namespace Volt{

struct DislocationSegment;

struct DislocationNode{
	DislocationSegment* segment = nullptr;
	DislocationNode* oppositeNode = nullptr;
	DislocationNode* junctionRing = nullptr;
	BurgersCircuit* circuit = nullptr;
	
	DislocationNode() = default;

	[[nodiscard]] inline ClusterVector burgersVector() const;
	[[nodiscard]] inline const Point3& position() const;
	[[nodiscard]] inline bool isForwardNode() const;
	[[nodiscard]] inline bool isBackwardNode() const;
	
	[[nodiscard]] inline bool isDangling() const{
		return junctionRing == this;
	}

	bool formsJunctionWith(DislocationNode* other) const{
		for(auto *node = junctionRing; ; node = node->junctionRing){
			if(node == other) return true;
			if(node == this) break;
		}

		return false;
	}

	void connectNodes(DislocationNode* other){
		assert(!formsJunctionWith(other));
		assert(!other->formsJunctionWith(this));
		std::swap(junctionRing, other->junctionRing);
		assert(formsJunctionWith(other));
		assert(other->formsJunctionWith(this));
	}

	void dissolveJunction(){
		for(auto *node = junctionRing; node != this;){
			auto *next = node->junctionRing;
			node->junctionRing = node;
			node = next;
		}

		junctionRing = this;
	}

	int countJunctionArms() const{
		int count = 1;
		for(auto *node = junctionRing; node != this; node = node->junctionRing){
			++count;
		}

		return count;
	}
};

struct DislocationSegment{
	int id;
	std::vector<Point3> line;
	std::vector<int> coreSize;
	ClusterVector burgersVector;
	DislocationNode* nodes[2];
	DislocationSegment* replacedWith;

	DislocationSegment(const ClusterVector& b, DislocationNode* forward, DislocationNode* backward)
		: id(-1), burgersVector(b), replacedWith(nullptr){

		nodes[0] = forward;
		nodes[1] = backward;

		forward->segment = this;
		forward->oppositeNode = backward;
		forward->junctionRing = forward;

		backward->segment = this;
		backward->oppositeNode = forward;
		backward->junctionRing = backward;
	}

	[[nodiscard]] DislocationNode &forwardNode() const{
		return *nodes[0];
	}

	[[nodiscard]] DislocationNode& backwardNode() const{
		return *nodes[1];
	}

	[[nodiscard]] bool isClosedLoop() const{
		return (nodes[0]->junctionRing == nodes[1]) && (nodes[1]->junctionRing == nodes[0]);
	}

	[[nodiscard]] bool isInfiniteLine() const{
		return isClosedLoop() && !line.front().equals(line.back(), CA_ATOM_VECTOR_EPSILON);
	}

	[[nodiscard]] bool isDegenerate() const{
		return line.size() <= 1;
	}

	[[nodiscard]] double calculateLength() const noexcept {
		double length = 0;
		if(line.size() < 2) return 0;
		auto next = line.begin();
		for(auto i1 = line.begin(), i2 = ++next; i2 != line.end(); i1 = i2, ++i2){
			length += (*i1 - *i2).length();
		}
		return length;
	}

	void flipOrientation(){
		burgersVector = -burgersVector;
		std::swap(nodes[0], nodes[1]);
		std::ranges::reverse(line);
		std::ranges::reverse(coreSize);
	}

	// line and coreSize hold one entry per point: these change both together.
	void addPoint(bool atForwardEnd, const Point3& point, int core){
		if(atForwardEnd){
			line.push_back(point);
			coreSize.push_back(core);
		}else{
			line.insert(line.begin(), point);
			coreSize.insert(coreSize.begin(), core);
		}
	}

	void trimEnds(std::size_t front, std::size_t back){
		line.erase(line.begin(), line.begin() + front);
		line.erase(line.end() - back, line.end());
		coreSize.erase(coreSize.begin(), coreSize.begin() + front);
		coreSize.erase(coreSize.end() - back, coreSize.end());
	}
};

inline bool DislocationNode::isForwardNode() const{
	return &segment->forwardNode() == this;
}

inline bool DislocationNode::isBackwardNode() const{
	return &segment->backwardNode() == this;
}

inline ClusterVector DislocationNode::burgersVector() const{
	return isForwardNode() ? segment->burgersVector : -segment->burgersVector;
}

inline const Point3& DislocationNode::position() const{
	return isForwardNode() ? segment->line.back() : segment->line.front();
}

class DislocationNetwork{
public:
	[[nodiscard]] const std::vector<DislocationSegment*>& segments() const{
		return _segments;
	}

	void smoothDislocationLines(double lineSmoothingLevel, double linePointInterval);

	// The number of arms of each junction, a ring of three or more segment ends, counted once per ring. A ring of two
	// ends is the seam of a closed loop, not a junction.
	[[nodiscard]] std::vector<int> junctionArmCounts() const;

	DislocationSegment* createSegment(const ClusterVector& burgersVector);
	void discardSegment(DislocationSegment* segment);

private:
	// Replaces runs of consecutive points with their mean. A run grows while its point count squared stays below
	// linePointInterval times its summed core size, so lines with wide cores keep fewer points. An open line keeps its ends.
	static std::vector<Point3> coarsenDislocationLine(
		double linePointInterval,
		const std::vector<Point3>& input,
		const std::vector<int>& coreSize,
		bool isClosedLoop,
		bool isInfiniteLine
	);
	static void smoothDislocationLine(double smoothingLevel, std::vector<Point3>& line, bool isLoop);

	MemoryPool<DislocationNode> _nodePool;
	std::vector<DislocationSegment*> _segments;
	MemoryPool<DislocationSegment> _segmentPool;
};

}
