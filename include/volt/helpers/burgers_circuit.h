#pragma once

#include <cassert>
#include <cstddef>
#include <vector>

#include <volt/core/volt.h>
#include <volt/helpers/cluster_vector.h>
#include <volt/pipeline/interface_mesh.h>

namespace Volt{

struct DislocationNode;

struct BurgersCircuit{
	InterfaceMesh::Edge* firstEdge = nullptr;

    InterfaceMesh::Edge* lastEdge = nullptr;

    std::vector<InterfaceMesh::Edge*> segmentMeshCap;

    std::size_t numPreliminaryPoints = 0;

    DislocationNode* dislocationNode = nullptr;

    std::size_t edgeCount = 0;

    bool isCompletelyBlocked = false;

    bool isDangling = true;

	BurgersCircuit() noexcept = default;

	[[nodiscard]]
	ClusterVector calculateBurgersVector() const noexcept{
		Vector3 b{};
		Matrix3 tm = Matrix3::Identity();
		auto* edge = firstEdge;
        if(!edge) return ClusterVector();

        size_t safety_counter = 0;
		do{
			b += tm * edge->clusterVector;

			if(edge->clusterTransition && !edge->clusterTransition->isSelfTransition()){
				tm = tm * edge->clusterTransition->reverse->tm;
			}

			edge = edge->nextCircuitEdge;
            if(++safety_counter > edgeCount + 100) break; 
		}while(edge != nullptr && edge != firstEdge);

        if(!firstEdge || !firstEdge->clusterTransition) return ClusterVector(b, nullptr);
		return ClusterVector(b, firstEdge->clusterTransition->cluster1);
	}

	[[nodiscard]]
	Point3 calculateCenter() const noexcept{
		Vector3 center{};
		Vector3 current{};
		auto* edge = firstEdge;
		if(!edge || edgeCount == 0) return Point3::Origin();

        size_t safety_counter = 0;
		do{
			center += current;
			current += edge->physicalVector;
			edge = edge->nextCircuitEdge;
            if(++safety_counter > edgeCount + 100) break;
		}while(edge != nullptr && edge != firstEdge);

		return firstEdge->vertex1()->pos() + (center / static_cast<double>(edgeCount));
	}

	[[nodiscard]]
	size_t countEdges() const noexcept{
		size_t cnt = 0;
		auto *edge = firstEdge;
        if(!edge) return 0;
		do{
			++cnt;
			edge = edge->nextCircuitEdge;
            if(cnt > 10000) break;
		}while(edge != nullptr && edge != firstEdge);
		return cnt;
	}

	[[nodiscard]]
	InterfaceMesh::Edge* getEdge(size_t idx) const noexcept{
        if(edgeCount == 0) return nullptr;
		auto *edge = firstEdge;
		while(idx-- && edge != nullptr){
			edge = edge->nextCircuitEdge;
		}
		return edge;
	}

	void storeCircuit() noexcept{
		assert(segmentMeshCap.empty());
		segmentMeshCap.reserve(edgeCount);
		auto* edge = firstEdge;
		do{
			segmentMeshCap.push_back(edge);
			edge = edge->nextCircuitEdge;
		}while(edge != firstEdge);
		assert(segmentMeshCap.size() >= 2);
	}
};

}
