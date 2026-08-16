#include <volt/pipeline/burgers_loop_builder.h>
#include <volt/pipeline/interface_mesh.h>
#include <volt/helpers/triangle_intersection.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/blocked_range.h>
#include <tbb/spin_mutex.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <ranges>
#include <array>
#include <atomic>
#include <algorithm>

namespace Volt{

BurgersCircuit* BurgersLoopBuilder::allocateCircuit(){
    BurgersCircuit* circuit = nullptr;
    tbb::spin_mutex::scoped_lock lock(_circuit_pool_mutex);
    
    if(_unusedCircuit != nullptr){
        circuit = _unusedCircuit;
        _unusedCircuit = nullptr;
    }
    
    if(circuit == nullptr){
        return _circuitPool.construct();
    }
    return circuit;
}

bool BurgersLoopBuilder::traceDislocationSegments(){
    mesh().clearFaceFlag(0);

    if(_markCoreAtoms){
        const double alpha = DefectCellIndex::kAlphaScale
            * _mesh.structureAnalysis().maximumNeighborDistance();
        _defectCells.emplace(_mesh.tessellation(), alpha);
        _coreCellClaims = std::vector<CoreCellClaim>(_defectCells->cellCount());
    }

    for(int L : std::views::iota(3, _maxExtendedBurgersCircuitSize + 1)){
        auto dangling = _danglingNodes;

        if(!dangling.empty()){
            for(auto* node : dangling){
                traceSegment(*node->segment, *node, L, L <= _maxBurgersCircuitSize);
            }
        }

        if((L & 1) && L <= _maxBurgersCircuitSize){
            if(!findPrimarySegments(L)){
                return false;
            }
        }
        (void) joinSegments(L);

        if(L >= _maxBurgersCircuitSize && !dangling.empty()){
            tbb::parallel_for_each(dangling.begin(), dangling.end(), [&](DislocationNode* node){
                auto* C = node->circuit;
                if(C->isDangling && C->segmentMeshCap.empty()){
                    C->storeCircuit();
                    C->numPreliminaryPoints = 0;
                }
            });
        }
    }

    return true;
}

void BurgersLoopBuilder::discardCircuit(BurgersCircuit* circuit){
    tbb::spin_mutex::scoped_lock lock(_circuit_pool_mutex);
    assert(_unusedCircuit == nullptr);
    _unusedCircuit = circuit;
}

void BurgersLoopBuilder::finishDislocationSegments(std::string_view referenceTopologyName){
    auto& segs = network().segments();

    tbb::parallel_for(tbb::blocked_range<size_t>(0, segs.size()), 
        [&](const tbb::blocked_range<size_t>& r){
            for(size_t i = r.begin(); i != r.end(); ++i){
                auto* s = segs[i];
                auto pre  = s->backwardNode().circuit->numPreliminaryPoints;
                auto post = s->forwardNode().circuit->numPreliminaryPoints;
                s->id = static_cast<int>(i);

                auto& line = s->line;
                auto& core = s->coreSize;

                line.erase(line.begin(), line.begin() + pre);
                line.erase(line.end() - post, line.end());
                core.erase(core.begin(), core.begin() + pre);
                core.erase(core.end() - post, core.end());
            }
	});

    tbb::parallel_for(tbb::blocked_range<size_t>(0, segs.size()), 
        [&](const tbb::blocked_range<size_t>& r){
            for(size_t i = r.begin(); i != r.end(); ++i){
                auto* s = segs[i];
                auto* orig = s->burgersVector.cluster();
                if(orig->topologyName != referenceTopologyName){
                    for(auto* t = orig->transitions; t && t->distance <= 1; t = t->next){
                        if(t->cluster2->topologyName == referenceTopologyName){
                            s->burgersVector = ClusterVector(
                                t->transform(s->burgersVector.localVec()),
                                t->cluster2
                            );
                            break;
                        }
                    }
                }
            }
	});

    tbb::parallel_for(tbb::blocked_range<size_t>(0, segs.size()), 
        [&](const tbb::blocked_range<size_t>& r){
            for(size_t i = r.begin(); i != r.end(); ++i){
                auto* s = segs[i];
                auto& line = s->line;
                if (line.empty()) continue;
                Vector3 dir = line.back() - line.front();
                if(dir.isZero(CA_ATOM_VECTOR_EPSILON)) continue;

                auto absx = std::abs(dir.x()), absy = std::abs(dir.y()), absz = std::abs(dir.z());
                if((absx >= absy && absx >= absz && dir.x() < 0) ||
                   (absy >= absx && absy >= absz && dir.y() < 0) ||
                   (absz >= absx && absz >= absy && dir.z() < 0))
                {
                    s->flipOrientation();
                }
            }
	});
}

struct BurgersCircuitSearchStruct{
    InterfaceMesh::Vertex* node;
    Point3 latticeCoord;
    Matrix3 tm;
    int recursiveDepth;
    InterfaceMesh::Edge* predecessorEdge;
    BurgersCircuitSearchStruct* nextToProcess;
};

using SearchNode = BurgersLoopBuilder::SearchNode;

bool BurgersLoopBuilder::findPrimarySegments(int maxBurgersCircuitSize){
    const int searchDepth = (maxBurgersCircuitSize - 1) / 2;
    assert(searchDepth >= 1);

    auto& verts = mesh().vertices();

    struct Candidate{
        InterfaceMesh::Edge* edge;
    };
    std::vector<Candidate> candidates;
    tbb::spin_mutex candidatesMutex;

    tbb::parallel_for(tbb::blocked_range<size_t>(0, verts.size(), 16),
        [&](const tbb::blocked_range<size_t>& range){

        MemoryPool<SearchNode> localPool;
        std::vector<SearchNode*> queue;
        queue.reserve(512);
        std::unordered_map<InterfaceMesh::Vertex*, SearchNode*> visited;
        visited.reserve(256);

        for(size_t vi = range.begin(); vi < range.end(); ++vi){
            auto* startVert = verts[vi];

            if(startVert->edges() == nullptr) continue;

            bool allEdgesOwned = true;
            for(auto* e = startVert->edges(); e != nullptr; e = e->nextVertexEdge()){
                if(!e->nextCircuitEdge && !e->face()->circuit &&
                   !e->_claimedForCircuit.load(std::memory_order_relaxed)){
                    allEdgesOwned = false;
                    break;
                }
            }
            if(allEdgesOwned) continue;

            queue.clear();
            visited.clear();

            auto* root = localPool.construct();
            root->node = startVert;
            root->coord = Point3::Origin();
            root->tm.setIdentity();
            root->depth = 0;
            root->viaEdge = nullptr;
            queue.push_back(root);
            visited[startVert] = root;

            bool found = false;

            for(size_t qi = 0; qi < queue.size() && !found; ++qi){
                auto* cur = queue[qi];

                for(auto* edge = cur->node->edges(); edge != nullptr && !found; edge = edge->nextVertexEdge()){
                    if(edge->nextCircuitEdge || edge->face()->circuit) continue;
                    if(edge->_claimedForCircuit.load(std::memory_order_relaxed)) continue;

                    auto* nbVert = edge->vertex2();
                    Point3 nbCoord = cur->coord + cur->tm * edge->clusterVector;

                    auto it = visited.find(nbVert);
                    if(it != visited.end()){
                        auto* prevNode = it->second;
                        Vector3 b = prevNode->coord - nbCoord;
                        if(!b.isZero(CA_LATTICE_VECTOR_EPSILON)){
                            Matrix3 R = cur->tm * edge->clusterTransition->reverse->tm;
                            if(R.equals(prevNode->tm, CA_TRANSITION_MATRIX_EPSILON)){
                                tbb::spin_mutex::scoped_lock lock(candidatesMutex);
                                candidates.push_back({edge});
                                found = true;
                            }
                        }
                    }else if(cur->depth < searchDepth){
                        auto* nb = localPool.construct();
                        nb->node = nbVert;
                        nb->coord = nbCoord;
                        nb->depth = cur->depth + 1;
                        nb->viaEdge = edge;
                        nb->tm = edge->clusterTransition->isSelfTransition()
                            ? cur->tm
                            : cur->tm * edge->clusterTransition->reverse->tm;
                        queue.push_back(nb);
                        visited[nbVert] = nb;
                    }
                }
            }

            localPool.clear(true);
        }
    });

    std::sort(candidates.begin(), candidates.end(),
        [](const Candidate& a, const Candidate& b){
            const int a1 = a.edge->vertex1()->index();
            const int b1 = b.edge->vertex1()->index();
            if(a1 != b1) return a1 < b1;
            return a.edge->vertex2()->index() < b.edge->vertex2()->index();
        });

    MemoryPool<SearchNode> pool;
    std::vector<SearchNode*> queue;
    queue.reserve(1024);

    for(auto& cand : candidates){
        auto* startVert = cand.edge->vertex1();

        if(cand.edge->circuit != nullptr || cand.edge->_claimedForCircuit.load(std::memory_order_relaxed))
            continue;

        queue.clear();

        auto* root = pool.construct();
        root->node = startVert;
        root->coord = Point3::Origin();
        root->tm.setIdentity();
        root->depth = 0;
        root->viaEdge = nullptr;
        startVert->burgersSearchStruct.store(reinterpret_cast<BurgersCircuitSearchStruct*>(root), std::memory_order_relaxed);
        queue.push_back(root);

        bool found = false;
        for(size_t qi = 0; qi < queue.size() && !found; ++qi){
            auto* cur = queue[qi];

            for(auto* edge = cur->node->edges(); edge != nullptr && !found; edge = edge->nextVertexEdge()){
                if(edge->nextCircuitEdge || edge->face()->circuit) continue;

                auto* nbVert = edge->vertex2();
                Point3 nbCoord = cur->coord + cur->tm * edge->clusterVector;

                auto* prevStruct = nbVert->burgersSearchStruct.load(std::memory_order_relaxed);
                if(prevStruct){
                    auto* prev = reinterpret_cast<SearchNode*>(prevStruct);
                    Vector3 b = prev->coord - nbCoord;
                    if(!b.isZero(CA_LATTICE_VECTOR_EPSILON)){
                        Matrix3 R = cur->tm * edge->clusterTransition->reverse->tm;
                        if(R.equals(prev->tm, CA_TRANSITION_MATRIX_EPSILON)){
                            found = createBurgersCircuitOriginal(edge, maxBurgersCircuitSize);
                        }
                    }
                }else if(cur->depth < searchDepth){
                    auto* nb = pool.construct();
                    nb->node = nbVert;
                    nb->coord = nbCoord;
                    nb->depth = cur->depth + 1;
                    nb->viaEdge = edge;
                    nb->tm = edge->clusterTransition->isSelfTransition()
                        ? cur->tm
                        : cur->tm * edge->clusterTransition->reverse->tm;
                    nbVert->burgersSearchStruct.store(reinterpret_cast<BurgersCircuitSearchStruct*>(nb), std::memory_order_relaxed);
                    queue.push_back(nb);
                }
            }
        }

        for(auto* sn : queue){
            sn->node->burgersSearchStruct.store(nullptr, std::memory_order_relaxed);
        }
        pool.clear(true);
    }

    return true;
}

bool BurgersLoopBuilder::createBurgersCircuitOriginal(InterfaceMesh::Edge* edge, int maxBurgersCircuitSize){
	if(edge->circuit != nullptr) return false;

	InterfaceMesh::Vertex* currentNode = edge->vertex1();
	InterfaceMesh::Vertex* neighborNode = edge->vertex2();
	auto* currentStruct = currentNode->burgersSearchStruct.load(std::memory_order_relaxed);
	auto* neighborStruct = neighborNode->burgersSearchStruct.load(std::memory_order_relaxed);
	assert(currentStruct != neighborStruct);

	BurgersCircuit* forwardCircuit = allocateCircuit();
	forwardCircuit->edgeCount = 1;
	forwardCircuit->firstEdge = forwardCircuit->lastEdge = edge->oppositeEdge();
	assert(forwardCircuit->firstEdge->circuit == nullptr);
	forwardCircuit->firstEdge->circuit = forwardCircuit;

	std::unordered_set<InterfaceMesh::Vertex*> firstBranchSet;
	for(BurgersCircuitSearchStruct* a = currentStruct; ; a = a->predecessorEdge->vertex1()->burgersSearchStruct.load(std::memory_order_relaxed)){
		firstBranchSet.insert(a->node);
		if(a->predecessorEdge == nullptr) break;
	}

	for(BurgersCircuitSearchStruct* a = neighborStruct; ; a = a->predecessorEdge->vertex1()->burgersSearchStruct.load(std::memory_order_relaxed)){
		if(firstBranchSet.count(a->node)) break;
		assert(a->predecessorEdge != nullptr);
		assert(a->predecessorEdge->circuit == nullptr);
		a->predecessorEdge->nextCircuitEdge = forwardCircuit->firstEdge;
		forwardCircuit->firstEdge = a->predecessorEdge;
		forwardCircuit->edgeCount++;
		forwardCircuit->firstEdge->circuit = forwardCircuit;
	}

	InterfaceMesh::Vertex* meetNode = nullptr;
	for(BurgersCircuitSearchStruct* a = neighborStruct; ; a = a->predecessorEdge->vertex1()->burgersSearchStruct.load(std::memory_order_relaxed)){
		if(firstBranchSet.count(a->node)){
			meetNode = a->node;
			break;
		}
		if(a->predecessorEdge == nullptr) break;
	}

	for(BurgersCircuitSearchStruct* a = currentStruct; a->node != meetNode; a = a->predecessorEdge->vertex1()->burgersSearchStruct.load(std::memory_order_relaxed)){
		assert(a->predecessorEdge != nullptr);
		assert(a->predecessorEdge->oppositeEdge()->circuit == nullptr);
		forwardCircuit->lastEdge->nextCircuitEdge = a->predecessorEdge->oppositeEdge();
		forwardCircuit->lastEdge = forwardCircuit->lastEdge->nextCircuitEdge;
		forwardCircuit->edgeCount++;
		forwardCircuit->lastEdge->circuit = forwardCircuit;
	}

	forwardCircuit->lastEdge->nextCircuitEdge = forwardCircuit->firstEdge;
	assert(forwardCircuit->firstEdge != forwardCircuit->firstEdge->nextCircuitEdge);
	assert(forwardCircuit->countEdges() == forwardCircuit->edgeCount);
	assert(forwardCircuit->edgeCount >= 3);

	InterfaceMesh::Edge* e = forwardCircuit->firstEdge;
	Vector3 edgeSum = Vector3::Zero();
	Matrix3 frankRotation = Matrix3::Identity();
	Vector3 b = Vector3::Zero();
	do{
		edgeSum += e->physicalVector;
		b += frankRotation * e->clusterVector;
		if(!e->clusterTransition->isSelfTransition())
			frankRotation = frankRotation * e->clusterTransition->reverse->tm;
		e = e->nextCircuitEdge;
	}while(e != forwardCircuit->firstEdge);
	assert(frankRotation.equals(Matrix3::Identity(), CA_TRANSITION_MATRIX_EPSILON));

	bool intersects = intersectsOtherCircuits(forwardCircuit);
	const bool zeroBurgers = b.isZero(CA_LATTICE_VECTOR_EPSILON);
	const bool badEdgeSum = !edgeSum.isZero(CA_ATOM_VECTOR_EPSILON);
	if(zeroBurgers || badEdgeSum || intersects){
		e = forwardCircuit->firstEdge;
		do{
			InterfaceMesh::Edge* nextEdge = e->nextCircuitEdge;
			e->nextCircuitEdge = nullptr;
			e->circuit = nullptr;
			e = nextEdge;
		}while(e != forwardCircuit->firstEdge);
		discardCircuit(forwardCircuit);
		return intersects;
	}

	createAndTraceSegment(ClusterVector(b, forwardCircuit->firstEdge->clusterTransition->cluster1), forwardCircuit, maxBurgersCircuitSize);
	return true;
}

void BurgersLoopBuilder::createAndTraceSegment(const ClusterVector& burgersVector, BurgersCircuit* forwardCircuit, int maxCircuitLength){
	BurgersCircuit* backwardCircuit = buildReverseCircuit(forwardCircuit);

	DislocationSegment* segment = network().createSegment(burgersVector);
	segment->forwardNode().circuit = forwardCircuit;
	segment->backwardNode().circuit = backwardCircuit;
	forwardCircuit->dislocationNode = &segment->forwardNode();
	backwardCircuit->dislocationNode = &segment->backwardNode();

	_danglingNodes.push_back(&segment->forwardNode());
	_danglingNodes.push_back(&segment->backwardNode());

	segment->line.push_back(backwardCircuit->calculateCenter());
	segment->coreSize.push_back(backwardCircuit->countEdges());

	appendLinePoint(segment->forwardNode());

	traceSegment(*segment, segment->forwardNode(), maxCircuitLength, true);
	traceSegment(*segment, segment->backwardNode(), maxCircuitLength, true);
}

bool BurgersLoopBuilder::intersectsOtherCircuits(BurgersCircuit* circuit){
    InterfaceMesh::Edge* startEdge1 = circuit->firstEdge;
    for(InterfaceMesh::Edge* edge1 = startEdge1; ; edge1 = edge1->nextCircuitEdge){
        InterfaceMesh::Edge* edge2 = edge1->nextCircuitEdge;

        if(edge1 != edge2->oppositeEdge()){
            InterfaceMesh::Edge* sentinel = edge1->oppositeEdge();
            InterfaceMesh::Edge* cur = sentinel;

            do{
                InterfaceMesh::Edge* prev = cur->prevFaceEdge();
                if(prev->circuit){
                    int goingOutside = 0;
					int goingInside = 0;
                    
                    if(prev->nextCircuitEdge && 
                       edge2->oppositeEdge() && 
                       edge1->oppositeEdge() &&
                       prev->nextCircuitEdge->vertex1() == prev->vertex2()) {
                        circuitCircuitIntersection(
                            edge2->oppositeEdge(),
                            edge1->oppositeEdge(),
                            prev,
                            prev->nextCircuitEdge,
                            goingOutside,
                            goingInside
                        );

                        if(goingOutside){
                            return true;
                        }
                    }
                }
                cur = prev->oppositeEdge();
            }while (cur != sentinel);
        }

        if(edge2 == startEdge1){
            break;
        }
    }

    return false;
}

BurgersCircuit* BurgersLoopBuilder::buildReverseCircuit(BurgersCircuit* forwardCircuit){
	BurgersCircuit* backwardCircuit = allocateCircuit();

	backwardCircuit->edgeCount = 0;
	backwardCircuit->firstEdge = nullptr;
	backwardCircuit->lastEdge = nullptr;
	InterfaceMesh::Edge* edge1 = forwardCircuit->firstEdge;
	do{
		InterfaceMesh::Edge* edge2 = edge1->nextCircuitEdge;
		InterfaceMesh::Edge* oppositeEdge1 = edge1->oppositeEdge();
		InterfaceMesh::Edge* oppositeEdge2 = edge2->oppositeEdge();
		InterfaceMesh::Face* facet1 = oppositeEdge1->face();
		InterfaceMesh::Face* facet2 = oppositeEdge2->face();
		assert(facet1 != nullptr && facet2 != nullptr);
		assert(facet1->circuit == nullptr || facet1->circuit == backwardCircuit);
		assert(facet2->circuit == nullptr || facet2->circuit == backwardCircuit);
		assert(edge1->vertex2() == edge2->vertex1());
		assert((edge1->clusterVector + oppositeEdge1->clusterTransition->tm * oppositeEdge1->clusterVector).isZero(CA_LATTICE_VECTOR_EPSILON));
		assert((edge2->clusterVector + oppositeEdge2->clusterTransition->tm * oppositeEdge2->clusterVector).isZero(CA_LATTICE_VECTOR_EPSILON));

		if(facet1 != facet2){
			InterfaceMesh::Edge* outerEdge1 = oppositeEdge1->nextFaceEdge()->oppositeEdge();
			InterfaceMesh::Edge* innerEdge1 = oppositeEdge1->prevFaceEdge()->oppositeEdge();
			InterfaceMesh::Edge* outerEdge2 = oppositeEdge2->prevFaceEdge()->oppositeEdge();
			InterfaceMesh::Edge* innerEdge2 = oppositeEdge2->nextFaceEdge()->oppositeEdge();
			assert(innerEdge1 != nullptr && innerEdge2 != nullptr);
			assert(innerEdge1->vertex1() == edge1->vertex2());
			assert(innerEdge2->vertex2() == edge1->vertex2());
			assert(innerEdge1->vertex1() == innerEdge2->vertex2());
			assert(innerEdge1->circuit == nullptr || innerEdge1->circuit == backwardCircuit);
			assert(innerEdge2->circuit == nullptr || innerEdge2->circuit == backwardCircuit);
			facet1->setFlag(1);
			facet1->circuit = backwardCircuit;
			facet2->setFlag(1);
			facet2->circuit = backwardCircuit;
			innerEdge1->circuit = backwardCircuit;
			innerEdge2->circuit = backwardCircuit;
			innerEdge2->nextCircuitEdge = innerEdge1;

			if(backwardCircuit->lastEdge == nullptr){
				assert(backwardCircuit->firstEdge == nullptr);
				assert(innerEdge1->nextCircuitEdge == nullptr);
				backwardCircuit->lastEdge = innerEdge1;
				backwardCircuit->firstEdge = innerEdge2;
				backwardCircuit->edgeCount += 2;
			}else if(backwardCircuit->lastEdge != innerEdge2){
				if(innerEdge1 != backwardCircuit->firstEdge){
					innerEdge1->nextCircuitEdge = backwardCircuit->firstEdge;
					backwardCircuit->edgeCount += 2;
				}else{
					backwardCircuit->edgeCount += 1;
				}
				backwardCircuit->firstEdge = innerEdge2;
			}else if(backwardCircuit->firstEdge != innerEdge1){
				innerEdge1->nextCircuitEdge = backwardCircuit->firstEdge;
				backwardCircuit->firstEdge = innerEdge1;
				backwardCircuit->edgeCount += 1;
			}

			assert(innerEdge1->vertex1() != innerEdge1->vertex2());
			assert(innerEdge2->vertex1() != innerEdge2->vertex2());
		}

		edge1 = edge2;
	}while(edge1 != forwardCircuit->firstEdge);

	assert(backwardCircuit->lastEdge->vertex2() == backwardCircuit->firstEdge->vertex1());
	assert(backwardCircuit->lastEdge->nextCircuitEdge == nullptr || backwardCircuit->lastEdge->nextCircuitEdge == backwardCircuit->firstEdge);

	backwardCircuit->lastEdge->nextCircuitEdge = backwardCircuit->firstEdge;

	assert(backwardCircuit->firstEdge != backwardCircuit->firstEdge->nextCircuitEdge);
	assert(backwardCircuit->countEdges() == backwardCircuit->edgeCount);
	assert(backwardCircuit->edgeCount >= 3);
	assert(!backwardCircuit->calculateBurgersVector().localVec().isZero(CA_LATTICE_VECTOR_EPSILON));

	return backwardCircuit;
}

void BurgersLoopBuilder::traceSegment(DislocationSegment& segment, DislocationNode& node, int maxCircuitLength, bool isPrimarySegment){
    BurgersCircuit& circuit = *node.circuit;
    assert(circuit.countEdges() == circuit.edgeCount);
    assert(circuit.isDangling);

    for(;;){
        size_t idx = _edgeStartIndex.fetch_add(1, std::memory_order_relaxed);
        int edgeIndex = (idx % circuit.edgeCount);

        InterfaceMesh::Edge* firstEdge = circuit.getEdge(edgeIndex);

        InterfaceMesh::Edge* edge0 = firstEdge;
        InterfaceMesh::Edge* edge1 = edge0->nextCircuitEdge;
        InterfaceMesh::Edge* edge2 = edge1->nextCircuitEdge;
        assert(edge1->circuit == &circuit);
        int counter = 0;
        do{
            assert(circuit.edgeCount >= 3);
            
            ClusterVector burgersVec = circuit.calculateBurgersVector();
            if(burgersVec.localVec().isZero(CA_LATTICE_VECTOR_EPSILON)) {
                std::cerr << "Warning: Burgers vector is zero for circuit with " << circuit.edgeCount 
                         << " edges. Attempting to continue..." << std::endl;
                if(circuit.edgeCount <= 3) {
                    std::cerr << "Error: Cannot recover circuit with only " << circuit.edgeCount << " edges" << std::endl;
                    return;
                }
            }
            
            assert(circuit.countEdges() == circuit.edgeCount);
            assert(edge0->circuit == &circuit && edge1->circuit == &circuit && edge2->circuit == &circuit);

            bool wasShortened = false;
            if(tryRemoveTwoCircuitEdges(edge0, edge1, edge2)){
                wasShortened = true;
			}else if(tryRemoveThreeCircuitEdges(edge0, edge1, edge2, isPrimarySegment)){
                wasShortened = true;
			}else if(tryRemoveOneCircuitEdge(edge0, edge1, edge2, isPrimarySegment)){
                wasShortened = true;
			}else if(trySweepTwoFacets(edge0, edge1, edge2, isPrimarySegment)){
                wasShortened = true;
			}

            if(wasShortened){
                appendLinePoint(node);
                counter = -1;
            }

            edge0 = edge1;
            edge1 = edge2;
            edge2 = edge2->nextCircuitEdge;
            counter++;
        }while(counter <= circuit.edgeCount);
        assert(circuit.edgeCount >= 3);
        assert(circuit.countEdges() == circuit.edgeCount);

        if(circuit.edgeCount >= maxCircuitLength) break;

        bool wasExtended = false;

        idx = _edgeStartIndex.fetch_add(1, std::memory_order_relaxed);
        edgeIndex = (idx % circuit.edgeCount);

        firstEdge = circuit.getEdge(edgeIndex);

        edge0 = firstEdge;
        edge1 = firstEdge->nextCircuitEdge;
        do{
            if(tryInsertOneCircuitEdge(edge0, edge1, isPrimarySegment)){
                wasExtended = true;
                appendLinePoint(node);
                break;
            }

            edge0 = edge1;
            edge1 = edge1->nextCircuitEdge;
        }while(edge0 != firstEdge);
        if(!wasExtended) break;
    }
}

bool BurgersLoopBuilder::tryRemoveTwoCircuitEdges(InterfaceMesh::Edge*& edge0, InterfaceMesh::Edge*& edge1, InterfaceMesh::Edge*& edge2){
	if(edge1 != edge2->oppositeEdge()) return false;

	BurgersCircuit* circuit = edge0->circuit;
	assert(circuit->edgeCount >= 4);
	edge0->nextCircuitEdge = edge2->nextCircuitEdge;
	
	if(edge0 == circuit->lastEdge){
		circuit->firstEdge = circuit->lastEdge->nextCircuitEdge;
	}else if(edge1 == circuit->lastEdge){
		circuit->lastEdge = edge0;
		circuit->firstEdge = edge0->nextCircuitEdge;
	}else if(edge2 == circuit->lastEdge){
		circuit->lastEdge = edge0;
	}

	circuit->edgeCount -= 2;
	assert(circuit->edgeCount >= 0);

	edge1 = edge0->nextCircuitEdge;
	edge2 = edge1->nextCircuitEdge;
	return true;
}

bool BurgersLoopBuilder::tryRemoveThreeCircuitEdges(
	InterfaceMesh::Edge*& edge0, 
	InterfaceMesh::Edge*& edge1, 
	InterfaceMesh::Edge*& edge2,
	bool isPrimarySegment
){
	InterfaceMesh::Face* facet1 = edge1->face();
	InterfaceMesh::Face* facet2 = edge2->face();

	if(facet2 != facet1 || facet1->circuit != nullptr) return false;

	BurgersCircuit* circuit = edge0->circuit;
	assert(circuit->edgeCount > 2);
	InterfaceMesh::Edge* edge3 = edge2->nextCircuitEdge;

	if(edge3->face() != facet1) return false;
	assert(circuit->edgeCount > 4);

	edge0->nextCircuitEdge = edge3->nextCircuitEdge;

	if(edge2 == circuit->firstEdge || edge3 == circuit->firstEdge){
		circuit->firstEdge = edge3->nextCircuitEdge;
		circuit->lastEdge = edge0;
	}else if(edge1 == circuit->firstEdge){
		circuit->firstEdge = edge3->nextCircuitEdge;
		assert(circuit->lastEdge == edge0);
	}else if(edge3 == circuit->lastEdge){
		circuit->lastEdge = edge0;
	}

	circuit->edgeCount -= 3;
	edge1 = edge3->nextCircuitEdge;
	edge2 = edge1->nextCircuitEdge;

	facet1->circuit = circuit;
	if(isPrimarySegment) facet1->setFlag(1);

	return true;
}

bool BurgersLoopBuilder::tryRemoveOneCircuitEdge(
	InterfaceMesh::Edge*& edge0, 
	InterfaceMesh::Edge*& edge1, 
	InterfaceMesh::Edge*& edge2, 
	bool isPrimarySegment
){
	InterfaceMesh::Face* facet1 = edge1->face();
	InterfaceMesh::Face* facet2 = edge2->face();
	if(facet2 != facet1 || facet1->circuit != nullptr) return false;

	BurgersCircuit* circuit = edge0->circuit;
	assert(circuit->edgeCount > 2);

	if(edge0->face() == facet1) return false;

	InterfaceMesh::Edge* shortEdge = edge1->prevFaceEdge()->oppositeEdge();
	assert(shortEdge->vertex1() == edge1->vertex1());
	assert(shortEdge->vertex2() == edge2->vertex2());

	if(shortEdge->circuit != nullptr) return false;
	if(!shortEdge->tryClaimForCircuit()) return false;

	assert(shortEdge->nextCircuitEdge == nullptr);
	shortEdge->nextCircuitEdge = edge2->nextCircuitEdge;
	assert(shortEdge != edge2->nextCircuitEdge->oppositeEdge());
	assert(shortEdge != edge0->oppositeEdge());
	edge0->nextCircuitEdge = shortEdge;
	if(edge0 == circuit->lastEdge){
		assert(circuit->lastEdge != edge2);
		assert(circuit->firstEdge == edge1);
		assert(shortEdge != circuit->lastEdge->oppositeEdge());
		circuit->firstEdge = shortEdge;
	}

	if(edge2 == circuit->lastEdge){
		circuit->lastEdge = shortEdge;
	}else if(edge2 == circuit->firstEdge){
		circuit->firstEdge = shortEdge->nextCircuitEdge;
		circuit->lastEdge = shortEdge;
	}

	circuit->edgeCount -= 1;
	edge1 = shortEdge;
	edge2 = shortEdge->nextCircuitEdge;
	shortEdge->circuit = circuit;

	facet1->circuit = circuit;
	if(isPrimarySegment) facet1->setFlag(1);

	return true;
}

bool BurgersLoopBuilder::trySweepTwoFacets(
	InterfaceMesh::Edge*& edge0, 
	InterfaceMesh::Edge*& edge1, 
	InterfaceMesh::Edge*& edge2, 
	bool isPrimarySegment
){
	InterfaceMesh::Face* facet1 = edge1->face();
	InterfaceMesh::Face* facet2 = edge2->face();

	if(facet1->circuit != nullptr || facet2->circuit != nullptr) return false;

	BurgersCircuit* circuit = edge0->circuit;
	if(facet1 == facet2 || circuit->edgeCount <= 2) return false;

	InterfaceMesh::Edge* outerEdge1 = edge1->prevFaceEdge()->oppositeEdge();
	InterfaceMesh::Edge* innerEdge1 = edge1->nextFaceEdge();
	InterfaceMesh::Edge* outerEdge2 = edge2->nextFaceEdge()->oppositeEdge();
	InterfaceMesh::Edge* innerEdge2 = edge2->prevFaceEdge();

	if(innerEdge1 != innerEdge2->oppositeEdge() || outerEdge1->circuit != nullptr || outerEdge2->circuit != nullptr)
		return false;

	if(!outerEdge1->tryClaimForCircuit()){
		return false;
	}
	if(!outerEdge2->tryClaimForCircuit()){
		outerEdge1->releaseCircuitClaim();
		return false;
	}

	assert(outerEdge1->nextCircuitEdge == nullptr);
	assert(outerEdge2->nextCircuitEdge == nullptr);
	outerEdge1->nextCircuitEdge = outerEdge2;
	outerEdge2->nextCircuitEdge = edge2->nextCircuitEdge;
	edge0->nextCircuitEdge = outerEdge1;

	if(edge0 == circuit->lastEdge){
		circuit->firstEdge = outerEdge1;
	}else if(edge1 == circuit->lastEdge){
		circuit->lastEdge = outerEdge1;
		circuit->firstEdge = outerEdge2;
	}else if(edge2 == circuit->lastEdge){
		circuit->lastEdge = outerEdge2;
	}

	outerEdge1->circuit = circuit;
	outerEdge2->circuit = circuit;

	facet1->circuit = circuit;
	facet2->circuit = circuit;
	if(isPrimarySegment){
		facet1->setFlag(1);
		facet2->setFlag(1);
	}

	edge0 = outerEdge1;
	edge1 = outerEdge2;
	edge2 = edge1->nextCircuitEdge;

	return true;
}

bool BurgersLoopBuilder::tryInsertOneCircuitEdge(
	InterfaceMesh::Edge*& edge0, 
	InterfaceMesh::Edge*& edge1, 
	bool isPrimarySegment
){
	assert(edge0 != edge1->oppositeEdge());

	InterfaceMesh::Face* facet = edge1->face();
	if(facet->circuit != nullptr) return false;

	InterfaceMesh::Edge* insertEdge1 = edge1->prevFaceEdge()->oppositeEdge();
	if(insertEdge1->circuit != nullptr) return false;

	InterfaceMesh::Edge* insertEdge2 = edge1->nextFaceEdge()->oppositeEdge();
	if(insertEdge2->circuit != nullptr) return false;

	if(!insertEdge1->tryClaimForCircuit()) return false;
	if(!insertEdge2->tryClaimForCircuit()){
		insertEdge1->releaseCircuitClaim();
		return false;
	}

	assert(insertEdge1->nextCircuitEdge == nullptr);
	assert(insertEdge2->nextCircuitEdge == nullptr);

	BurgersCircuit* circuit = edge0->circuit;
	
	insertEdge1->nextCircuitEdge = insertEdge2;
	insertEdge2->nextCircuitEdge = edge1->nextCircuitEdge;
	
	edge0->nextCircuitEdge = insertEdge1;
	if(edge0 == circuit->lastEdge){
		circuit->firstEdge = insertEdge1;
	}else if(edge1 == circuit->lastEdge){
		circuit->lastEdge = insertEdge2;
	}

	insertEdge1->circuit = circuit;
	insertEdge2->circuit = circuit;
	circuit->edgeCount++;

	assert(circuit->countEdges() == circuit->edgeCount);

	facet->circuit = circuit;
	if(isPrimarySegment) facet->setFlag(1);

	return true;
}

void BurgersLoopBuilder::appendLinePoint(DislocationNode& node){
	DislocationSegment& segment = *node.segment;
	assert(!segment.line.empty());

	int coreSize = node.circuit->edgeCount;

	const Point3& lastPoint = node.isForwardNode() ? segment.line.back() : segment.line.front();
	Point3 newPoint = lastPoint + cell().wrapVector(node.circuit->calculateCenter() - lastPoint);

	if(node.isForwardNode()){
		segment.line.push_back(newPoint);
		segment.coreSize.push_back(coreSize);
	}else{
		segment.line.insert(segment.line.begin(), newPoint);
		segment.coreSize.insert(segment.coreSize.begin(), coreSize);
	}

	node.circuit->numPreliminaryPoints++;

	if(_markCoreAtoms){
		markCoreCells(node, newPoint);
	}
}

void BurgersLoopBuilder::markCoreCells(DislocationNode& node, const Point3& center){
	if(node.circuit->firstEdge == nullptr) return;

	const DelaunayTessellation& tessellation = _mesh.tessellation();
	const bool capStored = !node.circuit->segmentMeshCap.empty();

	std::vector<TriangleIntersection::Triangle3> cap;
	cap.reserve(node.circuit->edgeCount);

	Box3 capBounds;
	auto* edge = node.circuit->firstEdge;
	do{
		const TriangleIntersection::Triangle3 triangle{
			center + cell().wrapVector(edge->vertex1()->pos() - center),
			center + cell().wrapVector(edge->vertex2()->pos() - center),
			center
		};
		for(const Point3& corner : triangle){
			capBounds.addPoint(corner);
		}
		cap.push_back(triangle);
		edge = edge->nextCircuitEdge;
	}while(edge != nullptr && edge != node.circuit->firstEdge);

	std::vector<DelaunayTessellation::CellHandle> candidates;
	_defectCells->queryOverlapping(capBounds, candidates);

	for(auto candidate : candidates){
		const int slot = tessellation.getUserField(candidate);
		if(slot < 0 || static_cast<std::size_t>(slot) >= _coreCellClaims.size()) continue;

		CoreCellClaim& claim = _coreCellClaims[slot];
		if(claim.node.load(std::memory_order_relaxed) != nullptr) continue;

		std::array<Point3, 4> tetrahedron;
		for(int corner = 0; corner < 4; ++corner){
			tetrahedron[corner] = tessellation.vertexPosition(tessellation.cellVertex(candidate, corner));
		}

		const bool intersects = std::ranges::any_of(cap, [&](const auto& triangle){
			for(int face = 0; face < 4; ++face){
				const TriangleIntersection::Triangle3 facet{
					tetrahedron[DelaunayTessellation::cellFacetVertexIndex(face, 0)],
					tetrahedron[DelaunayTessellation::cellFacetVertexIndex(face, 1)],
					tetrahedron[DelaunayTessellation::cellFacetVertexIndex(face, 2)]
				};
				if(TriangleIntersection::overlap(facet, triangle)) return true;
			}
			return false;
		});

		if(!intersects) continue;

		DislocationNode* unclaimed = nullptr;
		if(claim.node.compare_exchange_strong(unclaimed, &node)){
			claim.claimedWhileCapStored = capStored;
		}
	}
}

std::vector<int> BurgersLoopBuilder::coreAtomDislocationIds(std::size_t particleCount) const{
	std::vector<int> ids(particleCount, -1);
	if(_coreCellClaims.empty()) return ids;

	const DelaunayTessellation& tessellation = _mesh.tessellation();

	for(auto cellHandle : tessellation.cells()){
		const int slot = tessellation.getUserField(cellHandle);
		if(slot < 0 || static_cast<std::size_t>(slot) >= _coreCellClaims.size()) continue;

		const CoreCellClaim& claim = _coreCellClaims[slot];
		DislocationNode* node = claim.node.load(std::memory_order_relaxed);
		if(node == nullptr) continue;
		if(node->isDangling() && claim.claimedWhileCapStored) continue;

		const DislocationSegment* segment = node->segment;
		while(segment->replacedWith != nullptr){
			segment = segment->replacedWith;
		}

		for(int corner = 0; corner < 4; ++corner){
			const int particle = tessellation.vertexIndex(tessellation.cellVertex(cellHandle, corner));
			if(particle >= 0 && static_cast<std::size_t>(particle) < ids.size()){
				ids[particle] = segment->id;
			}
		}
	}

	return ids;
}

void BurgersLoopBuilder::circuitCircuitIntersection(
	InterfaceMesh::Edge* circuitAEdge1, 
	InterfaceMesh::Edge* circuitAEdge2, 
	InterfaceMesh::Edge* circuitBEdge1, 
	InterfaceMesh::Edge* circuitBEdge2, 
	int& goingOutside, 
	int& goingInside
){
	assert(circuitAEdge2->vertex1() == circuitBEdge2->vertex1());
	assert(circuitAEdge1->vertex2() == circuitBEdge2->vertex1());
	assert(circuitBEdge1->vertex2() == circuitBEdge2->vertex1());

	InterfaceMesh::Edge* edge = circuitBEdge2;
	bool contour1inside = false;
	bool contour2inside = false;
	int safetyCounter = 0;
	const int maxEdgeTraversal = 1000;
	
	for(;;){
		InterfaceMesh::Edge* oppositeEdge = edge->oppositeEdge();
		if(oppositeEdge == circuitBEdge1) break;
		if(edge != circuitBEdge2){
			if(oppositeEdge == circuitAEdge1) contour1inside = true;
			if(edge == circuitAEdge2) contour2inside = true;
		}

		edge = oppositeEdge->nextFaceEdge();
		if(edge->vertex1() != circuitBEdge2->vertex1() || edge == circuitBEdge2){
			break;
		}
		
		if(++safetyCounter > maxEdgeTraversal) break;
	}
	
	if(circuitAEdge2 != circuitBEdge2) {
	} else {
		contour2inside = false;
	}

	bool contour1outside = false;
	bool contour2outside = false;
	edge = circuitBEdge1;
	safetyCounter = 0;
	
	for(;;){
		InterfaceMesh::Edge* nextEdge = edge->nextFaceEdge();
		if(nextEdge == circuitBEdge2) break;
		InterfaceMesh::Edge* oppositeEdge = nextEdge->oppositeEdge();
		
		if(oppositeEdge->vertex2() != circuitBEdge2->vertex1()) {
			break;
		}
		
		edge = oppositeEdge;
		if(edge == circuitAEdge1) contour1outside = true;
		if(nextEdge == circuitAEdge2) contour2outside = true;
		
		if(++safetyCounter > 1000) break;
	}

	if(contour1outside && contour1inside) {
		contour1inside = false;
	}
	
	if(contour2outside && contour2inside) {
		contour2inside = false;
	}

	if(contour2outside == true && contour1outside == false){
		goingOutside += 1;
	}else if(contour2inside == true && contour1inside == false){
		goingInside += 1;
	}
}

size_t BurgersLoopBuilder::joinSegments(int maxCircuitLength){
	for(size_t nodeIndex = 0; nodeIndex < danglingNodes().size(); nodeIndex++){
		DislocationNode* node = danglingNodes()[nodeIndex];
		BurgersCircuit* circuit = node->circuit;
		assert(circuit->isDangling);

		InterfaceMesh::Edge* edge = circuit->firstEdge;
		do{
			assert(edge->circuit == circuit);
			BurgersCircuit* oppositeCircuit = edge->oppositeEdge()->circuit;
			if(oppositeCircuit == nullptr){
				assert(edge->oppositeEdge()->nextCircuitEdge == nullptr);

				createSecondarySegment(edge, circuit, maxCircuitLength);

				while(edge->oppositeEdge()->circuit == nullptr && edge != circuit->firstEdge){
					edge = edge->nextCircuitEdge;
				}
			}else{
				edge = edge->nextCircuitEdge;
			}
		}
		while(edge != circuit->firstEdge);
	}

	for(DislocationNode* node : danglingNodes()){
		BurgersCircuit* circuit = node->circuit;
		assert(circuit->isDangling);

		circuit->isCompletelyBlocked = true;
		InterfaceMesh::Edge* edge = circuit->firstEdge;
		do{
			assert(edge->circuit == circuit);
			BurgersCircuit* adjacentCircuit = edge->oppositeEdge()->circuit;
			if(adjacentCircuit == nullptr){
				circuit->isCompletelyBlocked = false;
				break;
			}else if(adjacentCircuit != circuit){
				assert(adjacentCircuit->isDangling);
				DislocationNode* adjacentNode = adjacentCircuit->dislocationNode;
				if(node->formsJunctionWith(adjacentNode) == false){
					node->connectNodes(adjacentNode);
				}
			}
			edge = edge->nextCircuitEdge;
		}while(edge != circuit->firstEdge);
	}

	size_t numJunctions = 0;

	for(DislocationNode* node : danglingNodes()){
		BurgersCircuit* circuit = node->circuit;

		if(circuit->isDangling == false) continue;

		if(!circuit->isCompletelyBlocked){
			node->dissolveJunction();
			continue;
		}

		if(node->junctionRing == node) continue;

		assert(node->segment->replacedWith == nullptr);

		Vector3 centerOfMassVector = Vector3::Zero();
		Point3 basePoint = node->position();
		int armCount = 1;
		bool allCircuitsCompletelyBlocked = true;
		DislocationNode* armNode = node->junctionRing;
		while(armNode != node){
			assert(armNode->segment->replacedWith == nullptr);
			assert(armNode->circuit->isDangling);
			if(armNode->circuit->isCompletelyBlocked == false) {
				allCircuitsCompletelyBlocked = false;
				break;
			}
			armCount++;
			centerOfMassVector += cell().wrapVector(armNode->position() - basePoint);
			armNode = armNode->junctionRing;
		}

		if(!allCircuitsCompletelyBlocked){
			node->dissolveJunction();
			continue;
		}

		assert(armCount >= 2);

		if(armCount >= 3){
			Point3 centerOfMass = basePoint + centerOfMassVector / armCount;

			armNode = node;
			do{
				armNode->circuit->isDangling = false;
				assert(armNode != armNode->junctionRing);

				std::vector<Point3>& line = armNode->segment->line;
				if(armNode->isForwardNode()){
					line.push_back(line.back() + cell().wrapVector(centerOfMass - line.back()));
					armNode->segment->coreSize.push_back(armNode->segment->coreSize.back());
				}else{
					line.insert(line.begin(), line.front() + cell().wrapVector(centerOfMass - line.front()));
					armNode->segment->coreSize.insert(armNode->segment->coreSize.begin(), armNode->segment->coreSize.front());
				}

				armNode->circuit->numPreliminaryPoints = 0;
				armNode = armNode->junctionRing;
			}while(armNode != node);
			numJunctions++;
		}else{
			DislocationNode* node1 = node;
			DislocationNode* node2 = node->junctionRing;
			assert(node1 != node2);
			assert(node2->junctionRing == node1);
			assert(node1->junctionRing == node2);

			BurgersCircuit* circuit1 = node1->circuit;
			BurgersCircuit* circuit2 = node2->circuit;
			circuit1->isDangling = false;
			circuit2->isDangling = false;
			circuit1->numPreliminaryPoints = 0;
			circuit2->numPreliminaryPoints = 0;

			if(node1->oppositeNode == node2){
				assert(node1->segment == node2->segment);
				DislocationSegment* loop = node1->segment;
				assert(loop->isClosedLoop());

				if(!cell().wrapVector(node1->position() - node2->position()).isZero(CA_ATOM_VECTOR_EPSILON)){
					loop->line.push_back(loop->line.back() + cell().wrapVector(loop->line.front() - loop->line.back()));
					assert(cell().wrapVector(node1->position() - node2->position()).isZero(CA_ATOM_VECTOR_EPSILON));
					loop->coreSize.push_back(loop->coreSize.back());
				}

				assert(loop->line.size() >= 3);
			}else{
				assert(node1->segment != node2->segment);

				DislocationNode* farEnd1 = node1->oppositeNode;
				DislocationNode* farEnd2 = node2->oppositeNode;
				DislocationSegment* segment1 = node1->segment;
				DislocationSegment* segment2 = node2->segment;
				if(node1->isBackwardNode()){
					segment1->nodes[1] = farEnd2;
					Vector3 shiftVector;
					if(node2->isBackwardNode()){
						shiftVector = calculateShiftVector(segment1->line.front(), segment2->line.front());
						segment1->line.insert(segment1->line.begin(), segment2->line.rbegin(), segment2->line.rend() - 1);
						segment1->coreSize.insert(segment1->coreSize.begin(), segment2->coreSize.rbegin(), segment2->coreSize.rend() - 1);
					}else{
						shiftVector = calculateShiftVector(segment1->line.front(), segment2->line.back());
						segment1->line.insert(segment1->line.begin(), segment2->line.begin(), segment2->line.end() - 1);
						segment1->coreSize.insert(segment1->coreSize.begin(), segment2->coreSize.begin(), segment2->coreSize.end() - 1);
					}

					if(shiftVector != Vector3::Zero()){
						for(auto p = segment1->line.begin(); p != segment1->line.begin() + segment2->line.size() - 1; ++p){
							*p -= shiftVector;
						}
					}
				}else{
					segment1->nodes[0] = farEnd2;
					Vector3 shiftVector;
					if(node2->isBackwardNode()){
						shiftVector = calculateShiftVector(segment1->line.back(), segment2->line.front());
						segment1->line.insert(segment1->line.end(), segment2->line.begin() + 1, segment2->line.end());
						segment1->coreSize.insert(segment1->coreSize.end(), segment2->coreSize.begin() + 1, segment2->coreSize.end());
					}else{
						shiftVector = calculateShiftVector(segment1->line.back(), segment2->line.back());
						segment1->line.insert(segment1->line.end(), segment2->line.rbegin() + 1, segment2->line.rend());
						segment1->coreSize.insert(segment1->coreSize.end(), segment2->coreSize.rbegin() + 1, segment2->coreSize.rend());
					}

					if(shiftVector != Vector3::Zero()){
						for(auto p = segment1->line.end() - segment2->line.size() + 1; p != segment1->line.end(); ++p){
							*p -= shiftVector;
						}
					}
				}
				
				farEnd2->segment = segment1;
				farEnd2->oppositeNode = farEnd1;
				farEnd1->oppositeNode = farEnd2;
				node1->oppositeNode = node2;
				node2->oppositeNode = node1;
				segment2->replacedWith = segment1;
				network().discardSegment(segment2);
			}
		}
	}

	_danglingNodes.erase(std::remove_if(_danglingNodes.begin(), _danglingNodes.end(),[](DislocationNode* node){
		return !node->isDangling();
	}), _danglingNodes.end());

	return numJunctions;
}

void BurgersLoopBuilder::createSecondarySegment(InterfaceMesh::Edge* firstEdge, BurgersCircuit* outerCircuit, int maxCircuitLength){
	assert(firstEdge->circuit == outerCircuit);

	int edgeCount = 1;
	Vector3 burgersVector = Vector3::Zero();
	Vector3 edgeSum = Vector3::Zero();
	Cluster* baseCluster = nullptr;
	Matrix3 frankRotation = Matrix3::Identity();
	int numCircuits = 1;
	InterfaceMesh::Edge* circuitStart = firstEdge->oppositeEdge();
	InterfaceMesh::Edge* circuitEnd = circuitStart;
	InterfaceMesh::Edge* edge = circuitStart;
	for(;;){
		for(;;){
			assert(edge->circuit == nullptr);
			InterfaceMesh::Edge* oppositeEdge = edge->oppositeEdge();
			InterfaceMesh::Face* oppositeFacet = oppositeEdge->face();
			InterfaceMesh::Edge* nextEdge = oppositeEdge->prevFaceEdge();
			assert(nextEdge->vertex2() == oppositeEdge->vertex1());
			assert(nextEdge->vertex2() == edge->vertex2());
			if(nextEdge->circuit != nullptr){
				if(nextEdge->circuit != outerCircuit){
					outerCircuit = nextEdge->circuit;
					numCircuits++;
				}
				
				edge = nextEdge->oppositeEdge();
				break;
			}
			edge = nextEdge;
		}

		circuitEnd->nextCircuitEdge = edge;
		edgeSum += edge->physicalVector;
		burgersVector += frankRotation * edge->clusterVector;
		if(!baseCluster) baseCluster = edge->clusterTransition->cluster1;
		if(!edge->clusterTransition->isSelfTransition()){
			frankRotation = frankRotation * edge->clusterTransition->reverse->tm;
		}

		if(edge == circuitStart) break;
		circuitEnd = edge;
		edgeCount++;

		if(edgeCount > maxCircuitLength) break;
	}

	const bool rejectNumCircuits = numCircuits == 1;
	const bool rejectTooLong = edgeCount > maxCircuitLength;
	const bool rejectZeroBurgers = burgersVector.isZero(CA_LATTICE_VECTOR_EPSILON);
	const bool rejectEdgeSum = edgeSum.isZero(CA_ATOM_VECTOR_EPSILON) == false;
	const bool rejectFrank = !frankRotation.equals(Matrix3::Identity(), CA_TRANSITION_MATRIX_EPSILON);
	if(rejectNumCircuits || rejectTooLong || rejectZeroBurgers || rejectEdgeSum || rejectFrank){

		edge = circuitStart;
		for(;;){
			assert(edge->circuit == nullptr);
			InterfaceMesh::Edge* nextEdge = edge->nextCircuitEdge;
			edge->nextCircuitEdge = nullptr;
			if(edge == circuitEnd) break;
			edge = nextEdge;
		}
		return;
	}
	assert(circuitStart != circuitEnd);

	BurgersCircuit* forwardCircuit = allocateCircuit();
	forwardCircuit->firstEdge = circuitStart;
	forwardCircuit->lastEdge = circuitEnd;
	forwardCircuit->edgeCount = edgeCount;
	edge = circuitStart;
	do{
		assert(edge->circuit == nullptr);
		edge->circuit = forwardCircuit;
		edge = edge->nextCircuitEdge;
	}while(edge != circuitStart);
	
	assert(forwardCircuit->countEdges() == forwardCircuit->edgeCount);

	createAndTraceSegment(ClusterVector(burgersVector, baseCluster), forwardCircuit, maxCircuitLength);
}

}
