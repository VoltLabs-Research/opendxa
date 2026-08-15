#include <volt/pipeline/interface_mesh.h>
#include <volt/pipeline/burgers_loop_builder.h>
#include <volt/helpers/manifold_construction_helper.h>
#include <algorithm>
#include <array>
#include <numeric>
#include <cassert>
#include <set> 
#include <boost/dynamic_bitset.hpp>

namespace Volt{

void InterfaceMesh::makeManifold(){
    auto original_vertices = vertices(); 

    for(auto* vertex : original_vertices){
        if(vertex->numEdges() < 3) continue;

        std::set<Edge*> visited_edges;
        Edge* start_edge = vertex->edges();
        Edge* current_edge = start_edge;
        do{
            visited_edges.insert(current_edge);
            current_edge = current_edge->oppositeEdge()->nextFaceEdge();
        }while(current_edge != start_edge);

        if(visited_edges.size() == vertex->numEdges()){
            continue;
        }

        while(visited_edges.size() < vertex->numEdges()){
            Vertex* new_vertex = createVertex(vertex->pos());
            Edge* fan_start_edge = nullptr;
            for(Edge* e = vertex->edges(); e != nullptr; e = e->nextVertexEdge()){
                if(visited_edges.find(e) == visited_edges.end()){
                    fan_start_edge = e;
                    break;
                }
            }
            
            if(fan_start_edge == nullptr) break;

            std::vector<Edge*> fan_to_transfer;
            current_edge = fan_start_edge;
            do{
                fan_to_transfer.push_back(current_edge);
                visited_edges.insert(current_edge);
                current_edge = current_edge->oppositeEdge()->nextFaceEdge();
            }while(current_edge != fan_start_edge);

            for(Edge* edge_to_move : fan_to_transfer){
                vertex->transferEdgeToVertex(edge_to_move, new_vertex);
            }
        }
    }
}

void InterfaceMesh::createMesh(double maxNeighborDist, double alphaScale){
    _isCompletelyGood = true;
    _isCompletelyBad  = true;

    auto tetraRegion = [&](auto cell) -> unsigned {
        if(!elasticMapping().isElasticMappingCompatible(cell)){
			_isCompletelyGood = false;

			return 0;
		}

		_isCompletelyBad = false;

		return 1;
    };

    auto prepareFace = [&](
		Face* face,
		std::array<int, 3> const& vIdx,
		std::array<decltype(tessellation().cellVertex(0,0)),3> const& vH,
		auto cell
	){
        std::array<Point3, 3> pos = {
            tessellation().vertexPosition(vH[0]),
            tessellation().vertexPosition(vH[1]),
            tessellation().vertexPosition(vH[2])
        };

        auto* e = face->edges();
        for(int i = 0; i < 3; ++i){
            int ni = (i + 1) % 3;
            e->physicalVector = pos[ni] - pos[i];

            std::tie(e->clusterVector, e->clusterTransition) = elasticMapping().getEdgeClusterVector(vIdx[i], vIdx[ni]);
            e = e->nextFaceEdge();
        }
    };

    double alpha = alphaScale * maxNeighborDist;
    ManifoldConstructionHelper<InterfaceMesh> helper{
        tessellation(),
        *this,
        alpha,
        structureAnalysis().context().positions
    };

    if(!helper.construct(tetraRegion, prepareFace)){
        throw std::runtime_error("Error building the faces and topology.");
	}
    makeManifold();
}

bool InterfaceMesh::generateDefectMesh(
    BurgersLoopBuilder const& tracer,
    InterfaceMeshTopology& outMesh
) const{
    std::vector<Vertex*> vertexMap(vertexCount(), nullptr);
    outMesh.reserveVertices(vertexCount());
    for(auto* v : vertices()){
        auto* nv = outMesh.createVertex(v->pos());
        if(v->index() >= 0 && static_cast<size_t>(v->index()) < vertexMap.size()){
            vertexMap[v->index()] = nv;
        }
    }

    auto mapped = [&](Vertex* v) -> Vertex*{
        const int idx = v ? v->index() : -1;
        return (idx >= 0 && static_cast<size_t>(idx) < vertexMap.size()) ? vertexMap[idx] : nullptr;
    };

    std::vector<Face*> faceMap(faces().size(), nullptr);
    outMesh.reserveFaces(faces().size());

    for(auto* f : faces()){
        if(f->circuit && (f->testFlag(1) || !f->circuit->isDangling)) continue;
        if(!f->edges()) continue;

        auto* nf = outMesh.createFace();
        if(f->index() >= 0 && static_cast<size_t>(f->index()) < faceMap.size()){
            faceMap[f->index()] = nf;
        }

        auto* e = f->edges();
        auto* start = e;
        do{
            auto* v1 = mapped(e->vertex1());
            auto* v2 = mapped(e->vertex2());
            if(v1 && v2) outMesh.createEdge(v1, v2, nf);
            e = e->nextFaceEdge();
        }while(e != start);
    }

    for(auto* of : faces()){
        if(of->index() < 0 || static_cast<size_t>(of->index()) >= faceMap.size()) continue;
        auto* nf = faceMap[of->index()];
        if(!nf || !of->edges() || !nf->edges()) continue;

        auto* eo = of->edges();
        auto* en = nf->edges();
        auto* startO = eo;
        do{
            auto* oppEo = eo->oppositeEdge();
            if(oppEo && oppEo->face() && !en->oppositeEdge()){
                const int oppIndex = oppEo->face()->index();
                auto* oppNF = (oppIndex >= 0 && static_cast<size_t>(oppIndex) < faceMap.size())
                    ? faceMap[oppIndex]
                    : nullptr;
                if(oppNF && oppNF->edges()){
                    auto* ec = oppNF->edges();
                    auto* startC = ec;
                    do{
                        if(!ec->oppositeEdge() && ec->vertex1() == en->vertex2() && ec->vertex2() == en->vertex1()){
                            en->linkToOppositeEdge(ec);
                            break;
                        }

                        ec = ec->nextFaceEdge();
                    }while(ec != startC);
                }
            }

            eo = eo->nextFaceEdge();
            en = en->nextFaceEdge();
        }while(eo != startO);
    }

    for(auto* dn : tracer.danglingNodes()){
        if(!dn) continue;
        auto* c = dn->circuit;
        if(!c || c->segmentMeshCap.size() < 2) continue;

        if(!dn->segment || dn->segment->line.empty()) continue;

        auto* capV = outMesh.createVertex(dn->position());

        for(auto* me : c->segmentMeshCap){
            if(!me) continue;
            auto* v1 = mapped(me->vertex2());
            auto* v2 = mapped(me->vertex1());
            if(!v1 || !v2) continue;

            auto* nf = outMesh.createFace();
            outMesh.createEdge(v1, v2, nf);
            outMesh.createEdge(v2, capV, nf);
            outMesh.createEdge(capV, v1, nf);
        }
    }

    return outMesh.connectOppositeHalfedges();
}

}
