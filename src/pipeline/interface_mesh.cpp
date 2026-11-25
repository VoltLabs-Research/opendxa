#include <volt/pipeline/interface_mesh.h>
#include <volt/pipeline/burgers_loop_builder.h>
#include <volt/helpers/manifold_construction_helper.h>
#include <algorithm>
#include <array>
#include <numeric>
#include <cassert>

namespace Volt{

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

// Links the copy of an interface edge to the copy of the edge's opposite, when the opposite's face was kept too.
void InterfaceMesh::linkCopiedOpposite(
    const InterfaceMeshTopology::Edge* edge,
    InterfaceMeshTopology::Edge* copy,
    const std::vector<InterfaceMeshTopology::Face*>& faceMap
){
    if(!edge->oppositeEdge() || copy->oppositeEdge()) return;

    InterfaceMeshTopology::Face* oppositeCopy = faceMap[edge->oppositeEdge()->face()->index()];
    if(!oppositeCopy) return;

    InterfaceMeshTopology::Edge* candidate = oppositeCopy->edges();
    do{
        if(!candidate->oppositeEdge() && candidate->vertex1() == copy->vertex2() && candidate->vertex2() == copy->vertex1()){
            copy->linkToOppositeEdge(candidate);
            return;
        }
        candidate = candidate->nextFaceEdge();
    }while(candidate != oppositeCopy->edges());
}

void InterfaceMesh::generateDefectMesh(
    BurgersLoopBuilder const& tracer,
    InterfaceMeshTopology& outMesh
) const{
    std::vector<Vertex*> vertexMap(vertexCount(), nullptr);
    outMesh.reserveVertices(vertexCount());
    for(auto* v : vertices()){
        vertexMap[v->index()] = outMesh.createVertex(v->pos());
    }

    auto mapped = [&](Vertex* v){
        return vertexMap[v->index()];
    };

    std::vector<Face*> faceMap(faces().size(), nullptr);
    outMesh.reserveFaces(faces().size());

    for(auto* f : faces()){
        if(f->circuit && (f->testFlag(InterfaceMeshFace::kSweptBeforeCapFlag) || !f->circuit->isDangling)) continue;
        if(!f->edges()) continue;

        auto* nf = outMesh.createFace();
        faceMap[f->index()] = nf;

        auto* e = f->edges();
        do{
            outMesh.createEdge(mapped(e->vertex1()), mapped(e->vertex2()), nf);
            e = e->nextFaceEdge();
        }while(e != f->edges());
    }

    for(auto* f : faces()){
        auto* nf = faceMap[f->index()];
        if(!nf) continue;

        auto* e = f->edges();
        auto* ne = nf->edges();
        do{
            linkCopiedOpposite(e, ne, faceMap);
            e = e->nextFaceEdge();
            ne = ne->nextFaceEdge();
        }while(e != f->edges());
    }

    for(auto* dn : tracer.danglingNodes()){
        const BurgersCircuit* c = dn->circuit;
        if(c->segmentMeshCap.size() < 2 || dn->segment->line.empty()) continue;

        auto* capV = outMesh.createVertex(dn->position());
        for(auto* me : c->segmentMeshCap){
            auto* v1 = mapped(me->vertex2());
            auto* v2 = mapped(me->vertex1());

            auto* nf = outMesh.createFace();
            outMesh.createEdge(v1, v2, nf);
            outMesh.createEdge(v2, capV, nf);
            outMesh.createEdge(capV, v1, nf);
        }
    }

    outMesh.connectOppositeHalfedges();
}

}
