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

// Build a watertight surface mesh over the interface where material properties may change
// (e.g, grain or cluster boundaries). We use a Delaunay Tessellation to generate tetrahedra,
// then carve out only those facets whose endpoints belongs to "compatible" clusters (as determined
// by the elastic mapping). Faces that bridge incompatible clusters get omitted, leaving
// a manifold of boundary faces.
void InterfaceMesh::createMesh(double maxNeighborDist, double alphaScale){
    _isCompletelyGood = true;
    _isCompletelyBad  = true;

    // Classify each tetrahedron. Return 1 to keep its faces if its interior is
    // "elastic-compatible" (no large strain or cell-size mismatch), otherwise 0.
    auto tetraRegion = [&](auto cell) -> unsigned {
        if(!elasticMapping().isElasticMappingCompatible(cell)){
            // Found at least one bad tetrahedra
			_isCompletelyGood = false;

            // Omit all faces of this cell
			return 0;
		}

        // At least one good tetrahedra exists.
		_isCompletelyBad = false;

        // Keep its faces
		return 1;
    };

    // Initialize each triangular face by copying its three vertex positions,
    // computing the physical (Cartesian) edge vectors, and looking up the
    // "ideal" cluster-to-cluster displacement and rotation that went into 
    // generating those vertices. Also detect if any PBC wrapping would
    // imply the simulation cell is too small for the neighbor distance.
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
            // Compute the actual displacement vector between the two vertices
            e->physicalVector = pos[ni] - pos[i];

            // Query the elastic mapping to find the cluster-to-cluster shift and rotation
            // that ideally produced this edge. We store both the local Burgers-vector-like "clusterVector"
            // and the small rotation "clusterTransition" that maps between the two grain orientations.
            std::tie(e->clusterVector, e->clusterTransition) = elasticMapping().getEdgeClusterVector(vIdx[i], vIdx[ni]);
            e = e->nextFaceEdge();
        }
    };

    // We pad the ghost layer size by several neighbor distances to ensure
    // all relevant tetrahedra are built across the periodic box. The manifold
    // helper will call back into our lambdas to decide which test and faces to keep.
    double alpha = alphaScale * maxNeighborDist;
    ManifoldConstructionHelper<InterfaceMesh> helper{
        tessellation(),
        *this,
        alpha,
        structureAnalysis().context().positions
    };

    // Build the faces and topology. If any step fails, bail out.
    if(!helper.construct(tetraRegion, prepareFace)){
        throw std::runtime_error("Error building the faces and topology.");
	}
    makeManifold();
}

// After tracing dislocation circuits on the interface mesh, drop the facets the tracing
// claimed for a dislocation and keep the rest, so what remains is the surface bounding the
// non-crystalline regions with a hole where each dislocation line passes through. Those holes
// are then capped, one fan of triangles per dangling Burgers circuit, making the result a
// closed surface again. The output is a separate half-edge mesh; this one is left untouched.
bool InterfaceMesh::generateDefectMesh(
    BurgersLoopBuilder const& tracer,
    InterfaceMeshTopology& outMesh
) const{
    // Copy every interface vertex, keyed by its original index rather than by iteration
    // order, so nothing downstream depends on the two meshes staying positionally aligned.
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

    // A facet belongs to a dislocation — and so is left out of the defect mesh — when it sits
    // on a Burgers circuit that either was swept by a primary segment (flag 1) or has already
    // been closed into a segment (!isDangling). Facets on still-dangling circuits are the rims
    // of the holes we cap further down, so they stay.
    std::vector<Face*> faceMap(faces().size(), nullptr);
    outMesh.reserveFaces(faces().size());

    for(auto* f : faces()){
        if(f->circuit && (f->testFlag(1) || !f->circuit->isDangling)) continue;
        if(!f->edges()) continue;

        auto* nf = outMesh.createFace();
        if(f->index() >= 0 && static_cast<size_t>(f->index()) < faceMap.size()){
            faceMap[f->index()] = nf;
        }

        // Walk its three edges in order and add half-edges
        auto* e = f->edges();
        auto* start = e;
        do{
            auto* v1 = mapped(e->vertex1());
            auto* v2 = mapped(e->vertex2());
            if(v1 && v2) outMesh.createEdge(v1, v2, nf);
            e = e->nextFaceEdge();
        }while(e != start);
    }

    // Link opposite half-edges between the pairs of facets that both survived. An edge whose
    // neighbour was dropped stays unpaired on purpose: that is a hole rim, closed by a cap.
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

    // Cap each dangling circuit's hole: one vertex at the dislocation line's endpoint, then a
    // triangle fan from it to every edge of the stored circuit loop. The triangle is wound
    // opposite to the rim edge so its normal agrees with the surrounding surface.
    for(auto* dn : tracer.danglingNodes()){
        if(!dn) continue;
        auto* c = dn->circuit;
        if(!c || c->segmentMeshCap.size() < 2) continue;

        // position() dereferences segment (via isForwardNode()) and then line.front()/back(),
        // so a null segment or a line that finishDislocationSegments() trimmed down to nothing
        // would be undefined behaviour. Such a node gets no cap, and the mesh comes back open
        // from connectOppositeHalfedges() below rather than the caller getting a crash.
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

    // Pair up whatever is still unpaired — the cap edges against the rims they close. The
    // return value is the watertightness check the old code left as a commented-out assert:
    // false means a hole survived, which the caller reports rather than aborting on.
    return outMesh.connectOppositeHalfedges();
}

}
