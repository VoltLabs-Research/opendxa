#pragma once

#include <volt/core/volt.h>
#include <volt/utilities/memory_pool.h>
#include <tbb/spin_mutex.h>
#include <mutex>
#include <set>
#include <vector>

namespace Volt{

struct EmptyHalfEdgeMeshStruct{};

template <
    class EdgeBase = EmptyHalfEdgeMeshStruct,
    class FaceBase = EmptyHalfEdgeMeshStruct,
    class VertexBase = EmptyHalfEdgeMeshStruct
>

class HalfEdgeMesh{
public:
    class Vertex;
    class Face;

    class Edge : public EdgeBase{
    public:
        Edge(Vertex* v2, Face* f) noexcept
            : _oppositeEdge(nullptr)
            , _vertex2(v2)
            , _face(f)
            , _nextVertexEdge(nullptr)
            , _nextFaceEdge(nullptr)
            , _prevFaceEdge(nullptr){}

        [[nodiscard]] Vertex* vertex1() const noexcept{
			return _prevFaceEdge->_vertex2;
		}

        [[nodiscard]] Vertex* vertex2() const noexcept{
			return _vertex2;
		}

        [[nodiscard]] Face* face() const noexcept{
			return _face;
		}

        [[nodiscard]] Edge* nextVertexEdge() const noexcept{
			return _nextVertexEdge;
		}

        [[nodiscard]] Edge* nextFaceEdge() const noexcept{
			return _nextFaceEdge;
		}

        [[nodiscard]] Edge* prevFaceEdge() const noexcept{
			return _prevFaceEdge;
		}

        [[nodiscard]] Edge* oppositeEdge() const noexcept{
			return _oppositeEdge;
		}

        void linkToOppositeEdge(Edge* other) noexcept{
            if(!other) return;
            _oppositeEdge = other;
            other->_oppositeEdge = this;
        }

    private:
        Edge* _oppositeEdge;
        Vertex* _vertex2;
        Face* _face;
        Edge* _nextVertexEdge;
        Edge* _nextFaceEdge;
        Edge* _prevFaceEdge;

        friend class HalfEdgeMesh;
        friend Vertex;
    };

    class Vertex : public VertexBase{
    public:
        Vertex(const Point3& p, int idx = -1) noexcept
            : _pos(p)
            , _edges(nullptr)
            , _numEdges(0)
            , _index(idx){}

        [[nodiscard]] Edge* edges() const noexcept{
			return _edges;
		}

        [[nodiscard]] const Point3& pos() const noexcept{
			return _pos;
		}

        [[nodiscard]] Point3& pos() noexcept{
			return _pos;
		}

        [[nodiscard]] int index() const noexcept{
			return _index;
		}

        [[nodiscard]] std::size_t numEdges() const noexcept{
			return _numEdges;
		}

        void transferEdgeToVertex(Edge* e, Vertex* v) noexcept{
            removeEdge(e);
            v->addEdge(e);
            e->oppositeEdge()->_vertex2 = v;
        }

    private:
        void addEdge(Edge* e) noexcept{
            tbb::spin_mutex::scoped_lock lock(_mutex);
            e->_nextVertexEdge = _edges;
            _edges = e;
            ++_numEdges;
        }

        void removeEdge(Edge* e) noexcept{
            tbb::spin_mutex::scoped_lock lock(_mutex);
            --_numEdges;
            if(e == _edges){
                _edges = e->nextVertexEdge();
                e->_nextVertexEdge = nullptr;
                return;
            }

            for(Edge* cur = _edges; cur; cur = cur->nextVertexEdge()){
                if(cur->nextVertexEdge() == e){
                    cur->_nextVertexEdge = e->_nextVertexEdge;
                    e->_nextVertexEdge = nullptr;
                    return;
                }
            }

            assert(false);
        }

        Point3 _pos;
        Edge* _edges;
        size_t _numEdges;
        int _index;
        mutable tbb::spin_mutex _mutex;

        friend class HalfEdgeMesh;
    };

    class Face : public FaceBase{
    public:
        [[nodiscard]] Edge* edges() const noexcept{
			return _edges;
		}

        [[nodiscard]] int index() const noexcept{
			return _index;
		}

        [[nodiscard]] bool testFlag(unsigned int f) const noexcept{
			return _flags & f;
		}

        void setFlag(unsigned int f) noexcept{
			_flags |= f;
		}

        void clearFlag(unsigned int f) const noexcept{
			_flags &= ~f;
		}

        [[nodiscard]] Edge* findEdge(Vertex* v1, Vertex* v2) const noexcept{
            if(!edges()) return nullptr;

            Edge* e = edges();
            do{
                if(e->vertex1() == v1 && e->vertex2() == v2){
					return e;
				}
                e = e->nextFaceEdge();
            }while(e != edges());

            return nullptr;
        }

    protected:
        // Protected so that faces are only made by createFace(), which numbers them.
        Face(int idx = -1) noexcept
            : _edges(nullptr)
            , _index(idx)
            , _flags(0){}

    private:
        Edge* _edges;
        int _index;
        mutable unsigned int _flags;

        friend class HalfEdgeMesh;
    };

public:
    HalfEdgeMesh() = default;
    HalfEdgeMesh(const HalfEdgeMesh&) = delete;
    HalfEdgeMesh& operator=(const HalfEdgeMesh&) = delete;

    void clear() noexcept{
        _vertices.clear();
        _faces.clear();
        _vertexPool.clear();
        _edgePool.clear();
        _facePool.clear();
    }

    [[nodiscard]] const std::vector<Vertex*>& vertices() const noexcept{
		return _vertices;
	}

    [[nodiscard]] const std::vector<Face*>& faces() const noexcept{
		return _faces;
	}

    [[nodiscard]] std::size_t vertexCount() const noexcept{
		return _vertices.size();
	}

    [[nodiscard]] std::size_t faceCount() const noexcept{
		return _faces.size();
	}

    [[nodiscard]] Vertex* vertex(int idx) const noexcept{
        assert(idx >= 0 && idx < static_cast<int>(vertexCount()));
        return _vertices[idx];
    }

    [[nodiscard]] Face* face(int idx) const noexcept{
        assert(idx >= 0 && idx < static_cast<int>(faceCount()));
        return _faces[idx];
    }

    void reserveVertices(std::size_t n){
		_vertices.reserve(n);
	}

    void reserveFaces(std::size_t n){
		_faces.reserve(n);
	}

    Vertex* createVertex(const Point3& p){
        tbb::spin_mutex::scoped_lock lock(_mutex);
        auto* v = _vertexPool.construct(p, static_cast<int>(vertexCount()));
        _vertices.push_back(v);
        return v;
    }

    Face* createFace(){
        tbb::spin_mutex::scoped_lock lock(_mutex);
        auto* f = _facePool.construct(static_cast<int>(faceCount()));
        _faces.push_back(f);
        return f;
    }

    Edge* createEdge(Vertex* v1, Vertex* v2, Face* f){
        tbb::spin_mutex::scoped_lock lock(_mutex);
        auto* e = _edgePool.construct(v2, f);
        v1->addEdge(e);

        if(f->_edges){
            e->_nextFaceEdge = f->_edges;
            e->_prevFaceEdge = f->_edges->_prevFaceEdge;
            f->_edges->_prevFaceEdge->_nextFaceEdge = e;
            f->_edges->_prevFaceEdge = e;
        }else{
            e->_nextFaceEdge = e;
            e->_prevFaceEdge = e;
            f->_edges = e;
        }

        return e;
    }

    template<typename It>
    Face* createFace(It b, It e){
        assert(std::distance(b, e) >= 2);
        Face* f = createFace();
        It v1 = b;
        for (It v2 = std::next(b); v2 != e; ++v1, ++v2){
            createEdge(*v1, *v2, f);
		}
        createEdge(*v1, *b, f);
        return f;
    }
	
    Face* createFace(std::initializer_list<Vertex*> v){
        return createFace(v.begin(), v.end());
    }

    void clearFaceFlag(unsigned int flag) const noexcept{
        for (Face* f : _faces){
            f->clearFlag(flag);
		}
    }

    // Whether every half-edge has its opposite, so the mesh has no hole.
    [[nodiscard]] bool isClosed() const noexcept{
        for(const Vertex* vertex : _vertices){
            for(const Edge* edge = vertex->edges(); edge; edge = edge->nextVertexEdge()){
                if(!edge->oppositeEdge()) return false;
            }
        }
        return true;
    }

    bool connectOppositeHalfedges(){
        bool closed = true;
        for(Vertex* v1 : _vertices){
            for(Edge* e = v1->edges(); e; e = e->nextVertexEdge()){
                if(e->oppositeEdge()) continue;
                for(Edge* oe = e->vertex2()->edges(); oe; oe = oe->nextVertexEdge()){
                    if(!oe->oppositeEdge() && oe->vertex2() == v1){
                        e->linkToOppositeEdge(oe);
                        break;
                    }
                }
                if(!e->oppositeEdge()) closed = false;
            }
        }

        return closed;
    }

    // Splits every vertex whose faces form more than one fan into one vertex per fan, so that the surface becomes a
    // manifold. Each new vertex takes the position and the attributes of the vertex it splits off.
    void makeManifold(){
        const std::vector<Vertex*> originalVertices = _vertices;
        for(Vertex* vertex : originalVertices){
            if(vertex->numEdges() < 3) continue;

            std::set<Edge*> visitedEdges;
            Edge* startEdge = vertex->edges();
            Edge* currentEdge = startEdge;
            do{
                visitedEdges.insert(currentEdge);
                currentEdge = currentEdge->oppositeEdge()->nextFaceEdge();
            }while(currentEdge != startEdge);

            while(visitedEdges.size() < vertex->numEdges()){
                Vertex* newVertex = createVertex(vertex->pos());
                static_cast<VertexBase&>(*newVertex) = static_cast<const VertexBase&>(*vertex);
                Edge* fanStartEdge = nullptr;
                for(Edge* e = vertex->edges(); e != nullptr; e = e->nextVertexEdge()){
                    if(visitedEdges.find(e) == visitedEdges.end()){
                        fanStartEdge = e;
                        break;
                    }
                }
                if(fanStartEdge == nullptr) break;

                std::vector<Edge*> fanToTransfer;
                currentEdge = fanStartEdge;
                do{
                    fanToTransfer.push_back(currentEdge);
                    visitedEdges.insert(currentEdge);
                    currentEdge = currentEdge->oppositeEdge()->nextFaceEdge();
                }while(currentEdge != fanStartEdge);

                for(Edge* edgeToMove : fanToTransfer){
                    vertex->transferEdgeToVertex(edgeToMove, newVertex);
                }
            }
        }
    }

private:
    // Gives the pool access to Face's protected constructor.
    class InternalFace : public Face{
    public:
        InternalFace(int idx) noexcept : Face(idx){}
    };

    std::vector<Vertex*> _vertices;
    MemoryPool<Vertex> _vertexPool;
    MemoryPool<Edge> _edgePool;
    std::vector<Face*> _faces;
    MemoryPool<InternalFace> _facePool;
    mutable tbb::spin_mutex _mutex;
};

}
