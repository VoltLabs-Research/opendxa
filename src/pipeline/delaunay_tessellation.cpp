#include <volt/core/volt.h>
#include <volt/pipeline/delaunay_tessellation.h>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <spdlog/spdlog.h>
#include <numeric>
#include <cstdint>
#include <algorithm>
#include <limits>
#include <mutex>

namespace Volt{

std::uint64_t DelaunayTessellation::splitMix64(std::uint64_t z){
	z += 0x9E3779B97F4A7C15ull;
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
	return z ^ (z >> 31);
}

void DelaunayTessellation::initializeGeogram(){
	static std::mutex geogramMutex;
	std::lock_guard<std::mutex> lock(geogramMutex);
	GEO::initialize(GEO::GEOGRAM_NO_HANDLER);
	GEO::set_assert_mode(GEO::ASSERT_ABORT);
}


void DelaunayTessellation::generateTessellation(
	const SimulationCell& simCell,
	const Point3* positions,
	size_t numPoints,
	double ghostLayerSize,
	bool coverDomainWithFiniteTets,
	const int* selectedPoints
){
	initializeGeogram();
	insertPrimaryPoints(simCell, positions, numPoints, selectedPoints);
	insertGhostImages(simCell, ghostLayerSize);
	if(coverDomainWithFiniteTets){
		addBoundingCorners(simCell, ghostLayerSize);
	}

	spdlog::info("  Geogram Delaunay: inserting {} points ({} primary + {} ghost)",
		_pointData.size(), _primaryVertexCount, _pointData.size() - _primaryVertexCount);
	triangulate();
	classifyCells();
}

void DelaunayTessellation::insertPrimaryPoints(
	const SimulationCell& simCell,
	const Point3* positions,
	size_t numPoints,
	const int* selectedPoints
){
	const double lengthScale = (simCell.matrix().column(0) + simCell.matrix().column(1) + simCell.matrix().column(2)).length();
	const double epsilon = kRelativeJitter * lengthScale;

	auto jitter = [epsilon](std::size_t atomIndex, unsigned axis) -> double {
		const std::uint64_t z = splitMix64(static_cast<std::uint64_t>(atomIndex) * 3ull + axis);
		// The top 53 bits as a double in [0, 1).
		const double unit = static_cast<double>(z >> 11) * (1.0 / 9007199254740992.0);
		return epsilon * (2.0 * unit - 1.0);
	};
	auto jitteredPoint = [&](std::size_t atomIndex, const Point3& position){
		const Point3 wp = simCell.wrapPoint(position);
		return Point3(
			wp.x() + jitter(atomIndex, 0),
			wp.y() + jitter(atomIndex, 1),
			wp.z() + jitter(atomIndex, 2)
		);
	};

	if(!selectedPoints){
		_pointData.resize(numPoints);
		_particleIndices.resize(numPoints);
		tbb::parallel_for(tbb::blocked_range<size_t>(0, numPoints, 4096),
			[&](const tbb::blocked_range<size_t>& r){
			for(size_t i = r.begin(); i < r.end(); ++i){
				_pointData[i] = jitteredPoint(i, positions[i]);
				_particleIndices[i] = i;
			}
		});
	}else{
		_particleIndices.clear();
		_pointData.clear();
		_particleIndices.reserve(numPoints);
		_pointData.reserve(numPoints);
		for(size_t i = 0; i < numPoints; i++){
			if(!selectedPoints[i]) continue;
			_pointData.push_back(jitteredPoint(i, positions[i]));
			_particleIndices.push_back(i);
		}
	}
	_primaryVertexCount = _particleIndices.size();
}

// Adds, for every periodic image of the cell within the ghost layer, the images of the primary points that fall
// inside the cell widened by the ghost layer. Images are counted first so that each one can be written in parallel
// at a fixed offset, in the same order as a serial pass.
void DelaunayTessellation::insertGhostImages(const SimulationCell& simCell, double ghostLayerSize){
	Vector3I stencilCount;
	double cuts[3][2];
	Vector3 cellNormals[3];
	for(size_t dim = 0; dim < 3; dim++){
		cellNormals[dim] = simCell.cellNormalVector(dim);
		cuts[dim][0] = cellNormals[dim].dot(simCell.reducedToAbsolute(Point3(0,0,0)) - Point3::Origin());
		cuts[dim][1] = cellNormals[dim].dot(simCell.reducedToAbsolute(Point3(1,1,1)) - Point3::Origin());
		stencilCount[dim] = simCell.hasPbc(dim)
			? (int)ceil(ghostLayerSize / simCell.matrix().column(dim).dot(cellNormals[dim]))
			: 0;
		cuts[dim][0] -= ghostLayerSize;
		cuts[dim][1] += ghostLayerSize;
	}

	auto isInsideGhostLayer = [&](const Point3& image){
		for(size_t dim = 0; dim < 3; ++dim){
			if(!simCell.hasPbc(dim)) continue;
			const double d = cellNormals[dim].dot(image - Point3::Origin());
			if(d < cuts[dim][0] || d > cuts[dim][1]) return false;
		}
		return true;
	};

	std::vector<Vector3> ghostShifts;
	ghostShifts.reserve(static_cast<size_t>(
		(2 * stencilCount[0] + 1) * (2 * stencilCount[1] + 1) * (2 * stencilCount[2] + 1) - 1));
	for(int ix = -stencilCount[0]; ix <= +stencilCount[0]; ix++)
		for(int iy = -stencilCount[1]; iy <= +stencilCount[1]; iy++)
			for(int iz = -stencilCount[2]; iz <= +stencilCount[2]; iz++){
				if(ix == 0 && iy == 0 && iz == 0) continue;
				ghostShifts.push_back(simCell.reducedToAbsolute(Vector3(ix, iy, iz)));
			}

	std::vector<size_t> ghostCounts(ghostShifts.size(), 0);
	tbb::parallel_for(tbb::blocked_range<size_t>(0, ghostShifts.size(), 8),
		[&](const tbb::blocked_range<size_t>& r){
		for(size_t shiftIdx = r.begin(); shiftIdx < r.end(); ++shiftIdx){
			size_t count = 0;
			for(size_t vertexIndex = 0; vertexIndex < _primaryVertexCount; ++vertexIndex){
				if(isInsideGhostLayer(_pointData[vertexIndex] + ghostShifts[shiftIdx])) ++count;
			}
			ghostCounts[shiftIdx] = count;
		}
	});

	std::vector<size_t> ghostOffsets(ghostCounts.size(), 0);
	size_t totalGhostPoints = 0;
	for(size_t i = 0; i < ghostCounts.size(); ++i){
		ghostOffsets[i] = totalGhostPoints;
		totalGhostPoints += ghostCounts[i];
	}

	const size_t basePointCount = _pointData.size();
	_pointData.resize(basePointCount + totalGhostPoints);
	_particleIndices.resize(basePointCount + totalGhostPoints);

	tbb::parallel_for(tbb::blocked_range<size_t>(0, ghostShifts.size(), 8),
		[&](const tbb::blocked_range<size_t>& r){
		for(size_t shiftIdx = r.begin(); shiftIdx < r.end(); ++shiftIdx){
			size_t writeIndex = basePointCount + ghostOffsets[shiftIdx];
			for(size_t vertexIndex = 0; vertexIndex < _primaryVertexCount; ++vertexIndex){
				const Point3 image = _pointData[vertexIndex] + ghostShifts[shiftIdx];
				if(!isInsideGhostLayer(image)) continue;
				_pointData[writeIndex] = image;
				_particleIndices[writeIndex] = _particleIndices[vertexIndex];
				++writeIndex;
			}
		}
	});
}

// The corners of a box around every point make the convex hull cover the whole domain with finite cells.
// They stand for no particle.
void DelaunayTessellation::addBoundingCorners(const SimulationCell& simCell, double ghostLayerSize){
	Box3 bb = Box3(Point3(0, 0, 0), Point3(1, 1, 1)).transformed(simCell.matrix());
	bb.addPoints(_pointData.data(), _pointData.size());
	bb = bb.padBox(ghostLayerSize);
	for(int corner = 0; corner < 8; corner++){
		_pointData.push_back(bb[corner]);
		_particleIndices.push_back(std::numeric_limits<size_t>::max());
	}
}

// Serial Geogram (BDEL) even when threads are available. The parallel backend (PDEL) builds the same tetrahedra but
// numbers them in an order that changes from run to run, and the interface mesh and the Burgers circuit search follow
// the cell order: with PDEL, the same input could give a different dislocation network on every run.
void DelaunayTessellation::triangulate(){
	_dt = GEO::Delaunay::create(3, "BDEL");
	_dt->set_keeps_infinite(true);
	_dt->set_reorder(true);
	_dt->set_vertices(_pointData.size(), reinterpret_cast<const double*>(_pointData.data()));

	spdlog::info("  Geogram Delaunay (BDEL, 1 threads) complete: {} cells", _dt->nb_cells());
}

// Marks the cells that belong to a ghost image and numbers the others in cell order.
void DelaunayTessellation::classifyCells(){
	const size_t numCells = _dt->nb_cells();
	_cellInfo.resize(numCells);

	tbb::parallel_for(tbb::blocked_range<size_t>(0, numCells, 4096),
		[&](const tbb::blocked_range<size_t>& r){
		for(size_t i = r.begin(); i < r.end(); ++i){
			_cellInfo[i].isGhost = classifyGhostCell(static_cast<CellHandle>(i));
			_cellInfo[i].index = -1;
		}
	});

	_numPrimaryTetrahedra = 0;
	for(size_t i = 0; i < numCells; ++i){
		if(!_cellInfo[i].isGhost){
			_cellInfo[i].index = _numPrimaryTetrahedra++;
		}
	}
}

bool DelaunayTessellation::classifyGhostCell(CellHandle cell) const{
	if(!isValidCell(cell)) return true;

	VertexHandle headVertex = cellVertex(cell, 0);
	int headVertexIndex = vertexIndex(headVertex);
	assert(headVertexIndex >= 0);
	for(int v = 1; v < 4; v++){
		VertexHandle p = cellVertex(cell, v);
		int vindex = vertexIndex(p);
		assert(vindex >= 0);
		if(vindex < headVertexIndex){
			headVertex = p;
			headVertexIndex = vindex;
		}
	}
	return isGhostVertex(headVertex);
}

std::optional<bool> DelaunayTessellation::alphaTest(CellHandle cell, double alpha) const{
    auto v0 = _dt->vertex_ptr(cellVertex(cell, 0));
    auto v1 = _dt->vertex_ptr(cellVertex(cell, 1));
    auto v2 = _dt->vertex_ptr(cellVertex(cell, 2));
    auto v3 = _dt->vertex_ptr(cellVertex(cell, 3));

    auto qpx = v1[0]-v0[0]; auto qpy = v1[1]-v0[1]; auto qpz = v1[2]-v0[2];
    auto qp2 = qpx*qpx + qpy*qpy + qpz*qpz;
    auto rpx = v2[0]-v0[0]; auto rpy = v2[1]-v0[1]; auto rpz = v2[2]-v0[2];
    auto rp2 = rpx*rpx + rpy*rpy + rpz*rpz;
    auto spx = v3[0]-v0[0]; auto spy = v3[1]-v0[1]; auto spz = v3[2]-v0[2];
    auto sp2 = spx*spx + spy*spy + spz*spz;

    auto num_x = determinant(qpy,qpz,qp2,rpy,rpz,rp2,spy,spz,sp2);
    auto num_y = determinant(qpx,qpz,qp2,rpx,rpz,rp2,spx,spz,sp2);
    auto num_z = determinant(qpx,qpy,qp2,rpx,rpy,rp2,spx,spy,sp2);
    auto den   = determinant(qpx,qpy,qpz,rpx,rpy,rpz,spx,spy,spz);

    double nomin = (num_x*num_x + num_y*num_y + num_z*num_z);
    double denom = (4 * den * den);

    if(std::abs(denom) < 1e-9 && std::abs(nomin) < 1e-9) return std::nullopt;
    return (nomin / denom) < alpha;
}

}
