#include <volt/core/volt.h>
#include <volt/pipeline/delaunay_tessellation.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <spdlog/spdlog.h>
#include <numeric>

namespace Volt{

void DelaunayTessellation::generateTessellation(
	const SimulationCell& simCell,
	const Point3* positions,
	size_t numPoints,
	double ghostLayerSize,
	bool coverDomainWithFiniteTets,
	const int* selectedPoints
){
	static std::mutex geogramMutex;
	{
		std::lock_guard<std::mutex> lock(geogramMutex);
		GEO::initialize(GEO::GEOGRAM_NO_HANDLER);
		GEO::set_assert_mode(GEO::ASSERT_ABORT);
	}

	double lengthScale = (simCell.matrix().column(0) + simCell.matrix().column(1) + simCell.matrix().column(2)).length();
	double epsilon = 1e-10 * lengthScale;

	_simCell = simCell;

	// Parallel point wrapping + jitter
	if(!selectedPoints){
		_pointData.resize(numPoints);
		_particleIndices.resize(numPoints);

		tbb::parallel_for(tbb::blocked_range<size_t>(0, numPoints, 4096),
			[&](const tbb::blocked_range<size_t>& r){
			std::mt19937 localRng(4 + static_cast<unsigned>(r.begin()));
			boost::random::uniform_real_distribution<double> disp(-epsilon, epsilon);
			for(size_t i = r.begin(); i < r.end(); ++i){
				Point3 wp = simCell.wrapPoint(positions[i]);
				_pointData[i] = Point3(
					(double)wp.x() + disp(localRng),
					(double)wp.y() + disp(localRng),
					(double)wp.z() + disp(localRng)
				);
				_particleIndices[i] = i;
			}
		});
		_primaryVertexCount = numPoints;
	}else{
		_particleIndices.clear();
		_pointData.clear();
		_particleIndices.reserve(numPoints);
		_pointData.reserve(numPoints);
		std::mt19937 rng(4);
		boost::random::uniform_real_distribution<double> displacement(-epsilon, epsilon);
		for(size_t i = 0; i < numPoints; i++, ++positions){
			if(!*selectedPoints++){ continue; }
			Point3 wp = simCell.wrapPoint(*positions);
			_pointData.emplace_back(
				(double)wp.x() + displacement(rng),
				(double)wp.y() + displacement(rng),
				(double)wp.z() + displacement(rng)
			);
			_particleIndices.push_back(i);
		}
		_primaryVertexCount = _particleIndices.size();
	}

	// Ghost layer (parallel count + fill)
	Vector3I stencilCount;
	double cuts[3][2];
	Vector3 cellNormals[3];
	for(size_t dim = 0; dim < 3; dim++){
		cellNormals[dim] = simCell.cellNormalVector(dim);
		cuts[dim][0] = cellNormals[dim].dot(simCell.reducedToAbsolute(Point3(0,0,0)) - Point3::Origin());
		cuts[dim][1] = cellNormals[dim].dot(simCell.reducedToAbsolute(Point3(1,1,1)) - Point3::Origin());
		if(simCell.hasPbc(dim)){
			stencilCount[dim] = (int)ceil(ghostLayerSize / simCell.matrix().column(dim).dot(cellNormals[dim]));
			cuts[dim][0] -= ghostLayerSize;
			cuts[dim][1] += ghostLayerSize;
		}else{
			stencilCount[dim] = 0;
			cuts[dim][0] -= ghostLayerSize;
			cuts[dim][1] += ghostLayerSize;
		}
	}

	std::vector<Vector3I> ghostShifts;
	ghostShifts.reserve(static_cast<size_t>(
		(2 * stencilCount[0] + 1) * (2 * stencilCount[1] + 1) * (2 * stencilCount[2] + 1) - 1));
	for(int ix = -stencilCount[0]; ix <= +stencilCount[0]; ix++)
		for(int iy = -stencilCount[1]; iy <= +stencilCount[1]; iy++)
			for(int iz = -stencilCount[2]; iz <= +stencilCount[2]; iz++){
				if(ix == 0 && iy == 0 && iz == 0) continue;
				ghostShifts.emplace_back(ix, iy, iz);
			}

	std::vector<size_t> ghostCounts(ghostShifts.size(), 0);
	tbb::parallel_for(tbb::blocked_range<size_t>(0, ghostShifts.size(), 8),
		[&](const tbb::blocked_range<size_t>& r){
		for(size_t shiftIdx = r.begin(); shiftIdx < r.end(); ++shiftIdx){
			const auto& imageShift = ghostShifts[shiftIdx];
			Vector3 shift = simCell.reducedToAbsolute(Vector3(imageShift.x(), imageShift.y(), imageShift.z()));
			size_t count = 0;
			for(size_t vertexIndex = 0; vertexIndex < _primaryVertexCount; ++vertexIndex){
				Point3 pimage = _pointData[vertexIndex] + shift;
				bool isClipped = false;
				for(size_t dim = 0; dim < 3; ++dim){
					if(simCell.hasPbc(dim)){
						double d = cellNormals[dim].dot(pimage - Point3::Origin());
						if(d < cuts[dim][0] || d > cuts[dim][1]){ isClipped = true; break; }
					}
				}
				if(!isClipped) ++count;
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
			const auto& imageShift = ghostShifts[shiftIdx];
			Vector3 shift = simCell.reducedToAbsolute(Vector3(imageShift.x(), imageShift.y(), imageShift.z()));
			size_t writeIndex = basePointCount + ghostOffsets[shiftIdx];
			for(size_t vertexIndex = 0; vertexIndex < _primaryVertexCount; ++vertexIndex){
				Point3 pimage = _pointData[vertexIndex] + shift;
				bool isClipped = false;
				for(size_t dim = 0; dim < 3; ++dim){
					if(simCell.hasPbc(dim)){
						double d = cellNormals[dim].dot(pimage - Point3::Origin());
						if(d < cuts[dim][0] || d > cuts[dim][1]){ isClipped = true; break; }
					}
				}
				if(!isClipped){
					_pointData[writeIndex] = pimage;
					_particleIndices[writeIndex] = _particleIndices[vertexIndex];
					++writeIndex;
				}
			}
		}
	});

	std::vector<Vector3I>().swap(ghostShifts);
	std::vector<size_t>().swap(ghostCounts);
	std::vector<size_t>().swap(ghostOffsets);

	if(coverDomainWithFiniteTets){
		Box3 bb = Box3(Point3(0, 0, 0), Point3(1, 1, 1)).transformed(simCell.matrix());
		bb.addPoints(_pointData.data(), _pointData.size());
		bb = bb.padBox(ghostLayerSize);
		for(size_t i = 0; i < 8; i++){
			Point3 corner = bb[static_cast<int>(i)];
			_pointData.push_back(corner);
			_particleIndices.push_back(std::numeric_limits<size_t>::max());
		}
	}

	spdlog::info("  Geogram Delaunay (serial BDEL): inserting {} points ({} primary + {} ghost)",
		_pointData.size(), _primaryVertexCount, _pointData.size() - _primaryVertexCount);

	// Use BDEL (serial Delaunay3d) instead of PDEL (ParallelDelaunay3d).
	// PDEL deadlocks/livelocks on high-core hosts: reproduced on a 160-core
	// EPYC where it hangs at point insertion regardless of the configured
	// thread count (160 threads -> spin at 99.9% CPU; 1 thread -> sleeping
	// deadlock). The serial backend is robust and fast enough for our point
	// counts (tens of thousands of points).
	_dt = GEO::Delaunay::create(3, "BDEL");
	_dt->set_keeps_infinite(true);
	_dt->set_reorder(true);
	_dt->set_vertices(_pointData.size(), reinterpret_cast<const double*>(_pointData.data()));

	spdlog::info("  Geogram Delaunay (serial BDEL) complete: {} cells", _dt->nb_cells());

	// Classify cells (parallel)
	const size_t numCells = _dt->nb_cells();
	_cellInfo.resize(numCells);

	tbb::parallel_for(tbb::blocked_range<size_t>(0, numCells, 4096),
		[&](const tbb::blocked_range<size_t>& r){
		for(size_t i = r.begin(); i < r.end(); ++i){
			CellHandle cell = static_cast<CellHandle>(i);
			_cellInfo[i].isGhost = classifyGhostCell(cell);
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

static inline double determinant(double a00, double a01, double a02,
                                 double a10, double a11, double a12,
                                 double a20, double a21, double a22){
    double m02 = a00*a21 - a20*a01;
    double m01 = a00*a11 - a10*a01;
    double m12 = a10*a21 - a20*a11;
    double m012 = m01*a22 - m02*a12 + m12*a02;
    return m012;
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

void DelaunayTessellation::releaseMemory() noexcept{
	_dt.reset();
	std::vector<Point3>().swap(_pointData);
	std::vector<CellInfo>().swap(_cellInfo);
	std::vector<size_t>().swap(_particleIndices);
	_primaryVertexCount = 0;
	_numPrimaryTetrahedra = 0;
}

}
