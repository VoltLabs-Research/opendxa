#pragma once

#include <volt/core/volt.h>
#include <volt/pipeline/delaunay_tessellation.h>

#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/strategies/relate/cartesian.hpp>

#include <cstddef>
#include <utility>
#include <vector>

namespace Volt{

class DefectCellIndex{
public:
    static constexpr double kAlphaScale = 3.5;

    DefectCellIndex(DelaunayTessellation& tessellation, double alpha);

    DefectCellIndex(const DefectCellIndex&) = delete;
    DefectCellIndex& operator=(const DefectCellIndex&) = delete;

    [[nodiscard]] std::size_t cellCount() const noexcept{
        return _cellCount;
    }

    [[nodiscard]] std::vector<DelaunayTessellation::CellHandle> queryOverlapping(const Box3& box) const;

private:
    using IndexPoint = boost::geometry::model::point<double, 3, boost::geometry::cs::cartesian>;
    using IndexBox = boost::geometry::model::box<IndexPoint>;
    using IndexValue = std::pair<IndexBox, DelaunayTessellation::CellHandle>;
    using CellTree = boost::geometry::index::rtree<IndexValue, boost::geometry::index::quadratic<128, 38>>;

    static IndexBox toIndexBox(const Box3& box);

    CellTree _tree;
    std::size_t _cellCount = 0;
};

}
