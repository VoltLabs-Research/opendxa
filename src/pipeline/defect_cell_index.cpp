#include <volt/pipeline/defect_cell_index.h>

#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/strategies/relate/cartesian.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/iterator/function_output_iterator.hpp>

#include <utility>

namespace Volt{

namespace{

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using IndexPoint = bg::model::point<double, 3, bg::cs::cartesian>;
using IndexBox = bg::model::box<IndexPoint>;
using IndexValue = std::pair<IndexBox, DelaunayTessellation::CellHandle>;
using CellTree = bgi::rtree<IndexValue, bgi::quadratic<128, 38>>;

IndexBox toIndexBox(const Box3& box){
    return IndexBox(
        IndexPoint(box.minc.x(), box.minc.y(), box.minc.z()),
        IndexPoint(box.maxc.x(), box.maxc.y(), box.maxc.z())
    );
}

}

struct DefectCellIndex::Impl{
    CellTree tree;

    explicit Impl(std::vector<IndexValue>&& values) : tree(std::move(values)){}
};

DefectCellIndex::DefectCellIndex(DelaunayTessellation& tessellation, double alpha){
    std::vector<IndexValue> values;

    for(auto cell : tessellation.cells()){
        const bool accepted = tessellation.isValidCell(cell) && tessellation.getUserField(cell) == 0
            && tessellation.alphaTest(cell, alpha).value_or(false);
        if(!accepted){
            tessellation.setUserField(cell, -1);
            continue;
        }

        Box3 bounds;
        for(int corner = 0; corner < 4; ++corner){
            bounds.addPoint(tessellation.vertexPosition(tessellation.cellVertex(cell, corner)));
        }

        values.emplace_back(toIndexBox(bounds), cell);
        tessellation.setUserField(cell, static_cast<int>(_cellCount++));
    }

    _impl = std::make_unique<Impl>(std::move(values));
}

DefectCellIndex::~DefectCellIndex() = default;

std::vector<DelaunayTessellation::CellHandle> DefectCellIndex::queryOverlapping(const Box3& box) const{
    std::vector<DelaunayTessellation::CellHandle> cells;
    _impl->tree.query(
        bgi::intersects(toIndexBox(box)),
        boost::make_function_output_iterator([&cells](const IndexValue& value){
            cells.push_back(value.second);
        })
    );
    return cells;
}

}
