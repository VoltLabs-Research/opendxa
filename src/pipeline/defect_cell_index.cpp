#include <volt/pipeline/defect_cell_index.h>

#include <boost/iterator/function_output_iterator.hpp>

#include <utility>

namespace Volt{

DefectCellIndex::IndexBox DefectCellIndex::toIndexBox(const Box3& box){
    return IndexBox(
        IndexPoint(box.minc.x(), box.minc.y(), box.minc.z()),
        IndexPoint(box.maxc.x(), box.maxc.y(), box.maxc.z())
    );
}

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

    _tree = CellTree(std::move(values));
}

std::vector<DelaunayTessellation::CellHandle> DefectCellIndex::queryOverlapping(const Box3& box) const{
    std::vector<DelaunayTessellation::CellHandle> cells;
    _tree.query(
        boost::geometry::index::intersects(toIndexBox(box)),
        boost::make_function_output_iterator([&cells](const IndexValue& value){
            cells.push_back(value.second);
        })
    );
    return cells;
}

}
