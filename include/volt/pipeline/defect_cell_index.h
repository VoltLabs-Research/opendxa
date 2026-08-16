#pragma once

#include <volt/core/volt.h>
#include <volt/pipeline/delaunay_tessellation.h>

#include <cstddef>
#include <memory>
#include <vector>

namespace Volt{

class DefectCellIndex{
public:
    static constexpr double kAlphaScale = 3.5;

    DefectCellIndex(DelaunayTessellation& tessellation, double alpha);
    ~DefectCellIndex();

    DefectCellIndex(const DefectCellIndex&) = delete;
    DefectCellIndex& operator=(const DefectCellIndex&) = delete;

    [[nodiscard]] std::size_t cellCount() const noexcept{
        return _cellCount;
    }

    void queryOverlapping(const Box3& box, std::vector<DelaunayTessellation::CellHandle>& out) const;

private:
    struct Impl;

    std::unique_ptr<Impl> _impl;
    std::size_t _cellCount = 0;
};

}
