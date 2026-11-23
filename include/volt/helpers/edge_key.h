#pragma once

#include <cstdint>
#include <utility>

namespace Volt{

// One 64-bit key per undirected edge, with the smaller vertex index in the high half: for non-negative indices,
// sorting the keys sorts the edges by their first vertex, then by their second.
inline std::uint64_t packEdgeKey(int a, int b) noexcept{
    if(a > b) std::swap(a, b);
    return (std::uint64_t{static_cast<std::uint32_t>(a)} << 32) | std::uint64_t{static_cast<std::uint32_t>(b)};
}

inline std::pair<int, int> unpackEdgeKey(std::uint64_t key) noexcept{
    return {static_cast<int>(key >> 32), static_cast<int>(key & 0xFFFFFFFFu)};
}

}
