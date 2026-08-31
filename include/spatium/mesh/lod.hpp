#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/mesh/mesh.hpp>
#  include <spatium/mesh/subdivision.hpp>
#  include <stdexcept>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::mesh {

// LOD chain: stores meshes at increasing subdivision levels.
// Level 0 = base mesh, level K = K subdivisions.

template<Surface S>
struct LodChain {
    std::vector<Mesh<S>> levels;

    static LodChain build(const Mesh<S>& base, const S& space, std::size_t max_level) {
        LodChain chain;
        chain.levels.reserve(max_level + 1);
        chain.levels.push_back(base);
        for (std::size_t i = 0; i < max_level; ++i)
            chain.levels.push_back(subdivide_once(chain.levels.back(), space));
        return chain;
    }

    [[nodiscard]] std::size_t level_count() const { return levels.size(); }

    [[nodiscard]] const Mesh<S>& at(std::size_t level) const {
        if (levels.empty()) throw std::out_of_range("LodChain is empty");
        return levels[std::min(level, levels.size() - 1)];
    }

    [[nodiscard]] const Mesh<S>& coarsest() const { return levels.front(); }
    [[nodiscard]] const Mesh<S>& finest() const { return levels.back(); }
};

} // namespace spatium::mesh
