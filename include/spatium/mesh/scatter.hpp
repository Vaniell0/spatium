#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/mesh/mesh.hpp>
#  include <spatium/mesh/topology.hpp>
#  include <spatium/mesh/voronoi.hpp>
#  include <cstdint>
#  include <limits>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::mesh {

// scatter_on_surface: place `count` points evenly across a meshed Surface
// via farthest-point sampling (repeated geodesic_voronoi() -- same
// technique examples/geodesic_procgen_demo.cpp uses for region seeding,
// here reused for point placement instead). Each returned point carries
// the surface's own analytic normal, ready to orient a small object
// dropped there (the "sprinkles" primitive: real geodesic coverage, not
// texture noise or independent RNG jitter that can cluster/gap).

template<Surface S>
struct ScatterPoint {
    typename S::PointType position;
    typename S::TangentVector normal;
};

template<Surface S>
std::vector<ScatterPoint<S>> scatter_on_surface(
    const Mesh<S>& mesh, const S& space, std::size_t count, uint32_t seed_vertex = 0)
{
    using T = typename S::ScalarType;

    auto topo = MeshTopology<S>::build(mesh);
    std::vector<uint32_t> sites{seed_vertex};
    for (std::size_t i = 1; i < count && i < mesh.vertex_count(); ++i) {
        auto vd = geodesic_voronoi(topo, space, sites);
        uint32_t farthest = 0;
        T best{-1};
        for (uint32_t v = 0; v < vd.distances.size(); ++v)
            if (vd.distances[v] > best) { best = vd.distances[v]; farthest = v; }
        sites.push_back(farthest);
    }

    std::vector<ScatterPoint<S>> out;
    out.reserve(sites.size());
    for (auto v : sites) {
        const auto& p = mesh.vertices[v];
        out.push_back({p, space.normal(p)});
    }
    return out;
}

} // namespace spatium::mesh
