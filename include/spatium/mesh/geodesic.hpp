#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/mesh/geodesic_types.hpp>
#  include <spatium/mesh/topology.hpp>
#  include <cstdint>
#  include <functional>
#  include <limits>
#  include <queue>
#  include <span>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::mesh {

// ── Single-source Dijkstra ────────────────────────────────────

template<Surface S>
DistanceField<S> geodesic_distances(
    const MeshTopology<S>& topo,
    const S& space,
    uint32_t source)
{
    return geodesic_distances(topo, space, std::span<const uint32_t>(&source, 1));
}

// ── Multi-source Dijkstra ─────────────────────────────────────

template<Surface S>
DistanceField<S> geodesic_distances(
    const MeshTopology<S>& topo,
    const S& space,
    std::span<const uint32_t> sources)
{
    using T = typename S::ScalarType;
    auto nv = topo.vertex_count();
    auto& mesh = topo.mesh();

    DistanceField<S> field;
    field.distances.assign(nv, std::numeric_limits<T>::max());
    field.predecessors.assign(nv, no_vertex);

    // Min-heap: (distance, vertex)
    using Entry = std::pair<T, uint32_t>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<>> pq;

    for (auto s : sources) {
        field.distances[s] = T{0};
        pq.push({T{0}, s});
    }

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();

        if (d > field.distances[u]) continue;  // stale entry

        for (auto v : topo.neighbors(u)) {
            T w = space.distance(mesh.vertices[u], mesh.vertices[v]);
            T new_dist = d + w;
            if (new_dist < field.distances[v]) {
                field.distances[v] = new_dist;
                field.predecessors[v] = u;
                pq.push({new_dist, v});
            }
        }
    }

    return field;
}

// ── Shortest path ─────────────────────────────────────────────

template<Surface S>
GeodesicPath<S> shortest_path(
    const MeshTopology<S>& topo,
    const S& space,
    uint32_t source,
    uint32_t target)
{
    auto field = geodesic_distances(topo, space, source);

    GeodesicPath<S> path;
    path.total_length = field.distances[target];

    if (field.distances[target] == std::numeric_limits<typename S::ScalarType>::max())
        return path;  // unreachable

    // Trace predecessors
    uint32_t cur = target;
    while (cur != no_vertex) {
        path.vertices.push_back(cur);
        if (cur == source) break;
        cur = field.predecessors[cur];
    }

    std::reverse(path.vertices.begin(), path.vertices.end());
    return path;
}

// ── Method selection ──────────────────────────────────────────

enum class GeodesicMethod { Dijkstra, Heat };

} // namespace spatium::mesh

// Include heat method after Dijkstra definitions (heat_geodesic.hpp
// uses geodesic_types.hpp for DistanceField, no circular dependency).
// In module mode heat_geodesic stays header-only (see mesh.cppm for
// rationale — Eigen/Sparse + vec_simd BMI conflict); consumers include
// `<spatium/mesh/heat_geodesic.hpp>` directly when they need it.
#if defined(SPATIUM_HAS_EIGEN) && SPATIUM_HAS_EIGEN && !defined(SPATIUM_BUILDING_MODULE)
#include <spatium/mesh/heat_geodesic.hpp>
#endif

namespace spatium::mesh {

template<Surface S>
DistanceField<S> geodesic_distances(
    const MeshTopology<S>& topo,
    const S& space,
    uint32_t source,
    GeodesicMethod method)
{
    if (method == GeodesicMethod::Dijkstra)
        return geodesic_distances(topo, space, source);
#if defined(SPATIUM_HAS_EIGEN) && SPATIUM_HAS_EIGEN
    return heat_geodesic_distances(topo, space, source);
#else
    static_assert(false, "Heat method requires SPATIUM_HAS_EIGEN");
#endif
}

} // namespace spatium::mesh
