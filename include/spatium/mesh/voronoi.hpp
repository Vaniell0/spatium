#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/mesh/topology.hpp>
#  include <cstdint>
#  include <functional>
#  include <limits>
#  include <queue>
#  include <span>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::mesh {

template<Surface S>
struct VoronoiDiagram {
    using ScalarT = typename S::ScalarType;

    std::vector<uint32_t> labels;      // per-vertex: index into sites
    std::vector<ScalarT> distances;    // per-vertex: distance to nearest site
    std::vector<uint32_t> sites;
};

template<Surface S>
VoronoiDiagram<S> geodesic_voronoi(
    const MeshTopology<S>& topo,
    const S& space,
    std::span<const uint32_t> sites)
{
    using T = typename S::ScalarType;
    auto nv = topo.vertex_count();
    auto& mesh = topo.mesh();

    VoronoiDiagram<S> vd;
    vd.sites.assign(sites.begin(), sites.end());
    vd.distances.assign(nv, std::numeric_limits<T>::max());
    vd.labels.assign(nv, no_vertex);

    // Min-heap: (distance, vertex)
    using Entry = std::pair<T, uint32_t>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<>> pq;

    for (uint32_t si = 0; si < sites.size(); ++si) {
        auto v = sites[si];
        vd.distances[v] = T{0};
        vd.labels[v] = si;
        pq.push({T{0}, v});
    }

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();

        if (d > vd.distances[u]) continue;

        for (auto v : topo.neighbors(u)) {
            T w = space.distance(mesh.vertices[u], mesh.vertices[v]);
            T new_dist = d + w;
            if (new_dist < vd.distances[v]) {
                vd.distances[v] = new_dist;
                vd.labels[v] = vd.labels[u];
                pq.push({new_dist, v});
            }
        }
    }

    return vd;
}

// Per-face label: majority vote of vertex labels.
// Returns no_vertex for faces with mixed labels (boundary faces).
template<Surface S>
std::vector<uint32_t> face_labels(
    const VoronoiDiagram<S>& vd,
    const Mesh<S>& mesh)
{
    std::vector<uint32_t> result(mesh.faces.size());
    for (std::size_t f = 0; f < mesh.faces.size(); ++f) {
        auto [a, b, c] = mesh.faces[f];
        auto la = vd.labels[a], lb = vd.labels[b], lc = vd.labels[c];
        if (la == lb || la == lc)
            result[f] = la;
        else if (lb == lc)
            result[f] = lb;
        else
            result[f] = no_vertex;  // boundary face
    }
    return result;
}

} // namespace spatium::mesh
