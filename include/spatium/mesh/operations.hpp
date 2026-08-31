#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/mesh/mesh.hpp>
#  include <spatium/geometry/transform.hpp>
#  include <cstdint>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::mesh {

// ── Merge two meshes ──────────────────────────────────────────

template<Surface S>
Mesh<S> merge(const Mesh<S>& a, const Mesh<S>& b) {
    Mesh<S> result;
    result.vertices.reserve(a.vertices.size() + b.vertices.size());
    result.vertices.insert(result.vertices.end(), a.vertices.begin(), a.vertices.end());
    result.vertices.insert(result.vertices.end(), b.vertices.begin(), b.vertices.end());

    auto offset = static_cast<uint32_t>(a.vertices.size());
    result.faces.reserve(a.faces.size() + b.faces.size());
    result.faces.insert(result.faces.end(), a.faces.begin(), a.faces.end());
    for (auto f : b.faces) {
        f[0] += offset;
        f[1] += offset;
        f[2] += offset;
        result.faces.push_back(f);
    }

    return result;
}

// ── Transform mesh vertices ───────────────────────────────────

template<Surface S>
    requires (S::dimension <= 3)
Mesh<S> transform(const Mesh<S>& m,
                  const geometry::AffineTransform<S::dimension + 1 <= 3 ? 3 : S::dimension + 1,
                                                   typename S::ScalarType>& t) {
    Mesh<S> result;
    result.vertices.reserve(m.vertices.size());
    for (auto& v : m.vertices)
        result.vertices.push_back(t.apply(v));
    result.faces = m.faces;
    return result;
}

// Simpler overload: apply a function to each vertex
template<Surface S, typename Fn>
Mesh<S> transform(const Mesh<S>& m, Fn&& fn) {
    Mesh<S> result;
    result.vertices.reserve(m.vertices.size());
    for (auto& v : m.vertices)
        result.vertices.push_back(fn(v));
    result.faces = m.faces;
    return result;
}

// ── Flip normals (reverse face winding) ───────────────────────

template<Surface S>
Mesh<S> flip_normals(const Mesh<S>& m) {
    Mesh<S> result;
    result.vertices = m.vertices;
    result.faces.reserve(m.faces.size());
    for (auto [a, b, c] : m.faces)
        result.faces.push_back({a, c, b});
    return result;
}

// ── Compute per-face normals ──────────────────────────────────

template<Surface S>
    requires (S::dimension + 1 == 3 || std::same_as<typename S::PointType, Vec<typename S::ScalarType, 3>>)
std::vector<typename S::PointType> compute_face_normals(const Mesh<S>& m) {
    std::vector<typename S::PointType> normals;
    normals.reserve(m.faces.size());
    for (auto& [a, b, c] : m.faces)
        normals.push_back(
            face_cross_3d(m.vertices[a], m.vertices[b], m.vertices[c])
                .normalized());
    return normals;
}

// ── Compute per-vertex normals (area-weighted average) ────────

template<Surface S>
    requires (S::dimension + 1 == 3 || std::same_as<typename S::PointType, Vec<typename S::ScalarType, 3>>)
std::vector<typename S::PointType> compute_vertex_normals(const Mesh<S>& m) {
    using T = typename S::ScalarType;
    using P = typename S::PointType;
    std::vector<P> normals(m.vertices.size());

    for (auto& [a, b, c] : m.faces) {
        // NOT normalised — area-weighted on purpose, so the per-vertex
        // sum below picks up bigger contributions from larger faces.
        P face_normal = Vec<T, 3>{
            face_cross_3d(m.vertices[a], m.vertices[b], m.vertices[c])};
        normals[a] = P{normals[a] + face_normal};
        normals[b] = P{normals[b] + face_normal};
        normals[c] = P{normals[c] + face_normal};
    }

    for (auto& n : normals)
        n = n.normalized();

    return normals;
}

// ── Center mesh at origin ─────────────────────────────────────

template<Surface S>
Mesh<S> centered(const Mesh<S>& m) {
    using T = typename S::ScalarType;
    using P = typename S::PointType;

    P sum{};
    for (auto& v : m.vertices)
        sum = P{sum + v};
    P centroid = P{sum / static_cast<T>(m.vertices.size())};

    return transform<S>(m, [&](const P& v) -> P { return P{v - centroid}; });
}

} // namespace spatium::mesh
