#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/mesh/mesh.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <algorithm>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium::mesh {

template<Surface S>
struct MeshQuality {
    using T = typename S::ScalarType;
    T min_angle{};
    T max_angle{};
    T avg_aspect_ratio{};
    std::size_t degenerate_count{};
};

// Aspect ratio of a single face (longest edge / shortest edge)
template<Surface S>
typename S::ScalarType face_aspect_ratio(const Mesh<S>& m, std::size_t fi) {
    using T = typename S::ScalarType;
    auto [a, b, c] = m.faces[fi];
    T e0 = (m.vertices[b] - m.vertices[a]).norm();
    T e1 = (m.vertices[c] - m.vertices[b]).norm();
    T e2 = (m.vertices[a] - m.vertices[c]).norm();
    T lo = std::min({e0, e1, e2});
    T hi = std::max({e0, e1, e2});
    return lo > epsilon<T>() ? hi / lo : std::numeric_limits<T>::max();
}

// Minimum angle of a triangle (radians)
template<Surface S>
typename S::ScalarType face_min_angle(const Mesh<S>& m, std::size_t fi) {
    using T = typename S::ScalarType;
    auto [a, b, c] = m.faces[fi];
    auto ab = m.vertices[b] - m.vertices[a];
    auto ac = m.vertices[c] - m.vertices[a];
    auto bc = m.vertices[c] - m.vertices[b];

    auto angle = [](const auto& u, const auto& v) -> T {
        using std::acos;
        T d = static_cast<T>(u.dot(v)) / (u.norm() * v.norm());
        return acos(std::clamp(d, T{-1}, T{1}));
    };

    using P = typename S::PointType;
    T a0 = angle(ab, ac);
    T a1 = angle(P{ab * T{-1}}, bc);
    T a2 = angle(P{ac * T{-1}}, P{bc * T{-1}});
    return std::min({a0, a1, a2});
}

template<Surface S>
MeshQuality<S> mesh_quality(const Mesh<S>& m) {
    using T = typename S::ScalarType;
    MeshQuality<S> q;
    q.min_angle = std::numbers::pi_v<T>;
    q.max_angle = T{0};

    T aspect_sum{0};
    for (std::size_t fi = 0; fi < m.faces.size(); ++fi) {
        T ar = face_aspect_ratio<S>(m, fi);
        aspect_sum += ar;

        T min_a = face_min_angle<S>(m, fi);
        q.min_angle = std::min(q.min_angle, min_a);
        q.max_angle = std::max(q.max_angle, std::numbers::pi_v<T> - T{2} * min_a);

        // Degenerate: near-zero area
        auto [a, b, c] = m.faces[fi];
        if constexpr (S::PointType::size == 3) {
            if (face_area_3d(m.vertices[a], m.vertices[b], m.vertices[c])
                    < T{0.5} * epsilon<T>())
                q.degenerate_count++;
        }
    }

    q.avg_aspect_ratio = m.faces.empty() ? T{0} : aspect_sum / static_cast<T>(m.faces.size());
    return q;
}

} // namespace spatium::mesh
