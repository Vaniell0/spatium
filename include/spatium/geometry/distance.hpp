#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/geometry/line.hpp>
#  include <spatium/geometry/hyperplane.hpp>
#  include <spatium/geometry/triangle.hpp>
#  include <spatium/geometry/box.hpp>
#  include <spatium/geometry/circle.hpp>
#  include <spatium/geometry/polygon.hpp>
#  include <cmath>
#  include <algorithm>
#endif

SPATIUM_EXPORT namespace spatium::geometry {

// ── Point-to-shape distances ───────────────────────────────────
// (Most shapes already have .distance(point) members.
//  These free functions provide uniform syntax.)

template<std::size_t N, Scalar T>
T distance(const Vec<T, N>& p, const Line<N, T>& line) {
    return line.distance(p);
}

template<std::size_t N, Scalar T>
T distance(const Vec<T, N>& p, const Ray<N, T>& ray) {
    return ray.distance(p);
}

template<std::size_t N, Scalar T>
T distance(const Vec<T, N>& p, const Segment<N, T>& seg) {
    return seg.distance(p);
}

template<std::size_t N, Scalar T>
T distance(const Vec<T, N>& p, const Hyperplane<N, T>& plane) {
    return plane.distance(p);
}

template<std::size_t N, Scalar T>
T distance(const Vec<T, N>& p, const Triangle<N, T>& tri) {
    return tri.distance(p);
}

template<std::size_t N, Scalar T>
T distance(const Vec<T, N>& p, const Box<N, T>& box) {
    return box.distance(p);
}

// ── Line-Line distance (3D) ───────────────────────────────────

template<Scalar T>
T distance(const Line<3, T>& a, const Line<3, T>& b) {
    auto w = a.origin - b.origin;
    auto u = a.direction;
    auto v = b.direction;
    auto a_val = u.dot(u);
    auto b_val = u.dot(v);
    auto c_val = v.dot(v);
    auto d_val = u.dot(w);
    auto e_val = v.dot(w);
    auto denom = a_val * c_val - b_val * b_val;

    if (std::abs(denom) < epsilon<T>()) {
        // Parallel lines: distance = distance from any point on a to line b
        return b.distance(a.origin);
    }

    auto sc = (b_val * e_val - c_val * d_val) / denom;
    auto tc = (a_val * e_val - b_val * d_val) / denom;
    auto dp = w + u * sc - v * tc;
    return dp.norm();
}

// ── Segment-Segment distance ──────────────────────────────────

template<std::size_t N, Scalar T>
T distance(const Segment<N, T>& s1, const Segment<N, T>& s2) {
    // Brute force for correctness: check all combinations
    T min_d = std::numeric_limits<T>::max();

    // Each endpoint of s1 to s2, and vice versa
    min_d = std::min(min_d, s2.distance(s1.a));
    min_d = std::min(min_d, s2.distance(s1.b));
    min_d = std::min(min_d, s1.distance(s2.a));
    min_d = std::min(min_d, s1.distance(s2.b));

    // For 3D: also check closest approach of interior
    if constexpr (N == 3) {
        auto d1 = s1.b - s1.a;
        auto d2 = s2.b - s2.a;
        auto r = s1.a - s2.a;
        auto a_val = d1.dot(d1);
        auto e_val = d2.dot(d2);
        auto f_val = d2.dot(r);
        auto b_val = d1.dot(d2);
        auto c_val = d1.dot(r);
        auto denom = a_val * e_val - b_val * b_val;
        if (std::abs(denom) > epsilon<T>()) {
            auto s = std::clamp((b_val * f_val - c_val * e_val) / denom, T{0}, T{1});
            auto t = std::clamp((a_val * f_val - b_val * c_val) / denom, T{0}, T{1});
            auto closest = r + d1 * s - d2 * t;
            min_d = std::min(min_d, closest.norm());
        }
    }

    return min_d;
}

// ── Point-to-shape (additional wrappers) ──────────────────────

template<std::size_t N, Scalar T>
T distance(const Vec<T, N>& p, const Circle<N, T>& c) {
    return c.distance(p);
}

template<std::size_t N, Scalar T>
T distance(const Vec<T, N>& p, const Disk<N, T>& d) {
    return d.distance(p);
}

template<std::size_t N, Scalar T>
T distance(const Vec<T, N>& p, const Polygon<N, T>& poly) {
    return poly.distance(p);
}

// ── Shape-to-shape distances ──────────────────────────────────

// Box-Box: per-axis gap distance
template<std::size_t N, Scalar T>
T distance(const Box<N, T>& a, const Box<N, T>& b) {
    T sum_sq{0};
    for (std::size_t i = 0; i < N; ++i) {
        T gap = std::max({T{0}, a.min_corner[i] - b.max_corner[i],
                                b.min_corner[i] - a.max_corner[i]});
        sum_sq += gap * gap;
    }
    return std::sqrt(sum_sq);
}

// Circle-Circle (2D): gap between two circles
template<Scalar T>
T distance(const Circle<2, T>& a, const Circle<2, T>& b) {
    auto center_dist = (a.center - b.center).norm();
    return std::max(T{0}, center_dist - a.radius - b.radius);
}

// Triangle-Triangle: min edge-pair distance + containment
template<std::size_t N, Scalar T>
T distance(const Triangle<N, T>& a, const Triangle<N, T>& b) {
    // Check containment (any vertex of one inside the other)
    for (std::size_t i = 0; i < 3; ++i) {
        if (a.contains(b[i]) || b.contains(a[i])) return T{0};
    }
    // Min distance across all edge pairs (3x3 = 9)
    T min_d = std::numeric_limits<T>::max();
    for (std::size_t i = 0; i < 3; ++i) {
        auto ea = Segment<N, T>{a[i], a[(i + 1) % 3]};
        for (std::size_t j = 0; j < 3; ++j) {
            auto eb = Segment<N, T>{b[j], b[(j + 1) % 3]};
            min_d = std::min(min_d, distance(ea, eb));
        }
    }
    return min_d;
}

// Polygon-Polygon: min edge-pair distance + containment
template<std::size_t N, Scalar T>
T distance(const Polygon<N, T>& a, const Polygon<N, T>& b) {
    // Containment check
    if (!a.vertices.empty() && b.contains(a.vertices[0])) return T{0};
    if (!b.vertices.empty() && a.contains(b.vertices[0])) return T{0};
    // Min distance across all edge pairs
    T min_d = std::numeric_limits<T>::max();
    for (std::size_t i = 0; i < a.vertices.size(); ++i) {
        auto ea = a.edge(i);
        for (std::size_t j = 0; j < b.vertices.size(); ++j) {
            auto eb = b.edge(j);
            min_d = std::min(min_d, distance(ea, eb));
        }
    }
    return min_d;
}

// Segment-Triangle: min of seg vs 3 edges + point-tri checks
template<std::size_t N, Scalar T>
T distance(const Segment<N, T>& seg, const Triangle<N, T>& tri) {
    // Check if segment endpoints are inside triangle
    if (tri.contains(seg.a)) return T{0};
    if (tri.contains(seg.b)) return T{0};
    // Check segment against each triangle edge
    T min_d = std::numeric_limits<T>::max();
    for (std::size_t i = 0; i < 3; ++i) {
        auto te = Segment<N, T>{tri[i], tri[(i + 1) % 3]};
        min_d = std::min(min_d, distance(seg, te));
    }
    // Also check triangle vertices against segment
    for (std::size_t i = 0; i < 3; ++i)
        min_d = std::min(min_d, seg.distance(tri[i]));
    return min_d;
}

} // namespace spatium::geometry
