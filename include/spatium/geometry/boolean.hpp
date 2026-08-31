#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/core/error.hpp>
#  include <spatium/geometry/concepts.hpp>
#  include <spatium/geometry/polygon.hpp>
#  include <spatium/geometry/triangle.hpp>
#  include <spatium/geometry/circle.hpp>
#  include <spatium/geometry/convex_hull.hpp>
#  include <algorithm>
#  include <cmath>
#  include <numbers>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::geometry {

// ── Helpers ───────────────────────────────────────────────────

// internal — do not use, no API stability. Houses the
// to_vertices()/Sutherland-Hodgman clipping helpers shared
// across the polygon boolean operators (intersection_region,
// difference_region, etc.). Public users should call those free
// functions; the helpers may be moved or rewritten any time.
namespace detail {

// Convert a ClosedShape to a vertex list (polygon approximation).
template<std::size_t N, Scalar T>
std::vector<Vec<T, N>> to_vertices(const Triangle<N, T>& tri) {
    return {tri[0], tri[1], tri[2]};
}

template<std::size_t N, Scalar T>
std::vector<Vec<T, N>> to_vertices(const Polygon<N, T>& poly) {
    return poly.vertices;
}

template<std::size_t N, Scalar T>
std::vector<Vec<T, N>> to_vertices(const Box<N, T>& box) requires (N == 2) {
    auto& lo = box.min_corner;
    auto& hi = box.max_corner;
    return {{lo[0], lo[1]}, {hi[0], lo[1]}, {hi[0], hi[1]}, {lo[0], hi[1]}};
}

// Approximate a disk as a regular polygon (for boolean ops).
template<Scalar T>
std::vector<Vec<T, 3>> to_vertices(const Disk<3, T>& disk, std::size_t segments = 32) {
    auto& c = disk.boundary.center;
    auto r = disk.boundary.radius;
    auto& n = disk.boundary.normal;

    // Build orthonormal basis in the disk's plane
    Vec<T, 3> u{};
    for (std::size_t i = 0; i < 3; ++i) {
        if (std::abs(n[i]) < T{0.9}) { u[i] = T{1}; break; }
    }
    u = u - n * n.dot(u);
    u = u / u.norm();
    auto v = n.cross(u);

    std::vector<Vec<T, 3>> verts;
    verts.reserve(segments);
    for (std::size_t i = 0; i < segments; ++i) {
        auto angle = T{2} * std::numbers::pi_v<T> * T(i) / T(segments);
        verts.push_back(c + u * (r * std::cos(angle)) + v * (r * std::sin(angle)));
    }
    return verts;
}

// 2D Sutherland-Hodgman: clip subject polygon by a convex clip polygon.
// All vertices assumed coplanar. Works in 2D or projected 3D.
template<std::size_t N, Scalar T>
std::vector<Vec<T, N>> sutherland_hodgman(
    const std::vector<Vec<T, N>>& subject,
    const std::vector<Vec<T, N>>& clip_poly,
    const Vec<T, N>& plane_normal = {})  // needed for 3D cross product orientation
{
    if (subject.empty() || clip_poly.empty()) return {};

    auto output = subject;
    auto cn = clip_poly.size();

    for (std::size_t i = 0; i < cn; ++i) {
        if (output.empty()) return {};
        auto input = std::move(output);
        output.clear();

        auto& edge_a = clip_poly[i];
        auto& edge_b = clip_poly[(i + 1) % cn];
        auto edge = edge_b - edge_a;

        // Determine inside/outside for each point
        auto is_inside = [&](const Vec<T, N>& p) -> bool {
            auto d = p - edge_a;
            if constexpr (N == 2) {
                return edge[0] * d[1] - edge[1] * d[0] >= -epsilon<T>();
            } else {
                // 3D: use cross product with plane normal to determine side
                auto cross = edge.cross(d);
                return cross.dot(plane_normal) >= -epsilon<T>();
            }
        };

        auto line_intersect = [&](const Vec<T, N>& p1, const Vec<T, N>& p2) -> Vec<T, N> {
            auto d1 = p1 - edge_a;
            auto d2 = p2 - edge_a;
            T c1, c2;
            if constexpr (N == 2) {
                c1 = edge[0] * d1[1] - edge[1] * d1[0];
                c2 = edge[0] * d2[1] - edge[1] * d2[0];
            } else {
                c1 = edge.cross(d1).dot(plane_normal);
                c2 = edge.cross(d2).dot(plane_normal);
            }
            auto t = c1 / (c1 - c2);
            return p1 + (p2 - p1) * t;
        };

        auto in_count = input.size();
        for (std::size_t j = 0; j < in_count; ++j) {
            auto& current = input[j];
            auto& prev = input[(j + in_count - 1) % in_count];
            bool curr_in = is_inside(current);
            bool prev_in = is_inside(prev);

            if (curr_in) {
                if (!prev_in)
                    output.push_back(line_intersect(prev, current));
                output.push_back(current);
            } else if (prev_in) {
                output.push_back(line_intersect(prev, current));
            }
        }
    }

    return output;
}

} // namespace detail

// ── Boolean operations for coplanar shapes ────────────────────

// Intersection region of two convex coplanar shapes → Polygon.
// Shapes must be convertible to vertex lists (Triangle, Polygon, Disk, Box2D).

template<std::size_t N, Scalar T, typename A, typename B>
    requires Shape<A> && Shape<B>
Result<Polygon<N, T>> intersection_region(const A& a, const B& b) {
    auto va = detail::to_vertices(a);
    auto vb = detail::to_vertices(b);

    if (va.size() < 3 || vb.size() < 3)
        return std::unexpected(Error{ErrorCode::DegenerateInput, "need 3+ vertices"});

    Vec<T, N> normal{};
    if constexpr (N >= 3) {
        // Compute plane normal from first shape
        auto e1 = va[1] - va[0];
        auto e2 = va[2] - va[0];
        normal = e1.cross(e2).normalized();
    }

    auto result = detail::sutherland_hodgman(va, vb, normal);

    if (result.size() < 3)
        return std::unexpected(Error{ErrorCode::NoIntersection, "no overlap region"});

    return Polygon<N, T>{std::move(result)};
}

// Convenience: Triangle ∩ Triangle → Polygon
template<Scalar T>
Result<Polygon<3, T>> intersection_region(const Triangle<3, T>& a, const Triangle<3, T>& b) {
    return intersection_region<3, T, Triangle<3, T>, Triangle<3, T>>(a, b);
}

template<Scalar T>
Result<Polygon<2, T>> intersection_region(const Triangle<2, T>& a, const Triangle<2, T>& b) {
    return intersection_region<2, T, Triangle<2, T>, Triangle<2, T>>(a, b);
}

// Convenience: Polygon ∩ Polygon
template<std::size_t N, Scalar T>
Result<Polygon<N, T>> intersection_region(const Polygon<N, T>& a, const Polygon<N, T>& b) {
    return intersection_region<N, T, Polygon<N, T>, Polygon<N, T>>(a, b);
}

// Convenience: Triangle ∩ Disk → Polygon (3D coplanar)
template<Scalar T>
Result<Polygon<3, T>> intersection_region(const Triangle<3, T>& a, const Disk<3, T>& b) {
    return intersection_region<3, T, Triangle<3, T>, Disk<3, T>>(a, b);
}

// Difference: points in A but not in B → Polygon.
// For convex shapes: subtract B from A using Sutherland-Hodgman against inverted B edges.
// Simplified: return the A vertices that aren't in B, plus intersection boundary points.
// Full implementation is complex for general polygons. Here we provide area-based check.
template<std::size_t N, Scalar T, typename A, typename B>
    requires Measurable<A> && Measurable<B>
Result<T> difference_area(const A& a, const B& b) {
    auto inter = intersection_region<N, T>(a, b);
    if (!inter) {
        // No intersection — difference is all of A
        Polygon<N, T> pa{detail::to_vertices(a)};
        return pa.area();
    }
    Polygon<N, T> pa{detail::to_vertices(a)};
    return pa.area() - inter->area();
}

// Symmetric difference area: area(A) + area(B) - 2 * area(A ∩ B)
template<std::size_t N, Scalar T, typename A, typename B>
    requires Measurable<A> && Measurable<B>
Result<T> symmetric_difference_area(const A& a, const B& b) {
    Polygon<N, T> pa{detail::to_vertices(a)};
    Polygon<N, T> pb{detail::to_vertices(b)};
    auto inter = intersection_region<N, T>(a, b);
    T inter_area = inter ? inter->area() : T{0};
    return pa.area() + pb.area() - T{2} * inter_area;
}

// ── Difference region (convex polygons) ───────────────────────
// A \ B for convex polygons. Walks the boundary of A, collecting
// vertices outside B and edge-edge intersection points.
// Result is a single convex polygon (or error if A ⊂ B).

template<std::size_t N, Scalar T>
Result<Polygon<N, T>> difference_region(const Polygon<N, T>& a, const Polygon<N, T>& b) {
    if (a.vertices.size() < 3 || b.vertices.size() < 3)
        return std::unexpected(Error{ErrorCode::DegenerateInput, "need 3+ vertices"});

    // Compute intersection first
    auto inter = intersection_region(a, b);

    if (!inter) {
        // No intersection — difference is all of A
        return a;
    }

    // Check if A is entirely inside B (intersection ≈ A)
    Polygon<N, T> pa{a.vertices};
    if (std::abs(inter->area() - pa.area()) < epsilon<T>() * pa.area() + epsilon<T>())
        return std::unexpected(Error{ErrorCode::NoIntersection, "A is entirely inside B"});

    // Build difference polygon: walk edges of A, keep parts outside B,
    // and add intersection boundary points where A exits/enters B.
    // For convex A \ convex B (partial overlap), the result is the
    // vertices of A that are outside B, plus the intersection boundary
    // points on the border of B, ordered correctly.

    std::vector<Vec<T, N>> result;

    auto an = a.vertices.size();
    auto bn = b.vertices.size();

    // For each edge of A, find intersection points with edges of B
    for (std::size_t i = 0; i < an; ++i) {
        auto& p1 = a.vertices[i];
        auto& p2 = a.vertices[(i + 1) % an];
        bool p1_in_b = b.contains(p1);

        if (!p1_in_b) result.push_back(p1);

        // Find intersections of edge (p1,p2) with all edges of B
        std::vector<std::pair<T, Vec<T, N>>> hits;
        auto edge_a = p2 - p1;
        for (std::size_t j = 0; j < bn; ++j) {
            auto& q1 = b.vertices[j];
            auto& q2 = b.vertices[(j + 1) % bn];
            auto edge_b = q2 - q1;
            auto denom = edge_a[0] * edge_b[1] - edge_a[1] * edge_b[0];
            if (std::abs(denom) < epsilon<T>()) continue;
            auto w = p1 - q1;
            auto ta = (edge_b[0] * w[1] - edge_b[1] * w[0]) / (-denom);
            auto tb = (edge_a[0] * w[1] - edge_a[1] * w[0]) / (-denom);
            if (ta > epsilon<T>() && ta < T{1} - epsilon<T>() &&
                tb > -epsilon<T>() && tb < T{1} + epsilon<T>()) {
                hits.emplace_back(ta, p1 + edge_a * ta);
            }
        }
        std::sort(hits.begin(), hits.end(), [](auto& a, auto& b) { return a.first < b.first; });
        for (auto& [t, pt] : hits) result.push_back(pt);
    }

    if (result.size() < 3)
        return std::unexpected(Error{ErrorCode::NoIntersection, "difference is degenerate"});

    return Polygon<N, T>{std::move(result)};
}

// ── Polygon boolean operators ─────────────────────────────────
// For convex polygons. operator^ omitted (symmetric difference
// can produce disjoint regions; use symmetric_difference_area()).

// Intersection: A ∩ B
template<std::size_t N, Scalar T>
Result<Polygon<N, T>> operator&(const Polygon<N, T>& a, const Polygon<N, T>& b) {
    return intersection_region(a, b);
}

// Difference: A \ B (2D only — `difference_region` evaluates a planar
// cross product internally, mirroring `operator+` which is also 2D).
template<Scalar T>
Result<Polygon<2, T>> operator-(const Polygon<2, T>& a, const Polygon<2, T>& b) {
    return difference_region(a, b);
}

// Union: convex hull of combined vertices (2D only)
// Uses operator+ (not operator|) to avoid ambiguity with intersection pipe.
// Convention: & = intersection, - = difference, + = union, | = intersect/pipe.
// Result<...> matches the signatures of operator& and operator- and
// surfaces the DegenerateInput from convex_hull when the combined
// vertex count is below three.
template<Scalar T>
Result<Polygon<2, T>> operator+(const Polygon<2, T>& a, const Polygon<2, T>& b) {
    auto combined = a.vertices;
    combined.insert(combined.end(), b.vertices.begin(), b.vertices.end());
    return convex_hull(std::move(combined));
}

template<Scalar T>
Result<Polygon<2, T>> union_of(const Polygon<2, T>& a, const Polygon<2, T>& b) {
    return a + b;
}

} // namespace spatium::geometry
