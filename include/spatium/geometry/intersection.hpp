#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/core/error.hpp>
#  include <spatium/geometry/box.hpp>
#  include <spatium/geometry/hyperplane.hpp>
#  include <spatium/geometry/line.hpp>
#  include <spatium/geometry/triangle.hpp>
#  include <spatium/geometry/clip.hpp>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium::geometry {

// ── Line-Hyperplane ────────────────────────────────────────────
// Returns the intersection point, or NoIntersection if parallel.

template<std::size_t N, Scalar T>
Result<Vec<T, N>> intersect(const Line<N, T>& line, const Hyperplane<N, T>& plane) {
    auto denom = line.direction.dot(plane.normal);
    if (std::abs(denom) < epsilon<T>())
        return std::unexpected(Error{ErrorCode::NoIntersection, "line parallel to plane"});
    auto t = (plane.offset - plane.normal.dot(line.origin)) / denom;
    return line.at(t);
}

// ── Ray-Hyperplane ─────────────────────────────────────────────

template<std::size_t N, Scalar T>
Result<Vec<T, N>> intersect(const Ray<N, T>& ray, const Hyperplane<N, T>& plane) {
    auto denom = ray.direction.dot(plane.normal);
    if (std::abs(denom) < epsilon<T>())
        return std::unexpected(Error{ErrorCode::NoIntersection, "ray parallel to plane"});
    auto t = (plane.offset - plane.normal.dot(ray.origin)) / denom;
    if (t < T{0})
        return std::unexpected(Error{ErrorCode::NoIntersection, "intersection behind ray"});
    return ray.origin + ray.direction * t;
}

// ── Segment-Hyperplane ─────────────────────────────────────────

template<std::size_t N, Scalar T>
Result<Vec<T, N>> intersect(const Segment<N, T>& seg, const Hyperplane<N, T>& plane) {
    auto dir = seg.b - seg.a;
    auto denom = dir.dot(plane.normal);
    if (std::abs(denom) < epsilon<T>())
        return std::unexpected(Error{ErrorCode::NoIntersection, "segment parallel to plane"});
    auto t = (plane.offset - plane.normal.dot(seg.a)) / denom;
    if (t < T{0} || t > T{1})
        return std::unexpected(Error{ErrorCode::NoIntersection, "intersection outside segment"});
    return seg.a + dir * t;
}

// ── Line-Line 2D ───────────────────────────────────────────────

template<Scalar T>
Result<Vec<T, 2>> intersect(const Line<2, T>& a, const Line<2, T>& b) {
    auto cross = a.direction[0] * b.direction[1] - a.direction[1] * b.direction[0];
    if (std::abs(cross) < epsilon<T>())
        return std::unexpected(Error{ErrorCode::NoIntersection, "parallel lines"});
    auto d = b.origin - a.origin;
    auto t = (d[0] * b.direction[1] - d[1] * b.direction[0]) / cross;
    return a.at(t);
}

// ── Segment-Segment 2D ────────────────────────────────────────

template<Scalar T>
Result<Vec<T, 2>> intersect(const Segment<2, T>& s1, const Segment<2, T>& s2) {
    auto d1 = s1.b - s1.a;
    auto d2 = s2.b - s2.a;
    auto cross = d1[0] * d2[1] - d1[1] * d2[0];
    if (std::abs(cross) < epsilon<T>())
        return std::unexpected(Error{ErrorCode::NoIntersection, "parallel segments"});

    auto d = s2.a - s1.a;
    auto t1 = (d[0] * d2[1] - d[1] * d2[0]) / cross;
    auto t2 = (d[0] * d1[1] - d[1] * d1[0]) / cross;

    if (t1 < T{0} || t1 > T{1} || t2 < T{0} || t2 > T{1})
        return std::unexpected(Error{ErrorCode::NoIntersection, "segments don't overlap"});

    return s1.a + d1 * t1;
}

// ── Ray-Triangle 3D (Moller-Trumbore) ──────────────────────────

template<Scalar T>
struct RayTriHit {
    Vec<T, 3> point;
    T t;
    T u;                // barycentric: vertex 1 weight
    T v;                // barycentric: vertex 2 weight (w = 1-u-v for vertex 0)
    Vec<T, 3> normal;   // face normal (right-hand rule, unnormalized length != 1)
};

template<Scalar T>
Result<RayTriHit<T>> ray_triangle(const Ray<3, T>& ray, const Triangle<3, T>& tri) {
    auto e1 = tri[1] - tri[0];
    auto e2 = tri[2] - tri[0];
    auto h = ray.direction.cross(e2);
    auto a = e1.dot(h);

    if (std::abs(a) < epsilon<T>())
        return std::unexpected(Error{ErrorCode::NoIntersection, "ray parallel to triangle"});

    auto f = T{1} / a;
    auto s = ray.origin - tri[0];
    auto u = f * s.dot(h);
    if (u < T{0} || u > T{1})
        return std::unexpected(Error{ErrorCode::NoIntersection});

    auto q = s.cross(e1);
    auto v = f * ray.direction.dot(q);
    if (v < T{0} || u + v > T{1})
        return std::unexpected(Error{ErrorCode::NoIntersection});

    auto t = f * e2.dot(q);
    if (t < T{0})
        return std::unexpected(Error{ErrorCode::NoIntersection, "intersection behind ray"});

    auto n = e1.cross(e2);
    auto inv_n = T{1} / n.norm();
    return RayTriHit<T>{ray.origin + ray.direction * t, t, u, v, n * inv_n};
}

template<Scalar T>
Result<Vec<T, 3>> intersect(const Ray<3, T>& ray, const Triangle<3, T>& tri) {
    auto hit = ray_triangle(ray, tri);
    if (!hit) return std::unexpected(hit.error());
    return hit->point;
}

// ── Line-Triangle 3D ───────────────────────────────────────────

template<Scalar T>
Result<Vec<T, 3>> intersect(const Line<3, T>& line, const Triangle<3, T>& tri) {
    auto e1 = tri[1] - tri[0];
    auto e2 = tri[2] - tri[0];
    auto h = line.direction.cross(e2);
    auto a = e1.dot(h);

    if (std::abs(a) < epsilon<T>())
        return std::unexpected(Error{ErrorCode::NoIntersection, "line parallel to triangle"});

    auto f = T{1} / a;
    auto s = line.origin - tri[0];
    auto u = f * s.dot(h);
    if (u < T{0} || u > T{1})
        return std::unexpected(Error{ErrorCode::NoIntersection});

    auto q = s.cross(e1);
    auto v = f * line.direction.dot(q);
    if (v < T{0} || u + v > T{1})
        return std::unexpected(Error{ErrorCode::NoIntersection});

    auto t = f * e2.dot(q);
    return line.at(t);
}

// ── Ray-Box (slab method) ──────────────────────────────────────

template<std::size_t N, Scalar T>
Result<std::pair<T, T>> intersect_parameters(const Ray<N, T>& ray, const Box<N, T>& box) {
    T tmin = T{0};
    T tmax = std::numeric_limits<T>::max();

    for (std::size_t i = 0; i < N; ++i) {
        if (std::abs(ray.direction[i]) < epsilon<T>()) {
            if (ray.origin[i] < box.min_corner[i] || ray.origin[i] > box.max_corner[i])
                return std::unexpected(Error{ErrorCode::NoIntersection});
        } else {
            auto inv_d = T{1} / ray.direction[i];
            auto t1 = (box.min_corner[i] - ray.origin[i]) * inv_d;
            auto t2 = (box.max_corner[i] - ray.origin[i]) * inv_d;
            if (t1 > t2) std::swap(t1, t2);
            tmin = std::max(tmin, t1);
            tmax = std::min(tmax, t2);
            if (tmin > tmax)
                return std::unexpected(Error{ErrorCode::NoIntersection});
        }
    }
    return std::pair{tmin, tmax};
}

template<std::size_t N, Scalar T>
Result<Vec<T, N>> intersect(const Ray<N, T>& ray, const Box<N, T>& box) {
    auto params = intersect_parameters(ray, box);
    if (!params) return std::unexpected(params.error());
    return ray.origin + ray.direction * params->first;
}

// ── Hyperplane-Hyperplane 3D (returns line) ────────────────────

template<Scalar T>
Result<Line<3, T>> intersect(const Hyperplane<3, T>& a, const Hyperplane<3, T>& b) {
    auto dir = a.normal.cross(b.normal);
    auto len = dir.norm();
    if (len < epsilon<T>())
        return std::unexpected(Error{ErrorCode::NoIntersection, "parallel planes"});
    dir = dir / len;

    // Find a point on the line: solve the 2-plane system
    // n1·p = d1, n2·p = d2
    // Pick the component with largest |dir| and set it to 0
    std::size_t max_i = 0;
    for (std::size_t i = 1; i < 3; ++i)
        if (std::abs(dir[i]) > std::abs(dir[max_i]))
            max_i = i;

    Vec<T, 3> point{};
    // Solve 2x2 system in the other two axes
    std::size_t i0 = (max_i + 1) % 3;
    std::size_t i1 = (max_i + 2) % 3;
    auto det = a.normal[i0] * b.normal[i1] - a.normal[i1] * b.normal[i0];
    point[i0] = (a.offset * b.normal[i1] - b.offset * a.normal[i1]) / det;
    point[i1] = (a.normal[i0] * b.offset - b.normal[i0] * a.offset) / det;
    point[max_i] = T{0};

    return Line<3, T>{point, dir};
}

// ── Triangle-Triangle 3D ───────────────────────────────────────
// Returns the intersection segment, or NoIntersection.

template<Scalar T>
Result<Segment<3, T>> intersect(const Triangle<3, T>& t1, const Triangle<3, T>& t2) {
    // Get planes of both triangles
    auto plane1 = t1.supporting_plane();
    auto plane2 = t2.supporting_plane();
    if (!plane1 || !plane2)
        return std::unexpected(Error{ErrorCode::DegenerateShape, "degenerate triangle"});

    // Intersect the two planes → line
    auto line_result = intersect(*plane1, *plane2);
    if (!line_result)
        return std::unexpected(Error{ErrorCode::NoIntersection, "coplanar triangles"});
    auto& line = *line_result;

    // Project triangle vertices onto the line → parameter intervals
    auto project_tri = [&](const Triangle<3, T>& tri) -> std::optional<std::pair<T, T>> {
        std::array<T, 3> params;
        std::array<T, 3> dists;
        const auto& plane_of_other = (&tri == &t1) ? *plane2 : *plane1;

        for (int i = 0; i < 3; ++i) {
            params[i] = line.parameter(tri[i]);
            dists[i] = plane_of_other.signed_distance(tri[i]);
        }

        // Find interval where triangle crosses the other plane
        T tmin = std::numeric_limits<T>::max();
        T tmax = std::numeric_limits<T>::lowest();

        for (int i = 0; i < 3; ++i) {
            int j = (i + 1) % 3;
            if ((dists[i] > T{0}) != (dists[j] > T{0}) &&
                std::abs(dists[i] - dists[j]) > epsilon<T>()) {
                T t_edge = params[i] + (params[j] - params[i]) * dists[i] / (dists[i] - dists[j]);
                tmin = std::min(tmin, t_edge);
                tmax = std::max(tmax, t_edge);
            }
            // Vertex on plane
            using std::abs;
            if (abs(dists[i]) < epsilon<T>()) {
                tmin = std::min(tmin, params[i]);
                tmax = std::max(tmax, params[i]);
            }
        }

        if (tmin > tmax) return std::nullopt;
        return std::pair{tmin, tmax};
    };

    auto interval1 = project_tri(t1);
    auto interval2 = project_tri(t2);
    if (!interval1 || !interval2)
        return std::unexpected(Error{ErrorCode::NoIntersection});

    // Intersect intervals
    T lo = std::max(interval1->first, interval2->first);
    T hi = std::min(interval1->second, interval2->second);
    if (lo > hi + epsilon<T>())
        return std::unexpected(Error{ErrorCode::NoIntersection});

    return Segment<3, T>{line.at(lo), line.at(hi)};
}

// ── Generic BoundedRegion intersect ───────────────────────────
// Decompose: intersect subspaces, then clip to both shapes.
// Only matches when no specialized overload exists (weaker constraint).

// Case 1: Both subspaces are Hyperplanes → intersection is a Line → clip to Segments.
template<typename A, typename B>
    requires BoundedRegion<A> && BoundedRegion<B>
          && std::same_as<subspace_t<A>, Hyperplane<A::ambient_dimension, typename A::ScalarType>>
          && std::same_as<subspace_t<B>, Hyperplane<B::ambient_dimension, typename B::ScalarType>>
          && (A::ambient_dimension == 3) && (B::ambient_dimension == 3)
          && std::same_as<typename A::ScalarType, typename B::ScalarType>
auto intersect_via_subspace(const A& a, const B& b)
    -> Result<Segment<3, typename A::ScalarType>>
{
    auto sa = a.subspace();
    auto sb = b.subspace();
    if (!sa || !sb)
        return std::unexpected(Error{ErrorCode::DegenerateShape, "degenerate shape"});

    auto line = intersect(*sa, *sb);
    if (!line)
        return std::unexpected(line.error());

    auto clipped_a = clip(*line, a);
    if (!clipped_a)
        return std::unexpected(clipped_a.error());

    auto clipped_b = clip(*clipped_a, b);
    if (!clipped_b)
        return std::unexpected(clipped_b.error());

    return *clipped_b;
}

// Case 2: One subspace is Line, other is Hyperplane → intersection is a Point → clip to both.
template<typename A, typename B>
    requires BoundedRegion<A> && BoundedRegion<B>
          && std::same_as<subspace_t<A>, Line<A::ambient_dimension, typename A::ScalarType>>
          && std::same_as<subspace_t<B>, Hyperplane<B::ambient_dimension, typename B::ScalarType>>
          && (A::ambient_dimension == B::ambient_dimension)
          && std::same_as<typename A::ScalarType, typename B::ScalarType>
auto intersect_via_subspace(const A& a, const B& b)
    -> Result<Vec<typename A::ScalarType, A::ambient_dimension>>
{
    auto sa = a.subspace();
    auto sb = b.subspace();
    if (!sa || !sb)
        return std::unexpected(Error{ErrorCode::DegenerateShape, "degenerate shape"});

    auto point = intersect(*sa, *sb);
    if (!point)
        return std::unexpected(point.error());

    auto ca = clip(*point, a);
    if (!ca)
        return std::unexpected(ca.error());

    auto cb = clip(*ca, b);
    if (!cb)
        return std::unexpected(cb.error());

    return *cb;
}

// Symmetric: Hyperplane-subspace first, Line-subspace second
template<typename A, typename B>
    requires BoundedRegion<A> && BoundedRegion<B>
          && std::same_as<subspace_t<A>, Hyperplane<A::ambient_dimension, typename A::ScalarType>>
          && std::same_as<subspace_t<B>, Line<B::ambient_dimension, typename B::ScalarType>>
          && (A::ambient_dimension == B::ambient_dimension)
          && std::same_as<typename A::ScalarType, typename B::ScalarType>
auto intersect_via_subspace(const A& a, const B& b)
    -> Result<Vec<typename A::ScalarType, A::ambient_dimension>>
{
    return intersect_via_subspace(b, a);
}

} // namespace spatium::geometry
