#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/core/error.hpp>
#  include <spatium/geometry/box.hpp>
#  include <spatium/geometry/circle.hpp>
#  include <spatium/geometry/concepts.hpp>
#  include <spatium/geometry/hyperplane.hpp>
#  include <spatium/geometry/line.hpp>
#  include <spatium/geometry/polygon.hpp>
#  include <spatium/geometry/triangle.hpp>
#  include <algorithm>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium::geometry {

// ── Point clip (containment gate) ─────────────────────────────
// Returns the point unchanged if inside the shape, NoIntersection otherwise.

template<std::size_t N, Scalar T>
Result<Vec<T, N>> clip(const Vec<T, N>& p, const Triangle<N, T>& tri) {
    if (tri.contains(p)) return p;
    return std::unexpected(Error{ErrorCode::NoIntersection, "point outside triangle"});
}

template<std::size_t N, Scalar T>
Result<Vec<T, N>> clip(const Vec<T, N>& p, const Segment<N, T>& seg) {
    auto ab = seg.b - seg.a;
    auto len_sq = ab.norm_squared();
    if (len_sq < epsilon<T>() * epsilon<T>()) {
        if ((p - seg.a).norm_squared() < epsilon<T>() * epsilon<T>())
            return p;
        return std::unexpected(Error{ErrorCode::NoIntersection, "point not on degenerate segment"});
    }
    auto t = (p - seg.a).dot(ab) / len_sq;
    if (t < -epsilon<T>() || t > T{1} + epsilon<T>())
        return std::unexpected(Error{ErrorCode::NoIntersection, "point outside segment"});
    return p;
}

template<std::size_t N, Scalar T>
Result<Vec<T, N>> clip(const Vec<T, N>& p, const Disk<N, T>& disk) {
    if (disk.contains(p)) return p;
    return std::unexpected(Error{ErrorCode::NoIntersection, "point outside disk"});
}

template<std::size_t N, Scalar T>
Result<Vec<T, N>> clip(const Vec<T, N>& p, const Box<N, T>& box) {
    if (box.contains(p)) return p;
    return std::unexpected(Error{ErrorCode::NoIntersection, "point outside box"});
}

template<std::size_t N, Scalar T>
Result<Vec<T, N>> clip(const Vec<T, N>& p, const Polygon<N, T>& poly) {
    if (poly.contains(p)) return p;
    return std::unexpected(Error{ErrorCode::NoIntersection, "point outside polygon"});
}

// ── Line clip (line → segment within shape bounds) ────────────

// Clip a line to a triangle: find where the line enters and exits the triangle.
// Works by intersecting the line with each triangle edge (as a segment).
template<Scalar T>
Result<Segment<3, T>> clip(const Line<3, T>& line, const Triangle<3, T>& tri) {
    // The line must lie in the triangle's plane for a segment result.
    // If not coplanar, intersection is at most a point — handled elsewhere.
    // For coplanar case: intersect line with each triangle edge,
    // collect the parameter range.

    auto plane = tri.supporting_plane();
    if (!plane) return std::unexpected(Error{ErrorCode::DegenerateShape, "degenerate triangle"});

    // Check line lies in the plane
    auto dist = std::abs(plane->signed_distance(line.origin));
    auto para = std::abs(line.direction.dot(plane->normal));
    if (dist > epsilon<T>() || para > epsilon<T>())
        return std::unexpected(Error{ErrorCode::NoIntersection, "line not in triangle plane"});

    // Find t range by clipping against each edge's half-plane (2D in triangle plane)
    T tmin = std::numeric_limits<T>::lowest();
    T tmax = std::numeric_limits<T>::max();

    for (std::size_t i = 0; i < 3; ++i) {
        auto& va = tri[i];
        auto& vb = tri[(i + 1) % 3];
        auto edge = vb - va;
        auto edge_normal = edge.cross(plane->normal); // inward-facing normal of edge
        auto en_len = edge_normal.norm();
        if (en_len < epsilon<T>()) continue;
        edge_normal = edge_normal / en_len;

        // Ensure edge normal points inward (toward opposite vertex)
        auto& vc = tri[(i + 2) % 3];
        if (edge_normal.dot(vc - va) < T{0})
            edge_normal = edge_normal * T{-1};

        auto denom = line.direction.dot(edge_normal);
        auto numer = edge_normal.dot(va - line.origin);

        if (std::abs(denom) < epsilon<T>()) {
            // Line parallel to edge — check which side
            // numer = en·(va - origin); if positive, origin is outside this edge
            if (numer > epsilon<T>())
                return std::unexpected(Error{ErrorCode::NoIntersection, "line outside triangle"});
        } else {
            auto t = numer / denom;
            if (denom > T{0})
                tmin = std::max(tmin, t); // entering
            else
                tmax = std::min(tmax, t); // exiting
        }
    }

    if (tmin > tmax + epsilon<T>())
        return std::unexpected(Error{ErrorCode::NoIntersection, "line misses triangle"});

    return Segment<3, T>{line.at(tmin), line.at(tmax)};
}

// Clip a line to a segment: project segment endpoints onto line, return overlap.
template<std::size_t N, Scalar T>
Result<Segment<N, T>> clip(const Line<N, T>& line, const Segment<N, T>& seg) {
    // Check collinearity
    auto d_a = (seg.a - line.project(seg.a)).norm();
    auto d_b = (seg.b - line.project(seg.b)).norm();
    if (d_a > epsilon<T>() || d_b > epsilon<T>())
        return std::unexpected(Error{ErrorCode::NoIntersection, "segment not on line"});

    auto t_a = line.parameter(seg.a);
    auto t_b = line.parameter(seg.b);
    if (t_a > t_b) std::swap(t_a, t_b);
    return Segment<N, T>{line.at(t_a), line.at(t_b)};
}

// Clip a line to a disk: find the chord where line crosses the disk.
template<Scalar T>
Result<Segment<3, T>> clip(const Line<3, T>& line, const Disk<3, T>& disk) {
    auto& c = disk.boundary.center;
    auto& n = disk.boundary.normal;
    auto r = disk.boundary.radius;

    // Line must lie in the disk's plane
    auto dist = std::abs(n.dot(line.origin - c));
    auto para = std::abs(line.direction.dot(n));
    if (dist > epsilon<T>() || para > epsilon<T>())
        return std::unexpected(Error{ErrorCode::NoIntersection, "line not in disk plane"});

    // Project center onto line, find half-chord length
    auto oc = c - line.origin;
    auto t_center = oc.dot(line.direction);
    auto closest = line.at(t_center);
    auto d_sq = (closest - c).norm_squared();
    auto r_sq = r * r;

    if (d_sq > r_sq + epsilon<T>())
        return std::unexpected(Error{ErrorCode::NoIntersection, "line misses disk"});

    auto half = std::sqrt(std::max(r_sq - d_sq, T{0}));
    return Segment<3, T>{line.at(t_center - half), line.at(t_center + half)};
}

// ── Segment reclip (segment → shorter segment) ────────────────
// Clip an existing segment to the bounds of a shape.

template<std::size_t N, Scalar T>
Result<Segment<N, T>> clip(const Segment<N, T>& seg, const Triangle<N, T>& tri) {
    // Parametrize seg as a + t*(b-a), t ∈ [0,1].
    // Find the sub-interval where points are inside the triangle.
    auto dir = seg.b - seg.a;
    auto len = dir.norm();
    if (len < epsilon<T>())
        return clip(seg.a, tri).transform([&](auto&&) { return seg; });

    // Use line clip, then restrict to [0,1]
    auto line = Line<N, T>{seg.a, dir / len};
    auto clipped = clip(line, tri);
    if (!clipped) return std::unexpected(clipped.error());

    auto t_a = line.parameter(clipped->a);
    auto t_b = line.parameter(clipped->b);
    auto seg_len = len;
    // Convert from line parameter to segment parameter [0, len] → [0, 1]
    T s_lo = std::max(t_a / seg_len, T{0});
    T s_hi = std::min(t_b / seg_len, T{1});

    if (s_lo > s_hi + epsilon<T>())
        return std::unexpected(Error{ErrorCode::NoIntersection, "segment outside triangle"});

    return Segment<N, T>{seg.at(s_lo), seg.at(s_hi)};
}

template<Scalar T>
Result<Segment<3, T>> clip(const Segment<3, T>& seg, const Disk<3, T>& disk) {
    auto dir = seg.b - seg.a;
    auto len = dir.norm();
    if (len < epsilon<T>())
        return clip(seg.a, disk).transform([&](auto&&) { return seg; });

    auto line = Line<3, T>{seg.a, dir / len};
    auto clipped = clip(line, disk);
    if (!clipped) return std::unexpected(clipped.error());

    auto t_a = line.parameter(clipped->a);
    auto t_b = line.parameter(clipped->b);
    T s_lo = std::max(t_a / len, T{0});
    T s_hi = std::min(t_b / len, T{1});

    if (s_lo > s_hi + epsilon<T>())
        return std::unexpected(Error{ErrorCode::NoIntersection, "segment outside disk"});

    return Segment<3, T>{seg.at(s_lo), seg.at(s_hi)};
}

} // namespace spatium::geometry
