#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/core/error.hpp>
#  include <spatium/geometry/concepts.hpp>
#  include <algorithm>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium::geometry {

// Forward declaration
template<std::size_t N, Scalar T> struct Box;

// ── Line ───────────────────────────────────────────────────────
// Infinite line: origin + t * direction, t ∈ (-∞, +∞)

template<std::size_t N, Scalar T = double>
struct Line {
    using ScalarType = T;
    using PointType = Vec<T, N>;
    static constexpr std::size_t ambient_dimension = N;

    PointType origin;
    PointType direction; // unit vector

    static Result<Line> from(PointType origin, PointType dir) {
        auto len = dir.norm();
        if (len < epsilon<T>())
            return std::unexpected(Error{ErrorCode::DegenerateInput, "zero direction"});
        return Line{origin, dir / len};
    }

    constexpr PointType at(T t) const { return origin + direction * t; }

    constexpr PointType centroid() const { return origin; }

    // Parameter t of the closest point on line to p
    constexpr T parameter(const PointType& p) const {
        return (p - origin).dot(direction);
    }

    constexpr PointType project(const PointType& p) const {
        return at(parameter(p));
    }

    T distance(const PointType& p) const {
        return (p - project(p)).norm();
    }
};

// ── Ray ────────────────────────────────────────────────────────
// Half-line: origin + t * direction, t ∈ [0, +∞)

template<std::size_t N, Scalar T = double>
struct Ray {
    using ScalarType = T;
    using PointType = Vec<T, N>;
    static constexpr std::size_t ambient_dimension = N;

    PointType origin;
    PointType direction; // unit vector

    static Result<Ray> from(PointType origin, PointType dir) {
        auto len = dir.norm();
        if (len < epsilon<T>())
            return std::unexpected(Error{ErrorCode::DegenerateInput, "zero direction"});
        return Ray{origin, dir / len};
    }

    constexpr PointType at(T t) const { return origin + direction * std::max(t, T{0}); }

    constexpr PointType centroid() const { return origin; }

    constexpr PointType project(const PointType& p) const {
        auto t = std::max((p - origin).dot(direction), T{0});
        return origin + direction * t;
    }

    T distance(const PointType& p) const {
        return (p - project(p)).norm();
    }

    // Subspace: the line this ray is a bounded region of.
    Result<Line<N, T>> subspace() const {
        return Line<N, T>::from(origin, direction);
    }
};

// ── Segment ────────────────────────────────────────────────────
// Finite line: a + t * (b - a), t ∈ [0, 1]

template<std::size_t N, Scalar T = double>
struct Segment {
    using ScalarType = T;
    using PointType = Vec<T, N>;
    static constexpr std::size_t ambient_dimension = N;

    PointType a, b;

    T measure() const { return (b - a).norm(); }

    constexpr T length_squared() const { return (b - a).norm_squared(); }

    constexpr PointType midpoint() const { return (a + b) * T{0.5}; }

    constexpr PointType at(T t) const { return a + (b - a) * t; }

    constexpr PointType centroid() const { return midpoint(); }

    T length() const { return measure(); }

    constexpr PointType project(const PointType& p) const {
        auto ab = b - a;
        auto len_sq = ab.norm_squared();
        if (len_sq < epsilon<T>() * epsilon<T>()) return a; // degenerate
        auto t = std::clamp((p - a).dot(ab) / len_sq, T{0}, T{1});
        return a + ab * t;
    }

    T distance(const PointType& p) const {
        return (p - project(p)).norm();
    }

    // Subspace: the line this segment is a bounded region of.
    Result<Line<N, T>> subspace() const {
        return Line<N, T>::from(a, b - a);
    }

    Box<N, T> bounding_box() const;
};

// Aliases
using Line2 = Line<2>;
using Line3 = Line<3>;
using Ray2 = Ray<2>;
using Ray3 = Ray<3>;
using Segment2 = Segment<2>;
using Segment3 = Segment<3>;

// Concept checks
static_assert(Shape<Line<3>>);
static_assert(DistanceQueryable<Line<3>>);
static_assert(Shape<Ray<3>>);
static_assert(DistanceQueryable<Ray<3>>);
static_assert(Shape<Segment<3>>);
static_assert(DistanceQueryable<Segment<3>>);
static_assert(BoundedRegion<Ray<3>>);
static_assert(BoundedRegion<Segment<3>>);

} // namespace spatium::geometry

// Deferred implementation (needs box.hpp)
#include <spatium/geometry/box.hpp>

namespace spatium::geometry {

template<std::size_t N, Scalar T>
Box<N, T> Segment<N, T>::bounding_box() const {
    Box<N, T> result;
    for (std::size_t i = 0; i < N; ++i) {
        result.min_corner[i] = std::min(a[i], b[i]);
        result.max_corner[i] = std::max(a[i], b[i]);
    }
    return result;
}

} // namespace spatium::geometry
