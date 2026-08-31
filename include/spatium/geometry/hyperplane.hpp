#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/core/error.hpp>
#  include <spatium/geometry/concepts.hpp>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium::geometry {

// Hyperplane in N-dimensional space: normal · p = offset
// In 3D this is a standard plane; in 2D a line; in N-D a codimension-1 surface.

template<std::size_t N, Scalar T = double>
struct Hyperplane {
    using ScalarType = T;
    using PointType = Vec<T, N>;
    static constexpr std::size_t ambient_dimension = N;

    Vec<T, N> normal; // unit normal
    T offset;         // signed distance from origin: normal · p = offset for points on plane

    // Factories
    static Result<Hyperplane> from_normal_and_point(Vec<T, N> n, PointType p) {
        auto len = n.norm();
        if (len < epsilon<T>())
            return std::unexpected(Error{ErrorCode::DegenerateInput, "zero normal"});
        auto unit = n / len;
        return Hyperplane{unit, unit.dot(p)};
    }

    // From 3 points (3D only)
    static Result<Hyperplane> from_points(PointType a, PointType b, PointType c)
        requires (N == 3)
    {
        auto ab = b - a;
        auto ac = c - a;
        auto n = ab.cross(ac);
        auto len = n.norm();
        if (len < epsilon<T>())
            return std::unexpected(Error{ErrorCode::DegenerateInput, "collinear points"});
        auto unit = n / len;
        return Hyperplane{unit, unit.dot(a)};
    }

    // From 2 points (2D only): the line through a and b
    static Result<Hyperplane> from_points(PointType a, PointType b)
        requires (N == 2)
    {
        auto ab = b - a;
        auto len = ab.norm();
        if (len < epsilon<T>())
            return std::unexpected(Error{ErrorCode::DegenerateInput, "coincident points"});
        // Normal is perpendicular to ab
        Vec<T, 2> n{-ab[1] / len, ab[0] / len};
        return Hyperplane{n, n.dot(a)};
    }

    // Signed distance: positive on normal side, negative on opposite
    constexpr T signed_distance(const PointType& p) const {
        return normal.dot(p) - offset;
    }

    T distance(const PointType& p) const {
        return std::abs(signed_distance(p));
    }

    constexpr PointType project(const PointType& p) const {
        return p - normal * signed_distance(p);
    }

    constexpr bool contains(const PointType& p, T eps = spatium::epsilon<T>()) const {
        return std::abs(signed_distance(p)) <= eps;
    }

    // Side test: +1 (normal side), -1 (opposite), 0 (on plane)
    constexpr int side(const PointType& p, T eps = spatium::epsilon<T>()) const {
        auto d = signed_distance(p);
        if (d > eps) return 1;
        if (d < -eps) return -1;
        return 0;
    }

    // Closest point on plane to origin
    constexpr PointType centroid() const {
        return normal * offset;
    }
};

// Convenience alias
template<Scalar T = double>
using Plane = Hyperplane<3, T>;

using Plane3 = Hyperplane<3>;
using Line2H = Hyperplane<2>; // 2D hyperplane = line

// Concept checks
static_assert(Shape<Hyperplane<3>>);
static_assert(DistanceQueryable<Hyperplane<3>>);

} // namespace spatium::geometry
