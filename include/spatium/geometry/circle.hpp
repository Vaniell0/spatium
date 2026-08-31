#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/geometry/box.hpp>
#  include <spatium/geometry/concepts.hpp>
#  include <spatium/geometry/hyperplane.hpp>
#  include <cmath>
#  include <numbers>
#endif

SPATIUM_EXPORT namespace spatium::geometry {

// Circle in N-dimensional space.
// In 2D: center + radius (no normal needed).
// In 3D+: center + radius + normal (defines orientation plane).

template<std::size_t N, Scalar T = double>
struct Circle {
    using ScalarType = T;
    using PointType = Vec<T, N>;
    static constexpr std::size_t ambient_dimension = N;

    PointType center;
    T radius;
    Vec<T, N> normal; // unit normal to the circle's plane (unused in 2D but still present)

    T circumference() const { return T{2} * std::numbers::pi_v<T> * radius; }
    T perimeter() const { return circumference(); }

    T measure() const { return circumference(); }

    constexpr PointType centroid() const { return center; }

    // Distance from point to the circle (curve, not disk)
    T distance(const PointType& p) const {
        if constexpr (N == 2) {
            auto r = (p - center).norm();
            return std::abs(r - radius);
        } else {
            // Project onto circle's plane, then compute distance to ring
            auto dp = p - center;
            auto along_normal = dp.dot(normal);
            auto in_plane = dp - normal * along_normal;
            auto r = in_plane.norm();
            auto radial_diff = r - radius;
            return std::sqrt(along_normal * along_normal + radial_diff * radial_diff);
        }
    }

    // Closest point on circle to p
    PointType project(const PointType& p) const {
        if constexpr (N == 2) {
            auto dp = p - center;
            auto r = dp.norm();
            if (r < epsilon<T>()) {
                // At center: pick arbitrary point on circle
                PointType result = center;
                result[0] += radius;
                return result;
            }
            return center + dp * (radius / r);
        } else {
            auto dp = p - center;
            auto along_normal = dp.dot(normal);
            auto in_plane = dp - normal * along_normal;
            auto r = in_plane.norm();
            if (r < epsilon<T>()) {
                // On the axis: pick arbitrary direction in plane
                // Find a vector perpendicular to normal
                Vec<T, N> perp{};
                for (std::size_t i = 0; i < N; ++i) {
                    if (std::abs(normal[i]) < T{0.9}) {
                        perp[i] = T{1};
                        break;
                    }
                }
                perp = perp - normal * normal.dot(perp);
                perp = perp / perp.norm();
                return center + perp * radius;
            }
            return center + in_plane * (radius / r);
        }
    }

    Box<N, T> bounding_box() const {
        // Conservative: axis-aligned box enclosing the circle
        Vec<T, N> half;
        if constexpr (N == 2) {
            half = Vec<T, N>{radius, radius};
        } else {
            // For each axis, the extent is radius * sin(angle between normal and axis)
            for (std::size_t i = 0; i < N; ++i)
                half[i] = radius * std::sqrt(T{1} - normal[i] * normal[i]);
        }
        return {center - half, center + half};
    }

    // Subspace: the plane this circle lies on.
    Result<Hyperplane<N, T>> subspace() const requires (N == 3) {
        return Hyperplane<N, T>::from_normal_and_point(normal, center);
    }
};

// Disk (filled circle)
template<std::size_t N, Scalar T = double>
struct Disk {
    using ScalarType = T;
    using PointType = Vec<T, N>;
    static constexpr std::size_t ambient_dimension = N;

    Circle<N, T> boundary;

    T measure() const { return std::numbers::pi_v<T> * boundary.radius * boundary.radius; }
    T perimeter() const { return boundary.circumference(); }

    T area() const { return measure(); }

    constexpr PointType centroid() const { return boundary.center; }

    bool contains(const PointType& p) const {
        if constexpr (N == 2) {
            return (p - boundary.center).norm_squared() <=
                   boundary.radius * boundary.radius + epsilon<T>();
        } else {
            auto dp = p - boundary.center;
            auto along = dp.dot(boundary.normal);
            if (std::abs(along) > epsilon<T>()) return false;
            auto in_plane = dp - boundary.normal * along;
            return in_plane.norm_squared() <=
                   boundary.radius * boundary.radius + epsilon<T>();
        }
    }

    T distance(const PointType& p) const {
        if constexpr (N == 2) {
            auto r = (p - boundary.center).norm();
            if (r <= boundary.radius) return T{0};
            return r - boundary.radius;
        } else {
            auto dp = p - boundary.center;
            auto along = dp.dot(boundary.normal);
            auto in_plane = dp - boundary.normal * along;
            auto r = in_plane.norm();
            if (r <= boundary.radius && std::abs(along) < epsilon<T>())
                return T{0}; // inside
            // Distance to nearest point on disk
            auto radial = std::min(r, boundary.radius);
            auto nearest = boundary.center + (r > epsilon<T>() ? in_plane * (radial / r) : in_plane);
            return (p - nearest).norm();
        }
    }

    PointType project(const PointType& p) const {
        if constexpr (N == 2) {
            auto dp = p - boundary.center;
            auto r = dp.norm();
            if (r <= boundary.radius) return p;
            return boundary.center + dp * (boundary.radius / r);
        } else {
            auto dp = p - boundary.center;
            auto along = dp.dot(boundary.normal);
            auto in_plane = dp - boundary.normal * along;
            auto r = in_plane.norm();
            if (r <= boundary.radius && std::abs(along) < epsilon<T>())
                return p; // already inside
            auto clamped_r = std::min(r, boundary.radius);
            if (r < epsilon<T>())
                return boundary.center;
            return boundary.center + in_plane * (clamped_r / r);
        }
    }

    Box<N, T> bounding_box() const { return boundary.bounding_box(); }

    // Subspace: the plane this disk lies on.
    Result<Hyperplane<N, T>> subspace() const requires (N == 3) {
        return boundary.subspace();
    }
};

// Concept checks
static_assert(Shape<Circle<3>>);
static_assert(Measurable<Circle<3>>);
static_assert(DistanceQueryable<Circle<3>>);
static_assert(Bounded<Circle<3>>);
static_assert(BoundedRegion<Circle<3>>);

static_assert(Shape<Disk<3>>);
static_assert(ClosedShape<Disk<3>>);
static_assert(Measurable<Disk<3>>);
static_assert(DistanceQueryable<Disk<3>>);
static_assert(Bounded<Disk<3>>);
static_assert(BoundedRegion<Disk<3>>);

// Aliases
using Circle2 = Circle<2>;
using Circle3 = Circle<3>;
using Disk2 = Disk<2>;
using Disk3 = Disk<3>;

} // namespace spatium::geometry
