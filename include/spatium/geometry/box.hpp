#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/geometry/concepts.hpp>
#  include <algorithm>
#  include <span>
#endif

SPATIUM_EXPORT namespace spatium::geometry {

template<std::size_t N, Scalar T = double>
struct Box {
    using ScalarType = T;
    using PointType = Vec<T, N>;
    static constexpr std::size_t ambient_dimension = N;

    PointType min_corner{};
    PointType max_corner{};

    // Factories
    static constexpr Box from_center_half_extents(PointType center, Vec<T, N> half) {
        return {center - half, center + half};
    }

    static Box from_points(std::span<const PointType> points) {
        if (points.empty()) return {};
        Box result{points[0], points[0]};
        for (std::size_t i = 1; i < points.size(); ++i)
            result.expand(points[i]);
        return result;
    }

    // Dimensions
    constexpr Vec<T, N> extents() const { return max_corner - min_corner; }

    constexpr PointType centroid() const {
        return (min_corner + max_corner) * T{0.5};
    }

    // Measure: product of extents (length in 1D, area in 2D, volume in 3D, ...)
    constexpr T measure() const {
        auto ext = extents();
        T result{1};
        for (std::size_t i = 0; i < N; ++i)
            result *= ext[i];
        return result;
    }

    // Convenience aliases
    constexpr T volume() const requires (N == 3) { return measure(); }
    constexpr T area() const requires (N == 2) { return measure(); }
    constexpr T perimeter() const requires (N == 2) { return surface_measure(); }
    constexpr T surface_area() const requires (N == 3) { return surface_measure(); }

    // Surface measure: perimeter (2D), surface area (3D)
    constexpr T surface_measure() const requires (N == 2) {
        auto ext = extents();
        return T{2} * (ext[0] + ext[1]);
    }

    constexpr T surface_measure() const requires (N == 3) {
        auto ext = extents();
        return T{2} * (ext[0] * ext[1] + ext[1] * ext[2] + ext[0] * ext[2]);
    }

    // Containment
    constexpr bool contains(const PointType& p) const {
        for (std::size_t i = 0; i < N; ++i)
            if (p[i] < min_corner[i] || p[i] > max_corner[i])
                return false;
        return true;
    }

    constexpr bool contains(const Box& other) const {
        for (std::size_t i = 0; i < N; ++i)
            if (other.min_corner[i] < min_corner[i] || other.max_corner[i] > max_corner[i])
                return false;
        return true;
    }

    // Intersection test
    constexpr bool intersects(const Box& other) const {
        for (std::size_t i = 0; i < N; ++i)
            if (min_corner[i] > other.max_corner[i] || max_corner[i] < other.min_corner[i])
                return false;
        return true;
    }

    // Intersection result (empty if no overlap)
    constexpr std::optional<Box> intersection(const Box& other) const {
        Box result;
        for (std::size_t i = 0; i < N; ++i) {
            result.min_corner[i] = std::max(min_corner[i], other.min_corner[i]);
            result.max_corner[i] = std::min(max_corner[i], other.max_corner[i]);
            if (result.min_corner[i] > result.max_corner[i])
                return std::nullopt;
        }
        return result;
    }

    // Union (bounding box of both)
    constexpr Box union_with(const Box& other) const {
        Box result;
        for (std::size_t i = 0; i < N; ++i) {
            result.min_corner[i] = std::min(min_corner[i], other.min_corner[i]);
            result.max_corner[i] = std::max(max_corner[i], other.max_corner[i]);
        }
        return result;
    }

    // Distance from point to box (0 if inside)
    T distance(const PointType& p) const {
        auto clamped = project(p);
        return (p - clamped).norm();
    }

    // Project point onto box (clamp to boundary)
    constexpr PointType project(const PointType& p) const {
        PointType result;
        for (std::size_t i = 0; i < N; ++i)
            result[i] = std::clamp(p[i], min_corner[i], max_corner[i]);
        return result;
    }

    // Bounding box of self (identity, for concept satisfaction)
    constexpr Box bounding_box() const { return *this; }

    constexpr bool operator==(const Box&) const = default;

private:
    constexpr void expand(const PointType& p) {
        for (std::size_t i = 0; i < N; ++i) {
            min_corner[i] = std::min(min_corner[i], p[i]);
            max_corner[i] = std::max(max_corner[i], p[i]);
        }
    }
};

// Concept verification
static_assert(Shape<Box<3>>);
static_assert(ClosedShape<Box<3>>);
static_assert(Measurable<Box<3>>);
static_assert(Bounded<Box<3>>);
static_assert(DistanceQueryable<Box<3>>);

// Aliases
using Box2 = Box<2>;
using Box3 = Box<3>;

} // namespace spatium::geometry
