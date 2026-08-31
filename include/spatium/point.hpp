#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/algebra/format.hpp>
#  include <ostream>
#endif

SPATIUM_EXPORT namespace spatium {

// Type-safe point wrapper. Prevents mixing points from different spaces.
// Lightweight: same layout as the raw point type.

template<Set S>
struct Point {
    using Space = S;
    using Raw = typename S::PointType;
    using ScalarType = typename S::ScalarType;

    Raw coords;

    constexpr Point() = default;
    constexpr explicit Point(Raw c) : coords(std::move(c)) {}

    // Implicit access to raw coords
    constexpr const Raw& raw() const { return coords; }
    constexpr Raw& raw() { return coords; }

    constexpr bool operator==(const Point&) const = default;

    // Distance to another point (requires MetricSpace and a space instance)
    ScalarType distance_to(const Point& other, const S& space) const
        requires MetricSpace<S>
    {
        return space.distance(coords, other.coords);
    }

    // Move along tangent vector (requires Manifold)
    Point exp(const typename S::TangentVector& v, ScalarType t, const S& space) const
        requires Manifold<S>
    {
        return Point{space.exp_map(coords, v, t)};
    }

    // Tangent vector to another point (requires Manifold)
    auto log(const Point& other, const S& space) const
        requires Manifold<S>
    {
        return space.log_map(coords, other.coords);
    }

    // Stream output
    friend std::ostream& operator<<(std::ostream& os, const Point& p) {
        return os << p.coords;
    }
};

// Factory: pt<Space>(coords)
template<Set S>
constexpr Point<S> pt(typename S::PointType coords) {
    return Point<S>{std::move(coords)};
}

} // namespace spatium
