#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <concepts>
#endif

SPATIUM_EXPORT namespace spatium::geometry {

// NOTE: All geometry types and operations (intersection, distance, area, etc.)
// assume Euclidean space. They use Vec dot product, cross product, and Euclidean
// norm internally. Do NOT use them with non-Euclidean metrics — results will be
// meaningless. For operations on curved spaces, use the Space interface
// (exp_map, log_map, distance) or ShapeSurface adapter instead.

// ── Shape sub-requirements ─────────────────────────────────────
// Each named concept below checks exactly one structural requirement.
// Composing them via && in Shape/ClosedShape/... gives precise
// diagnostics: failure reports the specific HasXxx<MyType> that did
// not hold, not "associated constraints not satisfied".

template<typename S>
concept HasShapeTypes = requires {
    typename S::ScalarType;
    typename S::PointType;
    { S::ambient_dimension } -> std::convertible_to<std::size_t>;
};

template<typename S>
concept HasCentroid = requires { typename S::PointType; }
    && requires(const S& s) {
        { s.centroid() } -> std::convertible_to<typename S::PointType>;
    };

template<typename S>
concept HasContains = requires { typename S::PointType; }
    && requires(const S& s, const typename S::PointType& p) {
        { s.contains(p) } -> std::convertible_to<bool>;
    };

// `measure()` is the dimension-generic name picked by the concept so that
// Segment (length), Triangle/Polygon/Disk (area), Box (volume) and
// arbitrary Simplex<N, K> all compose into the same generic algorithms
// (boolean ops, integration kernels, region builders).  Concrete 2D
// shapes additionally expose `area()` as a convenience alias because that
// is what users reach for when they write Triangle/Polygon code by hand;
// 3D shapes will mirror this with `volume()` when added.  The aliases
// must always forward to `measure()` — never duplicate the formula —
// so the two views stay in lock-step.
template<typename S>
concept HasMeasure = requires { typename S::ScalarType; }
    && requires(const S& s) {
        { s.measure() } -> std::convertible_to<typename S::ScalarType>;
    };

template<typename S>
concept HasBoundingBox = requires(const S& s) {
    { s.bounding_box() };
};

template<typename S>
concept HasPointDistance = requires { typename S::ScalarType; typename S::PointType; }
    && requires(const S& s, const typename S::PointType& p) {
        { s.distance(p) } -> std::convertible_to<typename S::ScalarType>;
    };

template<typename S>
concept HasProjectFromPoint = requires { typename S::PointType; }
    && requires(const S& s, const typename S::PointType& p) {
        { s.project(p) } -> std::convertible_to<typename S::PointType>;
    };

template<typename S>
concept HasSubspace = requires(const S& s) {
    { s.subspace() };  // returns Result<SubspaceType>
};

// A geometric shape with a centroid and known ambient dimension.
template<typename S>
concept Shape = HasShapeTypes<S> && HasCentroid<S>;

// Shape that supports point containment testing.
template<typename S>
concept ClosedShape = Shape<S> && HasContains<S>;

// Shape with a finite measure (length, area, volume, ...).
template<typename S>
concept Measurable = Shape<S> && HasMeasure<S>;

// Shape with a computable bounding box.
template<typename S>
concept Bounded = Shape<S> && HasBoundingBox<S>;

// Shape that supports distance and projection queries from a point.
template<typename S>
concept DistanceQueryable = Shape<S> && HasPointDistance<S> && HasProjectFromPoint<S>;

// Shape that is a bounded region of a lower-dimensional subspace.
// subspace() returns the containing subspace (e.g., Triangle → Plane, Segment → Line).
// Enables generic intersect: intersect(A.subspace(), B.subspace()) + clip to bounds.
template<typename S>
concept BoundedRegion = Shape<S> && HasSubspace<S>;

// Deduce the subspace type from a BoundedRegion.
template<BoundedRegion S>
using subspace_t = typename decltype(std::declval<const S&>().subspace())::value_type;

} // namespace spatium::geometry
