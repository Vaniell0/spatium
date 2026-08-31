#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/geometry/concepts.hpp>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium::geometry {

// Turns any geometric Shape into a Surface (space).
// Requirements: the Shape must support project(point) and normal(point).
//
// This bridges geometry/ and spaces/: a Triangle becomes a 2D manifold,
// a Disk becomes a flat 2D surface, any mesh face becomes navigable.
//
// Geodesics on flat surfaces are straight lines projected back onto the surface.
// For curved shapes (if project is nonlinear), geodesics are approximate.

// Shape with normal — either normal() or normal(point)
template<typename S>
concept HasPointNormal = requires(const S& s, const typename S::PointType& p) {
    s.normal(p);
};

template<typename S>
concept HasConstNormal = requires(const S& s) {
    s.normal();
};

template<typename S>
concept SurfaceCapable = DistanceQueryable<S>
    && (HasPointNormal<S> || HasConstNormal<S>);

template<SurfaceCapable S>
struct ShapeSurface {
    using T = typename S::ScalarType;
    using ScalarType = T;
    using PointType = typename S::PointType;
    using TangentVector = PointType;

    static constexpr std::size_t dimension = S::ambient_dimension - 1;
    static constexpr bool is_complete = false; // bounded shapes are not complete

    S shape;

    // TopologicalSpace
    bool contains(const PointType& p) const {
        return shape.distance(p) < epsilon<T>();
    }

    // MetricSpace: distance along surface
    // For flat shapes = euclidean distance between projected points
    // For curved = approximation
    ScalarType distance(const PointType& a, const PointType& b) const {
        auto pa = shape.project(a);
        auto pb = shape.project(b);
        return (pa - pb).norm();
    }

    // Manifold: exp_map — walk along tangent direction, project back
    PointType exp_map(const PointType& p, const TangentVector& v, ScalarType t) const {
        return shape.project(p + v * t);
    }

    // Manifold: log_map — tangent vector from p toward q, in tangent plane at p
    TangentVector log_map(const PointType& p, const PointType& q) const {
        auto pp = shape.project(p);
        auto pq = shape.project(q);
        TangentVector diff = pq - pp;
        auto n = get_normal(pp);
        auto n_len_sq = n.dot(n);
        if (n_len_sq > epsilon<T>())
            diff = TangentVector{diff - n * (diff.dot(n) / n_len_sq)};
        return diff;
    }

    // RiemannianManifold: flat metric on tangent space
    ScalarType metric_at(const PointType&,
                         const TangentVector& u,
                         const TangentVector& v) const {
        return u.dot(v);
    }

    // Surface
    PointType project(const PointType& p) const {
        return shape.project(p);
    }

    TangentVector normal(const PointType& p) const {
        return get_normal(p);
    }

private:
    PointType get_normal(const PointType& p) const {
        if constexpr (HasPointNormal<S>)
            return shape.normal(p);
        else
            return shape.normal();
    }
};

// Factory: shape → surface
template<SurfaceCapable S>
ShapeSurface<S> as_surface(S shape) {
    return {std::move(shape)};
}

} // namespace spatium::geometry
