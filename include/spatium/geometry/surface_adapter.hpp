#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/geometry/concepts.hpp>
#  include <spatium/spaces/sphere.hpp>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium::geometry {

// Turns any geometric Shape into a Surface (space).
// Requirements: the Shape must support project(point) and normal(point).
//
// This bridges geometry/ and spaces/: a Triangle becomes a 2D manifold,
// a Disk becomes a flat 2D surface, any mesh face becomes navigable.
//
// Geodesics on flat surfaces are straight lines projected back onto the
// surface -- exact, since for a flat shape the tangent-plane approximation
// IS the true geodesic. For curved shapes the same tangent-plane
// projection is only approximate, with one exception: a Quadric that is
// a round sphere (uniform radius, not a general ellipsoid) has an exact
// closed-form geodesic, and HasExactGeodesic below lets ShapeSurface
// detect that case and delegate exp_map/log_map to spaces::Sphere<2,T>
// instead. Torus, ellipsoid, cylinder and cone have no comparable simple
// closed form and stay on the tangent-plane approximation.

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

// Shape that can expose an exact spherical center + radius, letting
// ShapeSurface delegate exp_map/log_map to spaces::Sphere<2,T>'s
// closed-form great-circle formula instead of the generic tangent-plane
// approximation. is_round_sphere() is a runtime query rather than a
// type-level tag: Quadric<T> represents sphere/ellipsoid/cylinder/cone
// with the same C++ type (they differ only in the matrix it holds), so
// whether a *given instance* is a round sphere can't be decided at
// compile time -- only whether the type exposes the query at all.
template<typename S>
concept HasExactGeodesic = requires { typename S::PointType; typename S::ScalarType; }
    && requires(const S& s) {
        { s.is_round_sphere() } -> std::convertible_to<bool>;
        { s.sphere_center() } -> std::convertible_to<typename S::PointType>;
        { s.sphere_radius() } -> std::convertible_to<typename S::ScalarType>;
    };

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

    // Manifold: exp_map — walk along tangent direction, project back.
    // Exact for a round-sphere Quadric (delegates to spaces::Sphere<2,T>
    // after translating to the origin and scaling to unit radius);
    // tangent-plane-projection approximation otherwise.
    PointType exp_map(const PointType& p, const TangentVector& v, ScalarType t) const {
        if constexpr (HasExactGeodesic<S>) {
            if (shape.is_round_sphere()) {
                auto c = shape.sphere_center();
                auto r = shape.sphere_radius();
                spatium::Sphere<2, T> unit_sphere{};
                auto result0 = unit_sphere.exp_map((p - c) / r, v / r, t);
                return result0 * r + c;
            }
        }
        return shape.project(p + v * t);
    }

    // Manifold: log_map — tangent vector from p toward q, in tangent plane at p.
    // Exact for a round-sphere Quadric (delegates to spaces::Sphere<2,T>);
    // tangent-plane projection of the chord otherwise.
    TangentVector log_map(const PointType& p, const PointType& q) const {
        if constexpr (HasExactGeodesic<S>) {
            if (shape.is_round_sphere()) {
                auto c = shape.sphere_center();
                auto r = shape.sphere_radius();
                spatium::Sphere<2, T> unit_sphere{};
                return unit_sphere.log_map((p - c) / r, (q - c) / r) * r;
            }
        }
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
