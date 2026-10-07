#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/concepts.hpp>
#  include <spatium/core/concepts.hpp>
#  include <cmath>
#endif

// A Lie group as a space, so that an algorithm written against exp, log and a metric
// -- `frechet_mean`, `geodesic`, `midpoint`, `verify_exp_log` -- runs on rotations and
// rigid motions with no copy of it for groups.
//
// Points are group elements, tangent vectors are elements of the Lie algebra (the
// tangent at the identity, carried to p by left translation), and
//
//   exp_p(t v) = p exp(t v)            log_p(q) = log(p^-1 q)
//   <u, v>_p   = u . v                  (the same at every point: left-invariant)
//
// For SO(3) this is exactly the Riemannian geometry of the bi-invariant metric --
// its geodesics are the one-parameter subgroups, the distance is the rotation angle --
// so the Frechet mean is the rotation average of SLAM and sensor fusion, and it is
// equivariant under left and right translation. SE(3) has no bi-invariant metric: there
// exp and log are the group's own (the Cartan-Schouten connection), which invert each
// other and give the mean everyone computes for poses, but are not the geodesics of the
// left-invariant metric u . v; the distance below is that of the algebra's norm, a
// consistent choice rather than a Levi-Civita one.

SPATIUM_EXPORT namespace spatium {

template<class G>
    requires LieGroup<G>
struct LieGroupManifold {
    using ScalarType    = typename G::ScalarType;
    using PointType     = typename G::ElementType;
    using TangentVector = typename G::AlgebraType;

    static constexpr std::size_t dimension = TangentVector::size;
    static constexpr bool is_complete = true;

    G group{};

    constexpr bool contains(const PointType&) const { return true; }

    PointType exp_map(const PointType& p, const TangentVector& v, ScalarType t) const {
        return group.compose(p, group.exp(TangentVector{v * t}));
    }

    TangentVector log_map(const PointType& p, const PointType& q) const {
        return TangentVector{group.log(group.compose(group.inverse(p), q))};
    }

    ScalarType metric_at(const PointType&, const TangentVector& u, const TangentVector& v) const {
        return u.dot(v);
    }

    ScalarType distance(const PointType& p, const PointType& q) const {
        using std::sqrt;
        return sqrt(log_map(p, q).norm_squared());
    }
};

} // namespace spatium
