#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/spaces/cos_sinc.hpp>
#  include <algorithm>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium {

// N-dimensional hyperbolic space in the hyperboloid model.
// Points on upper sheet: -x0^2 + x1^2 + ... + xN^2 = -1, x0 > 0

template<std::size_t N, Scalar T = double>
struct Hyperbolic {
    using ScalarType    = T;
    using PointType     = Vec<T, N + 1>;
    using TangentVector = Vec<T, N + 1>;

    static constexpr std::size_t dimension = N;
    static constexpr bool is_complete = true;

    static constexpr T minkowski(const PointType& a, const PointType& b) {
        T sum = -a[0] * b[0];
        for (std::size_t i = 1; i <= N; ++i)
            sum += a[i] * b[i];
        return sum;
    }

    bool contains(const PointType& p) const {
        using std::abs;
        return abs(minkowski(p, p) + T{1}) < epsilon<T>() && p[0] > T{0};
    }

    // cosh d - 1 with no cancellation. Read off -<a,b> - 1 a point at distance 0
    // from itself, once exp has carried it a few units out, is 1 plus rounding
    // (acosh of that is 1e-8, not 0). On the sheet a0 b0 = sqrt((1+|a|^2)(1+|b|^2)),
    // and (a0 b0)^2 - (a.b + 1)^2 = |a - b|^2 + sum_{i<j} (a_i b_j - a_j b_i)^2 by
    // Lagrange's identity, so cosh d - 1 is that sum of squares over a0 b0 + a.b + 1,
    // a denominator of at least 1. Equal points give exactly 0.
    // A NaN coordinate stays NaN: every sum below carries it, and nothing is
    // raised by a comparison or max() (which Boost drops a NaN through).
    static T cosh_minus_one(const PointType& a, const PointType& b) {
        T diff{0}, wedge{0}, dot{0};
        for (std::size_t i = 1; i <= N; ++i) {
            const T d = a[i] - b[i];
            diff += d * d;
            dot += a[i] * b[i];
            for (std::size_t j = i + 1; j <= N; ++j) {
                const T w = a[i] * b[j] - a[j] * b[i];
                wedge += w * w;
            }
        }
        return (diff + wedge) / (a[0] * b[0] + dot + T{1});
    }

    // d = 2 asinh sqrt((cosh d - 1)/2) keeps the digits of a short distance as well.
    ScalarType distance(const PointType& a, const PointType& b) const {
        using std::asinh; using std::sqrt;
        return T{2} * asinh(sqrt(cosh_minus_one(a, b) / T{2}));
    }

    // As Sphere::exp_map: the arc length t |v|_L is signed and linear in t,
    // so exp is smooth through t = 0 and a Dual t carries its derivative v.
    // |v|_L^2 of a tangent at p with no cancellation. Far from the origin a
    // tangent's coordinates grow like e^d while its length does not, and
    // -v0^2 + |v_s|^2 loses e^{2d} of the arithmetic's digits (2.6e-4 of a
    // unit step from a point 8 away, all of them 12 away). With v0 = v_s.p_s / p0
    // the same number is (|v_s|^2 + |v_s|^2 |p_s|^2 - (v_s.p_s)^2) / p0^2, and
    // Lagrange's identity writes the difference as a sum of squares.
    static constexpr T tangent_norm_sq(const PointType& p, const TangentVector& v) {
        T spatial{0}, wedge{0};
        for (std::size_t i = 1; i <= N; ++i) {
            spatial += v[i] * v[i];
            for (std::size_t j = i + 1; j <= N; ++j) {
                const T w = v[i] * p[j] - v[j] * p[i];
                wedge += w * w;
            }
        }
        return (spatial + wedge) / (p[0] * p[0]);
    }

    // As Sphere::exp_map: the arc length t |v|_L is signed and linear in t,
    // so exp is smooth through t = 0 and a Dual t carries its derivative v.
    PointType exp_map(const PointType& p, const TangentVector& v, ScalarType t) const {
        // p cosh(s) + v t sinhc(s) with s^2 = t^2 |v|_L^2, smooth in v
        // through v = 0 (see spaces/cos_sinc.hpp).
        const T s2 = t * t * tangent_norm_sq(p, v);
        const auto [c, sinhc] = cos_sinc_of_square(s2, /*hyperbolic=*/true);
        // The time coordinate put back on the sheet as project() does, from the
        // spatial ones with nothing cancelling. A tangent that is tangent only to
        // rounding otherwise moves each point a little off the sheet, and an
        // iteration (the Frechet mean) carries it there step after step: 1e-6
        // off after fifty steps at H^8, where log, which reads sinh d off the
        // inner product, no longer agrees with the point.
        PointType r{p * c + v * (t * sinhc)};
        using std::sqrt;
        T spatial{0};
        for (std::size_t i = 1; i <= N; ++i) spatial += r[i] * r[i];
        r[0] = sqrt(T{1} + spatial);
        return r;
    }

    // proj = q + p <p,q> = (q - p) - p c with c = cosh d - 1 as distance() reads it, of
    // length sinh d = sqrt(c (c + 2)): nothing cancels in either, so a point's log at
    // itself is 0 and a short way keeps its digits however far out the two points are.
    // Read off acosh(-<p,q>), -<q,q> was 1 plus rounding 9 away from the origin, log_q q
    // a vector 8e-5 long and off the tangent space, and log_q of a point 1e-4 away 3.5e-2
    // wrong at 8 away (9 at 12).
    TangentVector log_map(const PointType& p, const PointType& q) const {
        using std::asinh; using std::sqrt;
        const T c = cosh_minus_one(p, q);
        if (!(c > T{0})) return c == T{0} ? TangentVector{} : TangentVector{q * c};   // NaN stays NaN
        const T d = T{2} * asinh(sqrt(c / T{2}));
        const TangentVector proj{(q - p) - p * c};
        return proj * (d / sqrt(c * (c + T{2})));
    }

    constexpr ScalarType metric_at(const PointType&,
                                   const TangentVector& u,
                                   const TangentVector& v) const {
        return minkowski(u, v);
    }

    PointType project(const PointType& p) const {
        using std::sqrt;
        T spatial_sq{0};
        for (std::size_t i = 1; i <= N; ++i)
            spatial_sq += p[i] * p[i];
        PointType result = p;
        result[0] = sqrt(T{1} + spatial_sq);
        return result;
    }

    TangentVector normal(const PointType& p) const { return p; }

    static constexpr PointType origin() {
        PointType p{};
        p[0] = T{1};
        return p;
    }
};

static_assert(RiemannianManifold<Hyperbolic<2>>);
static_assert(Surface<Hyperbolic<2>>);
static_assert(Complete<Hyperbolic<2>>);
static_assert(RiemannianManifold<Hyperbolic<1>>);
static_assert(Surface<Hyperbolic<3>>);

using H1 = Hyperbolic<1>;
using H2 = Hyperbolic<2>;
using H3 = Hyperbolic<3>;

} // namespace spatium
