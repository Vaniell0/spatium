#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/spaces/cos_sinc.hpp>
#  include <cmath>
#endif

// The spaces of constant curvature kappa as one family: the sphere for kappa > 0,
// the plane at kappa = 0, hyperbolic space for kappa < 0, in the kappa-
// stereographic model -- coordinates in R^N, no embedding -- whose formulas
// (Mobius addition and tan_kappa, artan_kappa) are written once for every sign.
//
// kappa is a value of the space, not a parameter of its type, so a `Dual` kappa
// gives the derivative of a distance, an exp or a log with respect to curvature
// -- how an answer changes as the space bends -- and a continuous family of
// spaces is one object, not three. The point kappa = 0 is a removable
// singularity of every formula (tan(sqrt(kappa) t) / sqrt(kappa) is t there); it
// is crossed by two functions analytic in w = kappa * (squared length), the
// `cos_sinc_signed` of cos_sinc.hpp and the `atan_sinc` below, each a series near
// zero and a closed form beyond, so the same code is right at kappa = 1e-30 and
// at 0, on double and on Real50, with no separate case for the flat space.
//
//   conformal factor      lambda_x = 2 / (1 + kappa |x|^2)
//   Mobius addition       x (+) y
//   exp_x(v)              x (+) tan_kappa(lambda_x |v| / 2) v / |v|
//   log_x(y)              (2 / lambda_x) artan_kappa(|-x (+) y|) (-x (+) y) / |-x (+) y|
//   distance              2 artan_kappa(|-x (+) y|)
//
// Held by the closed forms it must reproduce (Euclidean at 0, Sphere at 1 and
// Hyperbolic at -1 through their stereographic maps, the scaling law
// d_kappa(x, y) = d_1(sqrt(kappa) x, sqrt(kappa) y) / sqrt(kappa)) and by
// continuity through kappa = 0 (tests/test_constant_curvature.cpp).
//
// For kappa < 0 the points are the ball |x| < 1 / sqrt(-kappa); a point outside
// gives NaN, not a wrong number. For kappa > 0 the point at infinity is the one
// the model leaves out, and a geodesic through it is where the chart ends.

SPATIUM_EXPORT namespace spatium {

// atan(sqrt(w)) / sqrt(w), and atanh(sqrt(-w)) / sqrt(-w) for w < 0: artan_kappa(r) =
// r * atan_sinc(kappa r^2). Analytic through zero, 1 - w/3 + w^2/5 - ...; the series
// takes over below a threshold that follows the scalar's epsilon.
template<Scalar T>
T atan_sinc(const T& w) {
    using std::atan; using std::atanh; using std::sqrt; using std::abs;
    const double small = 0.5 * std::pow(machine_epsilon<T>(), 1.0 / 8.0);
    if (abs(primal_double(w)) < small) {
        T term{1}, acc{1};
        for (int k = 1; k <= 7; ++k) {
            term = T(-term * w);
            acc = T(acc + term / T(2 * k + 1));
        }
        return acc;
    }
    if (primal_double(w) > 0.0) {
        const T r = sqrt(w);
        return T(atan(r) / r);
    }
    const T r = sqrt(T(-w));
    return T(atanh(r) / r);
}

template<std::size_t N, Scalar T = double>
struct ConstantCurvature {
    using ScalarType    = T;
    using PointType     = Vec<T, N>;
    using TangentVector = Vec<T, N>;

    static constexpr std::size_t dimension = N;
    static constexpr bool is_complete = true;

    T kappa{};

    constexpr ConstantCurvature() = default;
    explicit constexpr ConstantCurvature(T k) : kappa(k) {}

    // R^N for kappa >= 0; the open ball of radius 1/sqrt(-kappa) for kappa < 0.
    bool contains(const PointType& x) const {
        return !(kappa < T{0}) || T(kappa * x.norm_squared()) > T{-1};
    }

    T conformal_factor(const PointType& x) const { return T(T{2} / (T{1} + kappa * x.norm_squared())); }

    // The metric at x on tangent vectors in coordinates: lambda_x^2 <u, v>.
    ScalarType metric_at(const PointType& x, const TangentVector& u, const TangentVector& v) const {
        const T lam = conformal_factor(x);
        return T(lam * lam * u.dot(v));
    }

    // Gyrovector (Mobius) addition: the translation of the model.
    PointType mobius_add(const PointType& x, const PointType& y) const {
        const T xy = x.dot(y), x2 = x.norm_squared(), y2 = y.norm_squared();
        const T a = T(T{1} - T{2} * kappa * xy - kappa * y2);
        const T b = T(T{1} + kappa * x2);
        const T den = T(T{1} - T{2} * kappa * xy + kappa * kappa * x2 * y2);
        return PointType{(x * a + y * b) / den};
    }

    ScalarType distance(const PointType& x, const PointType& y) const {
        using std::sqrt;
        const PointType u = mobius_add(-x, y);
        const T r2 = u.norm_squared();
        return T(T{2} * sqrt(r2) * atan_sinc(T(kappa * r2)));
    }

    // exp_x(t v): the step is smooth in t and v through zero, as Sphere's.
    PointType exp_map(const PointType& x, const TangentVector& v, ScalarType t) const {
        const PointType u = v * t;
        const T lam = conformal_factor(x);
        const T w = T(kappa * lam * lam * u.norm_squared() / T{4});     // kappa s^2, s = lambda |u| / 2
        const auto [c, sinc] = cos_sinc_signed(w);
        const T tanc = T(sinc / c);                                      // tan_kappa(s) / s
        return mobius_add(x, u * T(lam / T{2} * tanc));
    }

    TangentVector log_map(const PointType& x, const PointType& y) const {
        const PointType u = mobius_add(-x, y);
        const T lam = conformal_factor(x);
        return TangentVector{u * T(T{2} / lam * atan_sinc(T(kappa * u.norm_squared())))};
    }

    static constexpr PointType origin() { return PointType{}; }
};

static_assert(RiemannianManifold<ConstantCurvature<2>>);
static_assert(Complete<ConstantCurvature<3>>);

} // namespace spatium
