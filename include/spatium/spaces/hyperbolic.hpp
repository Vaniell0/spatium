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

    ScalarType distance(const PointType& a, const PointType& b) const {
        using std::acosh;
        // Raised to 1 by a comparison, not by max(): max(NaN, 1) is NaN for
        // double but 1 for a Boost number (its max drops the NaN), which made
        // a NaN point sit at distance 0 -- found by the infinity probe over
        // Real50. `NaN < 1` is false everywhere, so the NaN stays.
        T inner = -minkowski(a, b);
        if (inner < T{1}) inner = T{1};
        return acosh(inner);
    }

    // As Sphere::exp_map: the arc length t |v|_L is signed and linear in t,
    // so exp is smooth through t = 0 and a Dual t carries its derivative v.
    PointType exp_map(const PointType& p, const TangentVector& v, ScalarType t) const {
        using std::abs;
        // p cosh(s) + v t sinhc(s) with s^2 = t^2 |v|_L^2, smooth in v
        // through v = 0 (see spaces/cos_sinc.hpp).
        const T s2 = t * t * abs(minkowski(v, v));
        const auto [c, sinhc] = cos_sinc_of_square(s2, /*hyperbolic=*/true);
        return PointType{p * c + v * (t * sinhc)};
    }

    TangentVector log_map(const PointType& p, const PointType& q) const {
        using std::acosh; using std::sqrt; using std::abs;
        T inner = -minkowski(p, q);       // raised to 1 as in distance(): NaN stays NaN
        if (inner < T{1}) inner = T{1};
        auto d = acosh(inner);
        if (d < epsilon<T>()) return TangentVector{};
        auto proj = q + p * minkowski(p, q);
        auto proj_norm_sq = minkowski(proj, proj);
        if (abs(proj_norm_sq) < epsilon<T>() * epsilon<T>()) return TangentVector{};
        auto proj_norm = sqrt(abs(proj_norm_sq));
        return proj * (d / proj_norm);
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
