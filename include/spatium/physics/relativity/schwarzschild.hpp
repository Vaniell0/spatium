#pragma once

// Schwarzschild metric in standard (t, r, theta, phi) coordinates,
// signature (-+++), units G=c=1:
//   g_tt = -(1 - 2M/r),  g_rr = 1/(1 - 2M/r),
//   g_thth = r^2,        g_phph = r^2 sin^2(theta)
//
// operator() is templated on its OWN scalar parameter S, independent of
// the mass's stored type T -- deliberately not a std::function<...>, the
// exact shape ImplicitSurface's gradient() needed and didn't have,
// forcing it to finite differences (see docs/architecture.md's audit
// notes). This shape is what lets geodesic.hpp's Dual<T>-seeded
// evaluations substitute S=Dual<T> for S=T with zero change here, giving
// exact metric derivatives instead of approximate ones.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <cmath>
#  include <numbers>
#endif

SPATIUM_EXPORT namespace spatium::physics::relativity {

template<Scalar T = double>
struct SchwarzschildMetric {
    T mass;

    template<Scalar S>
    Matrix<S, 4, 4> operator()(const Vec<S, 4>& x) const {
        using std::sin;
        S r = x[1];
        S theta = x[2];
        S f = S{1} - S{2} * S(mass) / r;
        S s = sin(theta);

        Matrix<S, 4, 4> g{};
        g(0, 0) = -f;
        g(1, 1) = S{1} / f;
        g(2, 2) = r * r;
        g(3, 3) = r * r * s * s;
        return g;
    }
};

// Event horizon (r = 2M) and photon sphere (r = 3M, the unstable circular
// photon orbit) -- both plain radii, not approximations.
template<Scalar T>
constexpr T schwarzschild_horizon_radius(T mass) { return T{2} * mass; }

template<Scalar T>
constexpr T schwarzschild_photon_sphere_radius(T mass) { return T{3} * mass; }

// Critical impact parameter b_crit = 3*sqrt(3)*M: photons with b < b_crit
// are captured, b > b_crit escape (in the b -> b_crit limit, orbiting the
// photon sphere arbitrarily many times first). Closed-form, used to
// sanity-check the geodesic integrator against a known exact value.
template<Scalar T>
T schwarzschild_critical_impact_parameter(T mass) {
    using std::sqrt;
    return T{3} * sqrt(T{3}) * mass;
}

} // namespace spatium::physics::relativity
