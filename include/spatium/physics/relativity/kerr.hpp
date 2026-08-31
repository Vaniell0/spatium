#pragma once

// Kerr metric in Boyer-Lindquist coordinates (t, r, theta, phi),
// signature (-+++), units G=c=1, spin parameter a (0 <= |a| <= mass for
// a physical, non-naked-singularity hole):
//   Sigma = r^2 + a^2 cos^2(theta),  Delta = r^2 - 2Mr + a^2
//   g_tt    = -(1 - 2Mr/Sigma)
//   g_tphi  = g_phit = -2Mr a sin^2(theta) / Sigma
//   g_rr    = Sigma / Delta
//   g_thth  = Sigma
//   g_phph  = (r^2 + a^2 + 2Mr a^2 sin^2(theta)/Sigma) sin^2(theta)
// a=0 reduces exactly to SchwarzschildMetric (verified in
// tests/test_relativity.cpp, not just asserted). Same templated-callable
// shape as SchwarzschildMetric (see that file's header comment for why
// that shape matters), so geodesic.hpp's Dual<T>-seeded
// metric_derivatives()/christoffel() need zero changes -- this is the
// "future drop-in" docs/architecture.md and schwarzschild.hpp's own
// header comment promised: the metric inverse there already uses the
// general invert(), not a diagonal-only shortcut, specifically because
// g_tphi makes this metric non-diagonal.
//
// The BPT (Bardeen, Press & Teukolsky 1972) closed-form radii/orbit
// formulas below all take a `bool prograde` (co-rotating with the hole's
// spin, vs. counter-rotating) -- a single consistent convention across
// all four functions, even though the underlying textbook formulas don't
// agree with each other on which raw +-1 sign means "prograde" (the
// photon-orbit/ISCO formulas and the Omega/u^t formulas use opposite
// raw-sign conventions in the literature); each function flips its own
// internal sign so callers never have to know that. Every one of them
// collapses to its Schwarzschild value at a=0 regardless of prograde
// (cross-checked in tests/test_relativity.cpp), and at the extremal
// a=mass limit reproduces the textbook results (prograde photon
// orbit/ISCO both -> mass, i.e. the horizon; retrograde photon orbit ->
// 4*mass, retrograde ISCO -> 9*mass) -- also cross-checked.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/matrix.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium::physics::relativity {

template<Scalar T = double>
struct KerrMetric {
    T mass;
    T spin;

    template<Scalar S>
    Matrix<S, 4, 4> operator()(const Vec<S, 4>& x) const {
        using std::cos;
        using std::sin;
        S r = x[1];
        S theta = x[2];
        S a = S(spin);
        S c = cos(theta);
        S s = sin(theta);
        S sigma = r * r + a * a * c * c;
        S delta = r * r - S{2} * S(mass) * r + a * a;

        Matrix<S, 4, 4> g{};
        g(0, 0) = -(S{1} - S{2} * S(mass) * r / sigma);
        g(0, 3) = g(3, 0) = -S{2} * S(mass) * r * a * s * s / sigma;
        g(1, 1) = sigma / delta;
        g(2, 2) = sigma;
        g(3, 3) = (r * r + a * a + S{2} * S(mass) * r * a * a * s * s / sigma) * s * s;
        return g;
    }
};

// Outer event horizon r+ = M + sqrt(M^2 - a^2).
template<Scalar T>
T kerr_outer_horizon_radius(T mass, T spin) {
    using std::sqrt;
    return mass + sqrt(mass * mass - spin * spin);
}

// Equatorial circular photon orbit radius. a=0 gives 3*mass regardless of
// prograde -- the Schwarzschild photon sphere.
template<Scalar T>
T kerr_photon_orbit_radius(T mass, T spin, bool prograde) {
    using std::acos;
    using std::cos;
    T sign = prograde ? T{-1} : T{1};
    return T{2} * mass * (T{1} + cos(T{2} / T{3} * acos(sign * spin / mass)));
}

// Marginally-stable circular orbit radius (ISCO). a=0 gives 6*mass
// regardless of prograde -- the Schwarzschild ISCO.
template<Scalar T>
T kerr_isco_radius(T mass, T spin, bool prograde) {
    using std::cbrt;
    using std::sqrt;
    T sign = prograde ? T{-1} : T{1};
    T a_star = spin / mass;
    T z1 = T{1} + cbrt(T{1} - a_star * a_star) * (cbrt(T{1} + a_star) + cbrt(T{1} - a_star));
    T z2 = sqrt(T{3} * a_star * a_star + z1 * z1);
    return mass * (T{3} + z2 + sign * sqrt((T{3} - z1) * (T{3} + z1 + T{2} * z2)));
}

// Keplerian angular velocity for an equatorial circular orbit at radius
// r. a=0 recovers +-sqrt(mass/r^3) (positive for prograde by convention,
// since Schwarzschild itself has no preferred rotation direction).
template<Scalar T>
T kerr_equatorial_omega(T mass, T spin, T r, bool prograde) {
    using std::sqrt;
    T sign = prograde ? T{1} : T{-1};
    return sign * sqrt(mass) / (r * sqrt(r) + sign * spin * sqrt(mass));
}

// u^t for an equatorial circular orbit at radius r -- Bardeen, Press &
// Teukolsky 1972.
template<Scalar T>
T kerr_equatorial_ut(T mass, T spin, T r, bool prograde) {
    using std::pow;
    using std::sqrt;
    T sign = prograde ? T{1} : T{-1};
    T sqrt_r = sqrt(r);
    T sqrt_m = sqrt(mass);
    T num = r * sqrt_r + sign * spin * sqrt_m;
    T inner = r * sqrt_r - T{3} * mass * sqrt_r + sign * T{2} * spin * sqrt_m;
    return num / (pow(r, T{0.75}) * sqrt(inner));
}

// Emitted-to-observed frequency ratio for a photon with impact parameter
// b crossing the equatorial plane at radius r, from a circular-orbit
// emitter there: 1+z = u^t*(1 - b*Omega) -- the same general Killing-
// vector argument as schwarzschild.hpp's disk_redshift_factor(), which
// was never actually specific to a diagonal metric, only fed diagonal
// Omega/u^t.
template<Scalar T>
T kerr_disk_redshift_factor(T mass, T spin, T r, T b, bool prograde) {
    T omega = kerr_equatorial_omega(mass, spin, r, prograde);
    T ut = kerr_equatorial_ut(mass, spin, r, prograde);
    return ut * (T{1} - b * omega);
}

} // namespace spatium::physics::relativity
