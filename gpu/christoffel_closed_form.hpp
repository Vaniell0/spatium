#pragma once

// Closed-form Christoffel symbols for Schwarzschild and Kerr (Boyer-
// Lindquist), for the GPU geodesic kernel specifically -- deliberately
// NOT built on physics/relativity/geodesic.hpp's generic Dual<T>-seeded
// metric_derivatives() + general 4x4 invert(). That machinery is the
// CPU engine's own considered trade-off (pay indirection/an extra
// matrix inverse per step to support any future metric with zero new
// derivation work); a SIMT kernel that will only ever run these two
// published metrics pays for that indirection/branching without ever
// using the generality, so here the expressions are hand/symbolically
// derived once and dropped in as straight-line arithmetic instead.
//
// Derived via sympy directly from the exact metric forms in
// schwarzschild.hpp/kerr.hpp (not transcribed from a textbook table by
// hand -- less reliable than an independent symbolic computation, the
// same reasoning behind this project's other cross-checks), then
// verified two ways before ever being used here:
//   1. symbolically, Kerr at a=0 reduces term-by-term to the
//      Schwarzschild expressions (exact, sympy-checked);
//   2. numerically, both functions below are cross-checked against
//      geodesic.hpp's own christoffel() (the Dual<T>-exact ground
//      truth) at a spread of (r, theta, a) in
//      tests/test_relativity.cpp's "closed-form Christoffel" cases.
// Do not hand-edit these expressions without re-running that
// derivation -- see gpu/derive_christoffel.py.
//
// No Scalar concept constraint and no spatium/ includes: this header
// is meant to compile standalone under nvcc for __device__ code, so it
// stays free of anything the CPU-only generic engine depends on.
// Index convention throughout: 0=t, 1=r, 2=theta, 3=phi, matching
// geodesic.hpp's Vec<T,8> state layout and christoffel()'s
// Gamma[lambda](mu,nu).

#ifndef SPATIUM_BUILDING_MODULE
#  include <cmath>
#endif

#if defined(__CUDACC__)
#  define SPATIUM_HOST_DEVICE __host__ __device__
#else
#  define SPATIUM_HOST_DEVICE
#endif

namespace spatium::gpu {

// f = 1 - 2M/r. Gamma[lambda][mu][nu], both (mu,nu) orderings filled in
// (the tensor is symmetric in its lower two indices).
template<typename T>
SPATIUM_HOST_DEVICE void schwarzschild_christoffel_closed_form(T mass, T r, T theta,
                                                                T Gamma[4][4][4]) {
    using std::cos;
    using std::sin;
    for (int l = 0; l < 4; ++l)
        for (int m = 0; m < 4; ++m)
            for (int n = 0; n < 4; ++n) Gamma[l][m][n] = T{0};

    T M = mass;
    T s = sin(theta), c = cos(theta);
    T f = T{1} - T{2} * M / r;
    T r2 = r * r;

    T g_t_tr = M / (r2 * f);
    Gamma[0][0][1] = Gamma[0][1][0] = g_t_tr;

    T g_r_tt = M * f / r2;
    Gamma[1][0][0] = g_r_tt;

    T g_r_rr = -M / (r2 * f);
    Gamma[1][1][1] = g_r_rr;

    T g_r_thth = -r * f;
    Gamma[1][2][2] = g_r_thth;

    T g_r_phph = -r * f * s * s;
    Gamma[1][3][3] = g_r_phph;

    T g_th_rth = T{1} / r;
    Gamma[2][1][2] = Gamma[2][2][1] = g_th_rth;

    T g_th_phph = -s * c;
    Gamma[2][3][3] = g_th_phph;

    T g_ph_rph = T{1} / r;
    Gamma[3][1][3] = Gamma[3][3][1] = g_ph_rph;

    T g_ph_thph = c / s;
    Gamma[3][2][3] = Gamma[3][3][2] = g_ph_thph;
}

// Sigma = r^2 + a^2*cos^2(theta), Delta = r^2 - 2Mr + a^2.
template<typename T>
SPATIUM_HOST_DEVICE void kerr_christoffel_closed_form(T mass, T spin, T r, T theta,
                                                        T Gamma[4][4][4]) {
    using std::cos;
    using std::sin;
    for (int l = 0; l < 4; ++l)
        for (int m = 0; m < 4; ++m)
            for (int n = 0; n < 4; ++n) Gamma[l][m][n] = T{0};

    T M = mass, a = spin;
    T s = sin(theta), c = cos(theta);
    T s2 = s * s, c2 = c * c;
    T r2 = r * r, r3 = r2 * r, r4 = r2 * r2, r5 = r4 * r;
    T a2 = a * a, a4 = a2 * a2;
    T sigma = r2 + a2 * c2;
    T delta = r2 - T{2} * M * r + a2;
    T sigma2 = sigma * sigma, sigma3 = sigma2 * sigma;
    T sin2th = T{2} * s * c;  // sin(2*theta)

    T g_t_tr = M * (a4 * s2 - a4 + a2 * r2 * s2 + r4) / (delta * sigma2);
    Gamma[0][0][1] = Gamma[0][1][0] = g_t_tr;

    T g_t_tth = -M * a2 * r * sin2th / sigma2;
    Gamma[0][0][2] = Gamma[0][2][0] = g_t_tth;

    T g_t_rph = M * a * (a4 * c2 - a2 * r2 * c2 - a2 * r2 - T{3} * r4) * s2 / (delta * sigma2);
    Gamma[0][1][3] = Gamma[0][3][1] = g_t_rph;

    T g_t_thph = T{2} * M * a2 * a * r * s2 * s * c / sigma2;
    Gamma[0][2][3] = Gamma[0][3][2] = g_t_thph;

    T g_r_tt = M * delta * (r2 - a2 * c2) / sigma3;
    Gamma[1][0][0] = g_r_tt;

    T g_r_tph = M * a * delta * (a2 * c2 - r2) * s2 / sigma3;
    Gamma[1][0][3] = Gamma[1][3][0] = g_r_tph;

    T g_r_rr = M / delta + r / sigma - r / delta;
    Gamma[1][1][1] = g_r_rr;

    T g_r_rth = -a2 * sin2th / (T{2} * sigma);
    Gamma[1][1][2] = Gamma[1][2][1] = g_r_rth;

    T g_r_thth = -delta * r / sigma;
    Gamma[1][2][2] = g_r_thth;

    T g_r_phph = delta * (T{2} * M * a2 * r2 * s2 * s2 - M * a2 * sigma * s2 * s2 - r * sigma2 * s2) /
                 sigma3;
    Gamma[1][3][3] = g_r_phph;

    T g_th_tt = -M * a2 * r * sin2th / sigma3;
    Gamma[2][0][0] = g_th_tt;

    T g_th_tph = M * a * r * (a2 + r2) * sin2th / sigma3;
    Gamma[2][0][3] = Gamma[2][3][0] = g_th_tph;

    T g_th_rr = a2 * sin2th / (T{2} * delta * sigma);
    Gamma[2][1][1] = g_th_rr;

    T g_th_rth = r / sigma;
    Gamma[2][1][2] = Gamma[2][2][1] = g_th_rth;

    T g_th_thth = -a2 * sin2th / (T{2} * sigma);
    Gamma[2][2][2] = g_th_thth;

    T g_th_phph = -(T{2} * M * a2 * r * (a2 + r2) * s2 + sigma * (T{2} * M * a2 * r * s2 + sigma * (a2 + r2))) *
                  s * c / sigma3;
    Gamma[2][3][3] = g_th_phph;

    T g_ph_tr = -M * a * (a2 * c2 - r2) / (delta * sigma2);
    Gamma[3][0][1] = Gamma[3][1][0] = g_ph_tr;

    T g_ph_tth = -T{2} * M * a * r * c / (sigma2 * s);
    Gamma[3][0][2] = Gamma[3][2][0] = g_ph_tth;

    T g_ph_rph = (-M * a4 * c2 * c2 + M * a4 * c2 - M * a2 * r2 * c2 - M * a2 * r2 - T{2} * M * r4 +
                  a4 * r * c2 * c2 + T{2} * a2 * r3 * c2 + r5) /
                 (delta * sigma2);
    Gamma[3][1][3] = Gamma[3][3][1] = g_ph_rph;

    T g_ph_thph = (T{2} * M * a2 * r * s2 + a4 * s2 * s2 - T{2} * a4 * s2 + a4 - T{2} * a2 * r2 * s2 +
                   T{2} * a2 * r2 + r4) *
                  c / (sigma2 * s);
    Gamma[3][2][3] = Gamma[3][3][2] = g_ph_thph;
}

}  // namespace spatium::gpu
