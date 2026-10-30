#pragma once

// Closed-form Christoffel symbols for Schwarzschild and Kerr (Boyer-
// Lindquist), derived once by sympy from the exact metric forms in
// schwarzschild.hpp / kerr.hpp, checked symbolically (Kerr at a = 0 reduces
// term by term to Schwarzschild), and written out as straight-line
// arithmetic.
//
// They were the CUDA kernels' Christoffels; the kernels are gone and the
// device path takes its symbols from the metric's IR with exact
// derivatives. What remains is their use as a witness: an independent
// derivation that tests/test_relativity.cpp holds the generic Dual<T>
// christoffel() against at a spread of (r, theta, a). Do not edit the
// expressions by hand.
//
// Index convention throughout: 0=t, 1=r, 2=theta, 3=phi, matching
// geodesic.hpp's Vec<T,8> state layout and christoffel()'s
// Gamma[lambda](mu,nu).

#include <cmath>

namespace spatium::closed_form {

// f = 1 - 2M/r. Gamma[lambda][mu][nu], both (mu,nu) orderings filled in
// (the tensor is symmetric in its lower two indices).
template<typename T>
void schwarzschild_christoffel_closed_form(T mass, T r, T theta,
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
void kerr_christoffel_closed_form(T mass, T spin, T r, T theta,
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

}  // namespace spatium::closed_form
