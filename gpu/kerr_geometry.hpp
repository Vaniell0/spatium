#pragma once

// Device-portable transcription of blackhole_gr_demo.cpp's oblate-
// spheroidal tetrad (spheroidal_tetrad()) and the Kerr metric's diagonal
// components (needed at a ray's exit point to convert its 4-velocity
// into local orthonormal components) -- see that file's own comment
// above spheroidal_tetrad() for why the ambient-embedding tangent
// vectors ARE an exact orthonormal tetrad here (Boyer-Lindquist's
// g_r_theta = g_r_phi = g_theta_phi = 0 identically).
//
// Also carries kerr_equatorial_omega/kerr_equatorial_ut (from kerr.hpp)
// as device functions, needed for the disk redshift factor inside the
// render kernel.

#include "disk_physics.hpp"

#ifndef SPATIUM_BUILDING_MODULE
#  include <cmath>
#endif

#if defined(__CUDACC__)
#  define SPATIUM_HOST_DEVICE __host__ __device__
#else
#  define SPATIUM_HOST_DEVICE
#endif

namespace spatium::gpu {

struct SpheroidalTetradD {
    Vec3d e_r, e_theta, e_phi;
};

SPATIUM_HOST_DEVICE inline Vec3d normalize3(Vec3d v) {
    using std::sqrt;
    double n = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    return {v.x / n, v.y / n, v.z / n};
}

SPATIUM_HOST_DEVICE inline SpheroidalTetradD spheroidal_tetrad(double r, double theta, double phi,
                                                                 double a) {
    using std::cos;
    using std::sin;
    using std::sqrt;
    double s = sin(theta), c = cos(theta);
    double cp = cos(phi), sp = sin(phi);
    double rho = sqrt(r * r + a * a);
    Vec3d e_r{(r / rho) * s * cp, (r / rho) * s * sp, c};
    Vec3d e_th{rho * c * cp, rho * c * sp, -r * s};
    Vec3d e_ph{-rho * s * sp, rho * s * cp, 0.0};
    return {normalize3(e_r), normalize3(e_th), normalize3(e_ph)};
}

// Sigma/Delta as in christoffel_closed_form.hpp -- diagonal metric
// components only (g_rr, g_thth, g_phph); g_tt/g_tphi aren't needed at
// the exit point.
SPATIUM_HOST_DEVICE inline void kerr_metric_diag(double mass, double spin, double r, double theta,
                                                   double* g_rr, double* g_thth, double* g_phph) {
    using std::cos;
    using std::sin;
    double s = sin(theta), c = cos(theta);
    double sigma = r * r + spin * spin * c * c;
    double delta = r * r - 2.0 * mass * r + spin * spin;
    *g_rr = sigma / delta;
    *g_thth = sigma;
    *g_phph = (r * r + spin * spin + 2.0 * mass * r * spin * spin * s * s / sigma) * s * s;
}

SPATIUM_HOST_DEVICE inline double kerr_equatorial_omega(double mass, double spin, double r,
                                                          bool prograde) {
    using std::sqrt;
    double sign = prograde ? 1.0 : -1.0;
    return sign * sqrt(mass) / (r * sqrt(r) + sign * spin * sqrt(mass));
}

SPATIUM_HOST_DEVICE inline double kerr_equatorial_ut(double mass, double spin, double r,
                                                       bool prograde) {
    using std::pow;
    using std::sqrt;
    double sign = prograde ? 1.0 : -1.0;
    double sqrt_r = sqrt(r);
    double sqrt_m = sqrt(mass);
    double num = r * sqrt_r + sign * spin * sqrt_m;
    double inner = r * sqrt_r - 3.0 * mass * sqrt_r + sign * 2.0 * spin * sqrt_m;
    return num / (pow(r, 0.75) * sqrt(inner));
}

SPATIUM_HOST_DEVICE inline double kerr_disk_redshift_factor(double mass, double spin, double r,
                                                              double b, bool prograde) {
    double omega = kerr_equatorial_omega(mass, spin, r, prograde);
    double ut = kerr_equatorial_ut(mass, spin, r, prograde);
    return ut * (1.0 - b * omega);
}

}  // namespace spatium::gpu
