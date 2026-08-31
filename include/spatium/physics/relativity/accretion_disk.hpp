#pragma once

// Thin-disk emission-frequency shift for a circular equatorial orbit.
// The disk's radial temperature/brightness PROFILE is an artistic choice
// left to the caller (examples/blackhole_gr_demo.cpp); this file only
// has the physics -- how much a photon's frequency shifts between
// emission at the disk and reception at infinity, combining
// gravitational redshift and orbital Doppler into one formula.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium::physics::relativity {

// Marginally-stable circular orbit (ISCO): r = 6M. Real disk material
// spirals in rapidly once inside this radius rather than orbiting
// steadily, so a thin disk's conventional inner edge stops here, not at
// the horizon.
template<Scalar T>
T schwarzschild_isco_radius(T mass) { return T{6} * mass; }

// Keplerian angular velocity Omega = dphi/dt for a circular equatorial
// orbit at radius r (units G=c=1): Omega = sqrt(M/r^3). Only real for
// r > 3M -- inside the photon sphere no circular orbit exists at all,
// which is why disk rendering must stay outside it.
template<Scalar T>
T keplerian_omega(T mass, T r) {
    using std::sqrt;
    return sqrt(mass / (r * r * r));
}

// u^t of a circular-orbit observer at radius r, from normalizing
// u^mu = u^t*(1, 0, 0, Omega) against the Schwarzschild metric:
// u^t = 1/sqrt(1 - 3M/r). Diverges at the photon sphere for the same
// reason Omega's domain does -- a circular orbit there needs light speed.
template<Scalar T>
T circular_orbit_ut(T mass, T r) {
    using std::sqrt;
    return T{1} / sqrt(T{1} - T{3} * mass / r);
}

// Emitted-to-observed frequency ratio (1+z) for a photon with impact
// parameter b = L/E crossing the disk at radius r, received by a static
// observer at infinity:
//   1+z = u^t_disk * (1 - b*Omega(r))
// Standard thin-disk result (the same formula underlies Luminet 1979's
// original black hole imaging paper and the "Interstellar" visualization
// paper, James et al. 2015) -- a derivation to look up, not open
// research. b>0 (photon angular momentum aligned with the disk's own
// rotation) reduces (1+z) below the pure-gravitational baseline
// (approaching-side Doppler blueshift partially cancels it); b<0
// increases it -- the real asymmetric brightness pattern seen across a
// rotating disk's near and far limbs.
template<Scalar T>
T disk_redshift_factor(T mass, T r, T b) {
    T omega = keplerian_omega(mass, r);
    T ut = circular_orbit_ut(mass, r);
    return ut * (T{1} - b * omega);
}

} // namespace spatium::physics::relativity
