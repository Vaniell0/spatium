#pragma once

// IPC (Incremental Potential Contact) barrier potential.
//
// Li, Ferguson, Schneider, Langlois, Zorin, Panozzo, Jiang, Kaufman 2020.
// "Incremental Potential Contact: Intersection- and Inversion-free,
// Large-deformation Dynamics." SIGGRAPH 2020.
//
// The barrier energy
//
//     B(d, d̂) = -(d − d̂)² · log(d / d̂)         for 0 < d < d̂
//             = 0                                for d ≥ d̂
//
// blows up as d → 0⁺ and is C² across the threshold d̂. Integrated as a
// soft potential into the system Hamiltonian / Lagrangian, it provides
// **provably no-penetration** dynamics — the optimisation step refuses
// to cross d = 0 because the energy diverges. With a backtracking line
// search filter (Li et al. §5), the discrete trajectory stays
// intersection-free for all time.
//
// This header ships the energy and its derivatives as pure math; the full
// IPC toolkit (collision detection, line-search filter, friction)
// integrates later via the external `ipc-toolkit` dependency
// (https://github.com/ipc-sim/ipc-toolkit). For now we expose
// only the kernel so it can be combined with the variational framework
// (`mechanics/variational.hpp`) and the analytical-distance narrow-phase
// queries from `geometry/ray_*.hpp`.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium::physics::mechanics {

// ── IPC barrier energy B(d, d̂) ─────────────────────────────────
// Returns 0 if d ≥ d̂, otherwise the smooth log-barrier value. The energy
// scale grows like (log d̂/d) near d → 0 and vanishes smoothly at d = d̂
// with B(d̂) = B'(d̂) = 0 (C¹ across the threshold).
//
// Conventions: distance `d` should be ≥ 0 (signed/unsigned distance to
// nearest surface point). `d_hat` is the activation threshold > 0.
template<Scalar T>
constexpr T ipc_barrier(T d, T d_hat) {
    if (d >= d_hat) return T{0};
    if (d <= T{0}) {
        using std::numeric_limits;
        return numeric_limits<T>::infinity();
    }
    T diff = d - d_hat;
    using std::log;
    return -diff * diff * log(d / d_hat);
}

// First derivative: dB/dd. Negative when d → 0 (force points away from
// contact), zero at d = d̂, smooth in between.
template<Scalar T>
T ipc_barrier_grad(T d, T d_hat) {
    if (d >= d_hat) return T{0};
    if (d <= T{0}) {
        using std::numeric_limits;
        return -numeric_limits<T>::infinity();
    }
    T diff = d - d_hat;
    using std::log;
    return -T{2} * diff * log(d / d_hat) - diff * diff / d;
}

// Second derivative: d²B/dd². Used in IPC's Newton solver to produce the
// stiffness contribution of contact. Always positive in the active region
// (convex barrier).
template<Scalar T>
T ipc_barrier_hessian(T d, T d_hat) {
    if (d >= d_hat) return T{0};
    if (d <= T{0}) {
        using std::numeric_limits;
        return numeric_limits<T>::infinity();
    }
    T diff = d - d_hat;
    using std::log;
    T term1 = -T{2} * log(d / d_hat);
    T term2 = -T{4} * diff / d;          // combines two -2·diff/d cross terms
    T term3 =  diff * diff / (d * d);
    return term1 + term2 + term3;
}

// Stiffness scale κ — typical IPC choice that gives unit-magnitude
// contact forces near `d_hat / 2`. Multiply the barrier by κ to make
// the contact response comparable to other forces in the system.
template<Scalar T>
constexpr T ipc_default_stiffness(T d_hat) {
    // Empirical: κ ≈ 1 / B(d̂/2, d̂). Cf. Li et al. §5 default.
    using std::log;
    T half = d_hat * T{0.5};
    T diff = half - d_hat;
    T mag  = -diff * diff * log(half / d_hat);
    return (mag > T{0}) ? T{1} / mag : T{1};
}

// Convenience: full barrier potential including default stiffness.
template<Scalar T>
constexpr T ipc_potential(T d, T d_hat) {
    return ipc_default_stiffness(d_hat) * ipc_barrier(d, d_hat);
}

} // namespace spatium::physics::mechanics
