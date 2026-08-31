#pragma once

// Phase-space state for N-dimensional Euclidean mechanics.
//
// `PhaseState<N, T>` bundles configuration (position) with tangent (velocity).
// A future generalization to `PhaseState<Manifold>` via exp/log maps is
// possible; for now everything lives in flat R^N.
//
// Linear combinations (`a + b`, `s * a`) are defined to feed integrators
// (Euler, RK4) directly — they treat the state as a point in 2N-dim phase space.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <cstddef>
#endif

SPATIUM_EXPORT namespace spatium::physics::mechanics {

template<std::size_t N, Scalar T = double>
struct PhaseState {
    using ScalarType = T;
    static constexpr std::size_t dimension = N;

    Vec<T, N> position{};
    Vec<T, N> velocity{};

    constexpr PhaseState() = default;
    constexpr PhaseState(Vec<T, N> p, Vec<T, N> v) : position(p), velocity(v) {}

    constexpr PhaseState operator+(const PhaseState& rhs) const {
        return {Vec<T, N>{position + rhs.position}, Vec<T, N>{velocity + rhs.velocity}};
    }
    constexpr PhaseState operator-(const PhaseState& rhs) const {
        return {Vec<T, N>{position - rhs.position}, Vec<T, N>{velocity - rhs.velocity}};
    }
    constexpr PhaseState operator*(T s) const {
        return {Vec<T, N>{position * s}, Vec<T, N>{velocity * s}};
    }
    friend constexpr PhaseState operator*(T s, const PhaseState& st) { return st * s; }
};

// `Derivative<N>` — what an integrator gets from a force model: dq/dt = v,
// dv/dt = a. We just reuse PhaseState as a tagged tangent — `position` slot
// holds dq/dt, `velocity` slot holds dv/dt. Strict typing isn't worth the
// boilerplate at this scale; the convention is documented and tested.
template<std::size_t N, Scalar T = double>
using Derivative = PhaseState<N, T>;

} // namespace spatium::physics::mechanics
