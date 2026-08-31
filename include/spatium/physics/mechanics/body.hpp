#pragma once

// PointMass<N, T> — minimal body for flat-space mechanics. Holds mass +
// phase state. Rigid bodies (orientation + inertia tensor) are handled
// separately via Lie groups — see manifold_body.hpp/lie_integrator.hpp.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/physics/mechanics/state.hpp>
#  include <spatium/physics/mechanics/units.hpp>
#  include <cstddef>
#endif

SPATIUM_EXPORT namespace spatium::physics::mechanics {

template<std::size_t N, Scalar T = double>
struct PointMass {
    using ScalarType = T;
    static constexpr std::size_t dimension = N;

    T mass{T{1}};
    PhaseState<N, T> state{};

    constexpr PointMass() = default;
    constexpr PointMass(T m, Vec<T, N> position, Vec<T, N> velocity)
        : mass(m), state(position, velocity) {}

    [[nodiscard]] constexpr Vec<T, N> momentum() const {
        return Vec<T, N>{state.velocity * mass};
    }

    [[nodiscard]] constexpr T kinetic_energy() const {
        // KE = ½·m·v²
        return T{0.5} * mass * state.velocity.dot(state.velocity);
    }
};

// Total kinetic energy of a system of point masses.
template<typename Bodies>
constexpr auto total_kinetic_energy(const Bodies& bodies) {
    using T = typename Bodies::value_type::ScalarType;
    T sum{};
    for (const auto& b : bodies) sum += b.kinetic_energy();
    return sum;
}

// Total linear momentum.
template<typename Bodies>
constexpr auto total_momentum(const Bodies& bodies) {
    using B = typename Bodies::value_type;
    using T = typename B::ScalarType;
    constexpr std::size_t N = B::dimension;
    Vec<T, N> p{};
    for (const auto& b : bodies) p = Vec<T, N>{p + b.momentum()};
    return p;
}

} // namespace spatium::physics::mechanics
