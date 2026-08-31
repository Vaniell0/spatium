#pragma once

// Rigid-body-integrator method dispatch over PointMass<2,double>, the
// domain identified as RSC-ready-with-no-new-algorithm-work: five real,
// independently-implemented steppers already exist in
// spatium/physics/mechanics/integrator.hpp (euler_step,
// semi_implicit_euler_step, verlet_step, rk4_step, yoshida4_step),
// exactly the "several interchangeable candidates for one state-stepping
// problem" shape ode_ops.hpp/ode_task.hpp already proved out for the
// generic-ODE case -- this reuses that shape, not a new one.
//
// Three test families, all reduced to mass=1, reference length/speed=1 so
// a single characteristic frequency feature (sqrt(rate)) applies
// uniformly across all three -- unlike ode_ops.hpp's CircularOrbit-only
// sqrt() case:
//   - UniformGravity: F = -m*g*y_hat, closed form is exactly quadratic in
//     t. Only Euler and semi-implicit Euler have nonzero error here
//     (missing/approximating the O(dt^2) term); Verlet/RK4/Yoshida4 are
//     exact to floating-point precision regardless of dt, since they all
//     integrate a degree<=2 polynomial acceleration exactly. A real, if
//     narrower, discriminator -- not degenerate.
//   - Spring: isotropic 2D Hookean spring through the origin (rest_length
//     0), reduces to two independent 1D SHM axes. w = sqrt(k/m) = sqrt(k).
//   - PointGravity: 2D circular orbit, same shape as ode_ops.hpp's
//     CircularOrbit family. w = sqrt(GM/r0^3) = sqrt(GM) at r0=1.
//
// Candidates ordered cheapest-to-priciest by force evaluations per step
// (euler/semi-implicit: 1, verlet: 2, rk4: 4, yoshida4: 6 -- three Verlet
// sub-steps), matching ComparisonTaskGenerator's cheapest-first contract;
// yoshida4 last, playing the same "closest-to-exact, compared-against"
// role rk4 plays in ode_ops.hpp.

#include <registry.hpp>
#include <spatium/physics/mechanics/body.hpp>
#include <spatium/physics/mechanics/force.hpp>
#include <spatium/physics/mechanics/integrator.hpp>
#include <cmath>
#include <vector>

namespace rsc {

enum class IntegratorFamily { UniformGravity = 0, Spring = 1, PointGravity = 2 };

namespace integrator_detail {

using spatium::physics::mechanics::PointMass;
using spatium::physics::mechanics::UniformGravity;
using spatium::physics::mechanics::Spring;
using spatium::physics::mechanics::PointGravity;

// All three families: mass = 1, reference length/speed = 1 -- picked so
// every family's characteristic frequency reduces to the same sqrt(rate)
// form (see characteristic_frequency below).

inline UniformGravity<2, double> uniform_gravity_force(double rate) {
    return UniformGravity<2, double>{spatium::Vec<double, 2>{0.0, -rate}};
}
inline PointMass<2, double> uniform_gravity_y0() {
    return PointMass<2, double>{1.0, {0.0, 0.0}, {1.0, 0.0}};
}
inline std::vector<double> uniform_gravity_exact(double rate, double T) {
    return {T, -0.5 * rate * T * T, 1.0, -rate * T};
}

inline Spring<2, double> spring_force(double rate) {
    return Spring<2, double>{spatium::Vec<double, 2>{0.0, 0.0}, rate, 0.0};
}
inline PointMass<2, double> spring_y0() {
    return PointMass<2, double>{1.0, {1.0, 0.0}, {0.0, 0.5}};
}
inline std::vector<double> spring_exact(double rate, double T) {
    double w = std::sqrt(rate);
    return {std::cos(w * T), (0.5 / w) * std::sin(w * T),
            -w * std::sin(w * T), 0.5 * std::cos(w * T)};
}

inline PointGravity<2, double> point_gravity_force(double rate) {
    return PointGravity<2, double>{spatium::Vec<double, 2>{0.0, 0.0}, rate};
}
inline PointMass<2, double> point_gravity_y0(double rate) {
    return PointMass<2, double>{1.0, {1.0, 0.0}, {0.0, std::sqrt(rate)}};
}
inline std::vector<double> point_gravity_exact(double rate, double T) {
    double w = std::sqrt(rate);
    return {std::cos(w * T), std::sin(w * T), -w * std::sin(w * T), w * std::cos(w * T)};
}

// Same nondimensionalization (mass=1, reference length/speed=1) makes
// this uniform across all three families -- no per-family branching
// needed, unlike ode_ops.hpp's CircularOrbit-only case.
inline double characteristic_frequency(double rate) { return std::sqrt(rate); }

template<typename Stepper>
std::vector<double> run_family(IntegratorFamily family, double rate, double dt, int steps,
                                Stepper&& stepper) {
    switch (family) {
    case IntegratorFamily::UniformGravity: {
        auto body = uniform_gravity_y0();
        auto force = uniform_gravity_force(rate);
        double t = 0.0;
        for (int i = 0; i < steps; ++i) { stepper(body, force, dt, t); t += dt; }
        return {body.state.position[0], body.state.position[1],
                body.state.velocity[0], body.state.velocity[1]};
    }
    case IntegratorFamily::Spring: {
        auto body = spring_y0();
        auto force = spring_force(rate);
        double t = 0.0;
        for (int i = 0; i < steps; ++i) { stepper(body, force, dt, t); t += dt; }
        return {body.state.position[0], body.state.position[1],
                body.state.velocity[0], body.state.velocity[1]};
    }
    case IntegratorFamily::PointGravity: {
        auto body = point_gravity_y0(rate);
        auto force = point_gravity_force(rate);
        double t = 0.0;
        for (int i = 0; i < steps; ++i) { stepper(body, force, dt, t); t += dt; }
        return {body.state.position[0], body.state.position[1],
                body.state.velocity[0], body.state.velocity[1]};
    }
    }
    return {};
}

inline std::vector<double> exact(IntegratorFamily family, double rate, double T) {
    switch (family) {
    case IntegratorFamily::UniformGravity: return uniform_gravity_exact(rate, T);
    case IntegratorFamily::Spring: return spring_exact(rate, T);
    case IntegratorFamily::PointGravity: return point_gravity_exact(rate, T);
    }
    return {};
}

} // namespace integrator_detail

// in = [family_id, rate, dt, steps], out = [x, y, vx, vy] -- always 4
// slots, unlike ode_ops.hpp's padded variable-width case: every family
// here is a 2D PointMass, so the state is always 4 numbers, no padding
// convention needed.
inline Registry build_integrator_registry() {
    Registry reg;

    auto make_fn = [](auto stepper) {
        return [stepper](std::span<const double> in, std::span<double> out) {
            auto family = static_cast<IntegratorFamily>(static_cast<int>(in[0]));
            auto result = integrator_detail::run_family(family, in[1], in[2],
                                                          static_cast<int>(in[3]), stepper);
            for (std::size_t i = 0; i < out.size(); ++i) out[i] = result[i];
        };
    };

    auto add_stepper = [&](const char* name, auto stepper) {
        reg.add({.name = name,
                 .tier = Tier::General,
                 .in_size = 4,
                 .out_size = 4,
                 .input_names = {"family", "rate", "dt", "steps"},
                 .output_names = {"x", "y", "vx", "vy"}},
                make_fn(stepper));
    };

    // Cheapest (1 force eval/step) to priciest (6 evals/step) -- see the
    // per-method eval counts in the header comment above.
    add_stepper("integrator_euler", [](auto& body, auto&& force, double dt, double t) {
        spatium::physics::mechanics::euler_step(body, force, dt, t);
    });
    add_stepper("integrator_semi_implicit_euler",
                [](auto& body, auto&& force, double dt, double t) {
                    spatium::physics::mechanics::semi_implicit_euler_step(body, force, dt, t);
                });
    add_stepper("integrator_verlet", [](auto& body, auto&& force, double dt, double t) {
        spatium::physics::mechanics::verlet_step(body, force, dt, t);
    });
    add_stepper("integrator_rk4", [](auto& body, auto&& force, double dt, double t) {
        spatium::physics::mechanics::rk4_step(body, force, dt, t);
    });
    add_stepper("integrator_yoshida4", [](auto& body, auto&& force, double dt, double t) {
        spatium::physics::mechanics::yoshida4_step(body, force, dt, t);
    });

    return reg;
}

} // namespace rsc
