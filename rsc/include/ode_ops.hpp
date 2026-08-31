#pragma once

// Domain 4: Cauchy/IVP method dispatch (Euler vs. RK4), subsuming what the
// settled pipeline called "orbital mechanics" -- not a separate domain,
// one more test-problem family below. Correction made while building
// this: physics/orbital.hpp (referenced by that pipeline entry) turned
// out to be *atomic* orbital theory (Legendre polynomials, Slater's
// rules for electron shells), not celestial mechanics -- there was no
// existing closed-form two-body solution to reuse. Rather than pull in
// eccentric-orbit Kepler's-equation solving (itself a root-finding
// problem, ironically), this uses the circular special case, whose
// closed form is elementary -- a real, if restricted, orbital-mechanics
// family, not a placeholder.
//
// Three test families sharing one dispatch mechanism, matching
// spatium/algebra/ode.hpp's genericity over Vec<T,N>:
//   - Decay:    y' = -k*y            (N=1, exact: y0*exp(-k*t))
//   - Oscillator: [x,v]' = [v,-w^2*x] (N=2, exact: closed-form sinusoid)
//   - CircularOrbit: 2D gravity, circular initial velocity (N=4, exact:
//     closed-form circular motion at w = sqrt(GM/r0^3))
// All three share the same qualitative signal: a characteristic
// frequency/rate times the step size (freq*dt) predicts whether Euler's
// bounded stability region still covers this step, or RK4's cost is
// actually needed -- see spatium/algebra/ode.hpp's euler_step doc comment
// for the boundary itself, verified directly in tests/test_ode.cpp.

#include <registry.hpp>
#include <spatium/algebra/ode.hpp>
#include <cmath>
#include <vector>

namespace rsc {

enum class OdeFamily { Decay = 0, Oscillator = 1, CircularOrbit = 2 };

namespace ode_detail {

constexpr double kOrbitRadius = 1.0;

inline auto decay_rhs(double k) {
    return [k](double, spatium::Vec<double, 1> y) { return spatium::Vec<double, 1>{-k * y[0]}; };
}
inline spatium::Vec<double, 1> decay_y0() { return {1.0}; }
inline std::vector<double> decay_exact(double k, double t) { return {std::exp(-k * t)}; }

inline auto osc_rhs(double w) {
    return [w](double, spatium::Vec<double, 2> y) {
        return spatium::Vec<double, 2>{y[1], -w * w * y[0]};
    };
}
inline spatium::Vec<double, 2> osc_y0() { return {1.0, 0.0}; }
inline std::vector<double> osc_exact(double w, double t) {
    return {std::cos(w * t), -w * std::sin(w * t)};
}

inline auto orbit_rhs(double GM) {
    return [GM](double, spatium::Vec<double, 4> y) {
        double x = y[0], py = y[1];
        double r = std::sqrt(x * x + py * py);
        double r3 = r * r * r;
        return spatium::Vec<double, 4>{y[2], y[3], -GM * x / r3, -GM * py / r3};
    };
}
inline spatium::Vec<double, 4> orbit_y0(double GM) {
    double v0 = std::sqrt(GM / kOrbitRadius);
    return {kOrbitRadius, 0.0, 0.0, v0};
}
inline std::vector<double> orbit_exact(double GM, double t) {
    double w = std::sqrt(GM / (kOrbitRadius * kOrbitRadius * kOrbitRadius));
    return {kOrbitRadius * std::cos(w * t), kOrbitRadius * std::sin(w * t),
            -kOrbitRadius * w * std::sin(w * t), kOrbitRadius * w * std::cos(w * t)};
}

// Characteristic frequency/rate per family -- what freq*dt means for the
// dispatch feature. For CircularOrbit, GM is the sampled parameter but
// w=sqrt(GM/r0^3) (=sqrt(GM) since r0=1) is the actual timescale.
inline double characteristic_frequency(OdeFamily family, double rate) {
    switch (family) {
    case OdeFamily::Decay: return rate;
    case OdeFamily::Oscillator: return rate;
    case OdeFamily::CircularOrbit: return std::sqrt(rate);
    }
    return rate;
}

template<typename Stepper>
std::vector<double> run_family(OdeFamily family, double rate, double dt, int steps,
                                Stepper&& stepper) {
    switch (family) {
    case OdeFamily::Decay: {
        auto f = decay_rhs(rate);
        auto y = spatium::integrate_fixed(f, 0.0, decay_y0(), dt, steps, stepper);
        return {y[0]};
    }
    case OdeFamily::Oscillator: {
        auto f = osc_rhs(rate);
        auto y = spatium::integrate_fixed(f, 0.0, osc_y0(), dt, steps, stepper);
        return {y[0], y[1]};
    }
    case OdeFamily::CircularOrbit: {
        auto f = orbit_rhs(rate);
        auto y = spatium::integrate_fixed(f, 0.0, orbit_y0(rate), dt, steps, stepper);
        return {y[0], y[1], y[2], y[3]};
    }
    }
    return {};
}

inline std::vector<double> exact(OdeFamily family, double rate, double t) {
    switch (family) {
    case OdeFamily::Decay: return decay_exact(rate, t);
    case OdeFamily::Oscillator: return osc_exact(rate, t);
    case OdeFamily::CircularOrbit: return orbit_exact(rate, t);
    }
    return {};
}

} // namespace ode_detail

// Registered for real usable dispatch (the "апишку не забыть" API
// surface), not just internal task-generator closures. in/out shape:
// in = [family_id, rate, dt, steps], out = final state, padded to 4
// slots (CircularOrbit's max) -- unused trailing slots left at 0 for the
// smaller families, same padding convention Dispatcher::features() uses.
inline Registry build_ode_registry() {
    Registry reg;

    auto make_fn = [](auto stepper) {
        return [stepper](std::span<const double> in, std::span<double> out) {
            auto family = static_cast<OdeFamily>(static_cast<int>(in[0]));
            auto result = ode_detail::run_family(family, in[1], in[2], static_cast<int>(in[3]),
                                                  stepper);
            for (std::size_t i = 0; i < out.size(); ++i)
                out[i] = i < result.size() ? result[i] : 0.0;
        };
    };

    auto euler_stepper = [](auto&& f, double t, auto y, double dt) {
        return spatium::euler_step(f, t, y, dt);
    };
    auto rk4_stepper = [](auto&& f, double t, auto y, double dt) {
        return spatium::rk4_step(f, t, y, dt);
    };

    reg.add({.name = "ode_euler",
             .tier = Tier::General,
             .in_size = 4,
             .out_size = 4,
             .input_names = {"family", "rate", "dt", "steps"},
             .output_names = {"y0", "y1", "y2", "y3"}},
            make_fn(euler_stepper));

    reg.add({.name = "ode_rk4",
             .tier = Tier::General,
             .in_size = 4,
             .out_size = 4,
             .input_names = {"family", "rate", "dt", "steps"},
             .output_names = {"y0", "y1", "y2", "y3"}},
            make_fn(rk4_stepper));

    return reg;
}

} // namespace rsc
