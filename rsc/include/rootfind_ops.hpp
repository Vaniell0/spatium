#pragma once

// Domain 3: root-finding method dispatch (Newton vs. bisection). Cheapest
// of the domains found still open this session -- reuses algebra/dual.hpp's
// Dual<T> exactly as calculus.hpp does (a generic-over-scalar-type functor,
// Dual-seeded for the derivative), no new Spatium primitive needed.
//
// Test function: f(x) = x^3 - a, root = cbrt(a). Chosen deliberately, not
// arbitrarily -- it has an exact closed-form reference (std::cbrt) and a
// single real inflection at x=0 where f'(x)=3x^2 vanishes, which is
// Newton's actual, textbook failure mode (huge first step, many iterations
// to recover), not a contrived edge case. See rootfind_task.hpp for how
// the resulting two-regime behavior becomes the dispatch signal.

#include <registry.hpp>
#include <spatium/algebra/dual.hpp>
#include <cmath>

namespace rsc {

template<typename S>
S cube_minus_a(S x, S a) { return x * x * x - a; }

// Fixed iteration budget, not adaptive convergence checking -- deliberately
// small enough (20) that a badly-conditioned start (x0 near the f'=0
// inflection) measurably fails to converge in time, giving the dispatcher
// a real signal instead of both candidates always succeeding.
inline double newton_root(double x0, double a, int max_iters = 20) {
    double x = x0;
    for (int i = 0; i < max_iters; ++i) {
        spatium::Dual<double> dx = spatium::Dual<double>::variable(x);
        spatium::Dual<double> da(a);
        spatium::Dual<double> d = cube_minus_a(dx, da);
        if (d.deriv == 0.0 || !std::isfinite(d.value) || !std::isfinite(d.deriv)) break;
        x = x - d.value / d.deriv;
        if (!std::isfinite(x)) break;
    }
    return x;
}

// Bracket must satisfy cube_minus_a(lo,a) and cube_minus_a(hi,a) having
// opposite signs -- rootfind_task.hpp's sample_problem is responsible for
// that, not this function.
inline double bisection_root(double lo, double hi, double a, int max_iters = 50) {
    double flo = cube_minus_a(lo, a);
    for (int i = 0; i < max_iters; ++i) {
        double mid = 0.5 * (lo + hi);
        double fmid = cube_minus_a(mid, a);
        if ((fmid < 0.0) == (flo < 0.0)) {
            lo = mid;
            flo = fmid;
        } else {
            hi = mid;
        }
    }
    return 0.5 * (lo + hi);
}

inline Registry build_rootfind_registry() {
    Registry reg;

    reg.add({.name = "newton_root",
             .tier = Tier::General,
             .in_size = 2,
             .out_size = 1,
             .input_names = {"x0", "a"},
             .output_names = {"root"}},
            [](std::span<const double> in, std::span<double> out) {
                out[0] = newton_root(in[0], in[1]);
            });

    reg.add({.name = "bisection_root",
             .tier = Tier::General,
             .in_size = 3,
             .out_size = 1,
             .input_names = {"lo", "hi", "a"},
             .output_names = {"root"}},
            [](std::span<const double> in, std::span<double> out) {
                out[0] = bisection_root(in[0], in[1], in[2]);
            });

    return reg;
}

} // namespace rsc
