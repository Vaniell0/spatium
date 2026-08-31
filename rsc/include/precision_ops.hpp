#pragma once

// The first Tier-2 domain: polynomial/precision-critical dispatch. Unlike
// Tier-1, both registered ops compute (nearly) the same mathematical answer
// -- solve the same cubic -- at different precision. Dispatch here isn't
// "which op explains this input/output pair" (Tier-1's shape); it reads a
// property of the input alone and decides whether double is accurate enough
// or Real50 needs to be paid for, matching the design doc's domain-1
// description directly. See precision_task.hpp for how "correct" is defined.
//
// Required fixing a real bug in Spatium itself first: solve_cubic<Real50>
// didn't compile at all before this session (qualified std:: math calls
// bypassing ADL, std::numbers::pi_v<T> restricted to std:: float types, and
// Boost.Multiprecision's number<> arithmetic returning lazy expression
// types that broke template deduction across solve_quadratic/solve_cubic/
// solve_quartic's cross-calls) -- see polynomial.hpp/complex.hpp and
// tests/test_polynomial.cpp's Real50 cases.

#include <registry.hpp>
#include <spatium/algebra/polynomial.hpp>
#include <spatium/core/precision.hpp>
#include <array>

namespace rsc {

inline void flatten_cubic_roots(const std::array<spatium::Complex<double>, 3>& roots,
                                 std::span<double> out) {
    for (int i = 0; i < 3; ++i) {
        out[2 * i] = roots[i].re;
        out[2 * i + 1] = roots[i].im;
    }
}

inline Registry build_precision_registry() {
    Registry reg;

    reg.add({.name = "solve_cubic_f64",
             .tier = Tier::Library,
             .in_size = 4,
             .out_size = 6,
             .input_names = {"a", "b", "c", "d"},
             .output_names = {"re0", "im0", "re1", "im1", "re2", "im2"}},
            [](std::span<const double> in, std::span<double> out) {
                auto roots = spatium::solve_cubic<double>(in[0], in[1], in[2], in[3]);
                flatten_cubic_roots(roots, out);
            });

    reg.add({.name = "solve_cubic_real50",
             .tier = Tier::Library,
             .in_size = 4,
             .out_size = 6,
             .input_names = {"a", "b", "c", "d"},
             .output_names = {"re0", "im0", "re1", "im1", "re2", "im2"}},
            [](std::span<const double> in, std::span<double> out) {
                using spatium::Real50;
                auto roots = spatium::solve_cubic<Real50>(Real50(in[0]), Real50(in[1]),
                                                            Real50(in[2]), Real50(in[3]));
                std::array<spatium::Complex<double>, 3> narrowed;
                for (int i = 0; i < 3; ++i) {
                    narrowed[i].re = roots[i].re.convert_to<double>();
                    narrowed[i].im = roots[i].im.convert_to<double>();
                }
                flatten_cubic_roots(narrowed, out);
            });

    return reg;
}

} // namespace rsc
