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
#include <cassert>

namespace rsc {

// The registry's pinned out_size is 6 -- three roots, re and im each --
// so this copies exactly three and asserts it has three to copy.
//
// An assert rather than a Result, and the reasoning is worth keeping
// because the first answer was the wrong one. The op receives
// `std::span<const double> in` with `in[0]` as the leading coefficient,
// which *looks* like a boundary with untrusted input. It is not. A
// boundary is where data arrives from outside the process and may be
// anything -- load_obj, json::parse. Here it comes from a harness living
// in the same binary, generating inputs from a distribution it defines
// itself. That is an internal pipeline invariant, and conventions.md
// answers those with an assert, the same as physics/mechanics/, spaces/
// and mesh/. A monic cubic always has three roots with multiplicity; if
// it has two, either the input or the solver is broken, and both are
// bugs rather than failure modes.
inline void flatten_cubic_roots(const spatium::UpTo<spatium::Complex<double>, 3>& roots,
                                 std::span<double> out) {
    assert(roots.size() == 3 && "flatten_cubic_roots: a monic cubic has three roots");
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
                assert(in[0] != 0.0 && "harness contract: monic cubic");
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
                assert(in[0] != 0.0 && "harness contract: monic cubic");
                auto roots = spatium::solve_cubic<Real50>(Real50(in[0]), Real50(in[1]),
                                                            Real50(in[2]), Real50(in[3]));
                spatium::UpTo<spatium::Complex<double>, 3> narrowed;
                for (std::size_t i = 0; i < roots.size(); ++i)
                    narrowed.push_back(spatium::Complex<double>{
                        roots[i].re.template convert_to<double>(),
                        roots[i].im.template convert_to<double>()});
                flatten_cubic_roots(narrowed, out);
            });

    return reg;
}

} // namespace rsc
