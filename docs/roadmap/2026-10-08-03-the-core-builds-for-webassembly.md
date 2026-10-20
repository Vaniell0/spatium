---
date: 2026-10-08
title: The core builds for WebAssembly, and gives the host's digits
stage: 7
files: [include/spatium/algebra/vec_simd.hpp, .github/workflows/ci.yml]
open:
  - the probe was a small program (the cubic and quartic solvers with the polished roots, and `integrate` over a finite interval and the whole line); the viewer, Boost.Multiprecision (so `Real50`), Eigen and the whole test suite have not been built for WebAssembly
  - the comparison was made by hand, host against `node`; it is not yet a pair in the registry, and no CI job builds with Emscripten
  - the size is one measurement (30 KB of .wasm for that program, 7.7 KB of JavaScript); what a growing menu costs is not budgeted
  - a floating-point result is identical here because both use IEEE doubles and the same rounding; a function that reaches for a libm routine whose result differs between the host's library and Emscripten's would not be, and no pair has asked about one
---

The web page of the solver is gated on one thing: that the core runs in the browser. It does. Emscripten 6.0.9 (from nixpkgs) compiles a program that calls `solve_cubic`, `solve_quartic` and `integrate` to 30 KB of WebAssembly and 7.7 KB of JavaScript, and under `node` it prints the same text as the native build to the last digit, the quartic with roots 1/4, 2, 64 and 4096 and the integral of sin over [0, pi] (1.9999999999999998 in 112 evaluations) and of exp(-x^2) over the line (1.7724538509055159) included.

It did not compile at first. `Vec<float,4>` names `simd::add_f4`, `sub_f4`, `mul_scalar_f4`, `div_scalar_f4`, `dot_f4` and `cross_f4` in a branch that an `if constexpr` never takes, and they existed only under `__SSE2__`; the doubles had stubs for the other targets, the floats did not. So the library did not build on any target without SSE2: WebAssembly, and ARM, which is every Apple-silicon laptop. The stubs are there now, and a CI step compiles the umbrella header with `-U__SSE2__` on the x86 runner, which reproduces the failure (five errors on the previous header) and holds the fix.
