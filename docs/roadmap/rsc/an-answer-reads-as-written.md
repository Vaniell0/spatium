---
date: 2026-11-03
title: a synthesised answer reads as a person writes it, where it still holds
files: [rsc/spatium/synthesis.hpp, tests/test_rsc_synthesis.cpp, site/wasm/questions.txt]
open:
  - only a line in x is rewritten; a quadratic argument, or a product of terms, keeps the interval's variable
  - a coefficient goes under a square root only; under a logarithm it would be a sum, not a factor
---

The root of x^2 + p = 0 for p in [-1, 0] used to come out `0.7071067811865477 * sqrt(1 - (p - (-0.5)) * 2)`, which is sqrt(-p). Synthesis builds its terms over the interval's own variable, t = (x - centre) * 2 / width, which keeps powers apart on a narrow interval. An accepted answer was printed in that variable, with its coefficient left outside.

Once an answer is accepted, it is now rewritten one term at a time:
- a term that is a line in x, or a function whose argument is one, is written b x + a, its numbers as fractions where they are;
- a coefficient is taken into the line, or under a square root as its square when what it leaves there is short; otherwise it stays outside, so π√(-p) does not become √(-9.87 p).

Each rewrite is kept only if the whole answer still holds on the check grid and is no longer than it was. This is the rule the fraction coefficients already followed. The question above now answers `sqrt(-p)`, and sqrt(p + 1) on [-1, 0] reads `\sqrt{p + 1}`. The gallery's answers were already short and did not change. The WebAssembly core and the host agree on 39 of 39 questions, the question above among them.
