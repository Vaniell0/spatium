---
date: 2026-11-03
title: in fractions, a closed form's roots are what they are
files: [include/spatium/core/rational.hpp, rsc/spatium/closed_form.hpp, rsc/spatium/query.hpp, rsc/spatium/evaluate.hpp, rsc/spatium/respond.hpp, tests/test_rsc_evaluate.cpp, site/wasm/questions.txt]
open:
  - the steps before the roots (coefficients, discriminant) still print the double's numbers
  - a linear equation's step prints its root as a double; its exact root is in `exact`
  - no exact route for a cubic or a quartic
---

With the engine's arithmetic set to fractions, `2x^2 - 3x - 5 = 0` answers 5/2 and -1, and `x^2 - 3x - 2 = 0` answers (3 ± √17)/2 instead of 3.5615528128088303 and -0.5615528128088303.

Nothing new solves it. `Rational` gained the square root it has: of 9/4 it is 3/2, and of 2 it is spoiled. The library's own `solve_quadratic_steps`, a template that finds `sqrt` by ADL, run on fractions then gives exact roots wherever the discriminant is a square. Where it is not, the roots are written (p ± q√m)/d, with m free of square factors (±2√2, (1 ± √7)/3, -1 ± 2i). The coefficients are computed again as fractions: the unknown is a slot of the same lowered program, run in `Rational`, and the degree is checked by equality at five more points, not to a tolerance.

One trap was found and closed on the way. A root's equation is simplified before it is solved, and simplifying folds every constant subtree into a double: 1/3 became 0.3333333333333333, and sqrt(2) became 1.4142135623730951. A fraction read from those is a decimal approximation that looks exact. So the exact route reads the equation as written (`Unfolded::written`). Now 1/3 is a third (`3x^2 - 1/3 = 0` gives 1/3 and -1/3), and an equation that reads `sqrt` is refused by name ("the equation reads Sqrt, which has no fraction"), with the double's closed form standing.

The pair: p ± q√m computed as doubles equals the double solver's roots within 4e-16 relative. The WebAssembly core and the host agree on 38 of 38 questions.
