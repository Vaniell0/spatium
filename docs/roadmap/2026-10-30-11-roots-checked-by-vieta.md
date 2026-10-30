---
date: 2026-10-30
title: the roots of a quadratic checked by Vieta's formulas, the last step of the chain
files: [include/spatium/algebra/polynomial.hpp, include/spatium/io/steps_text.hpp, rsc/spatium/closed_form.hpp, rsc/spatium/respond.hpp, tests/test_solution_steps.cpp]
open:
  - the solver still loses the small root to cancellation when b^2 dwarfs 4ac; the check now says so, the stable form q = -(b + sign(b) sqrt(D)) / 2 is not taken
---

`check_vieta` is a second route to the roots of a quadratic, cheaper than putting them back into the equation: they sum to -b/a and multiply to c/a, each measured relative to the size of what it adds or multiplies. It holds to 1e-12 on 20000 integer quadratics, and it fails where it must: for x^2 + 1e8 x + 1 the small root comes out -7.45e-9 against -1e-8 and the product of the roots is 0.745 instead of 1, the cancellation the solver's comment records as not addressed. The check is the last step of an equation's chain (`x_{1} + x_{2} = 1.5 = -\frac{b}{a} = 1.5, x_{1} x_{2} = -2.5 = \frac{c}{a} = -2.5`, or `\ne` with the relative error) and its verdict is in the answer's JSON.
