---
date: 2026-10-30
title: a synthesised coefficient that is a small fraction is written as one, and certified again
files: [rsc/spatium/synthesis.hpp, tests/test_rsc_synthesis.cpp]
open:
  - fractions up to a denominator of 1000 and to nine digits; a coefficient such as pi or sqrt 2 stays a number, since the operations have no constant for it
  - each coefficient is tried alone and the others are not refitted around it
  - a piece of a piecewise answer is still written in its interval's variable
---

A fitted coefficient is a number, and a small fraction near it is a different number, so writing 0.16666666666666663 as 1/6 is not a rewrite of the answer but a new answer. The synthesis now makes it one: after an answer is accepted, each coefficient that is a fraction p/q (q up to 1000) to nine digits -- the first such convergent of its continued fraction, the best approximation of its size -- is replaced by that fraction, and the replacement is kept only if the whole answer still holds on the check grid. One at a time, so a fraction that fails costs the others nothing.

On the gallery's questions: the root of x^2 - p reads `\sqrt{p}` (was `1 \cdot \sqrt{p}`), the integral of x^p over [0, 1] reads `\frac{1}{p + 1}` (was `1.0000000000000004 \cdot \frac{1}{p + 1}`), and the prism's volume reads `\frac{1}{6} \cdot p + 1`, certified to 1.33e-15 against 1.78e-15 for the fitted form. A coefficient of sqrt 2 keeps its digits.
