---
date: 2026-11-16
title: an identity is not a quadratic of noise
files: [include/spatium/answer/closed_form.hpp, tests/test_answer_evaluate.cpp]
open:
  - an identity is refused with its reason, not answered as "every a": the evaluation has no answer of that kind yet
  - an equation whose values are all below 1e-12 in doubles, 1e-20 x^2 = 1e-20, would read as one; the exact path answers it first
---

`sin^2(a - pi/6) + sin^2(a + pi/6) - sin^2(a) = 1/2` holds for every a, and the closed form answered a = 1, a = -1: its values at 0, 1 and -1 are rounding, so the quadratic through them was 1.1e-16 a^2 - 1.1e-16, which passed the degree check at five more points and Vieta's formulas. Zero at all eight points to 1e-12 is now an identity, and the closed form refuses it by name: "an identity, true for every a where it is defined". `x^2 - 1 = 0` is answered as before (`test_answer_evaluate`).
