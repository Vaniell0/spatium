---
date: 2026-11-27
title: a law over data takes the powers its caller names, and gives no law for no data
files: [include/spatium/answer/law.hpp, include/spatium/answer/synthesis.hpp, tests/test_answer_law.cpp]
open:
  - the exponent of a named power is written x^((-5) / 3), the fraction built from two constants, not x^(-5/3)
  - a damped sinusoid, exp(-x / tau) cos(2 pi f x), has two parameters that sit inside its operations; the terms are built from x, 1 and the rules, so f and tau can be only small whole numbers in a chosen unit, and a resonance in the data is not found
  - the standard errors take the points as independent; a frequency track read from overlapping windows has a residual of 0.3 to 0.5 sigma, and the length, which depends on the residual, then favours whichever of log(x), sqrt(x) and x fits the smoother curve
---

`synthesize_law` read the end of an empty vector when it was given no data: an empty file or a first argument that is no file ended in a segmentation fault. It now returns a `Law` with no terms and an infinite length, as it does for data that do not line up (y or sigma of another length than x) and for a point whose noise is not a positive finite number, which would weigh it infinitely, or whose x or y is not finite.

A law over data in a frequency, a phase in f^(-5/3), f^(-2/3), is a sum of powers with a third in the exponent, and no expression of x and 1 under + - * / and log, sqrt, sin, cos, exp in five nodes is one. `LawSpec::powers` (and `SynthesisSpec::powers`) lets the caller who knows which exponents can occur name them: each is a term of two nodes, the exponent built as the fraction it is where it is a small one. Held on 3 x^(-5/3) + 0.5 x^(-2/3) from 20 to 400 with noise 1e-5: with the powers named the law is found, the residual is within 1.5 sigma and its values are within 1e-3 of the truth at three points; with none, no two terms of the other rules reach the noise (residual over 5 sigma) and the length is larger.
