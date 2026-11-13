---
date: 2026-11-13
title: Hyperbolic distance with no cancellation, a point is 0 from itself
files:
  - include/spatium/spaces/hyperbolic.hpp
  - tests/test_hyperbolic.cpp
open:
  - log_map still reads acosh of the inner product to get d, so the length of log_p(q) for q close to p and far from the origin keeps the rounding of -<p,q> - 1 (it returns the zero tangent below epsilon, as before)
  - ConstantCurvature at kappa < 0 computes its distance in the Poincare ball and was not touched
---

The law finder reported that d(p,p) = d(exp(p,w), exp(p,w)) holds on the Euclidean, spherical and product spaces and not on Hyperbolic<1>, <2> and <4>, with the two further laws d(p,q) = d(p,q) + d(q,q) and d(q,r) = d(q,q) + d(q,r) that follow from it. The left side is 0 exactly; the right side was not.

`Hyperbolic::distance` was acosh(-<a,b>). For a point carried a few units from the origin by `exp_map`, -<a,a> is 1 plus the rounding of terms of size e^{2d} that cancel, and acosh near 1 turns an excess of 2e-16 into 2e-8 at best. Measured on H^2 with a unit step repeated along one direction, d(p,p) was 0 at 0.5 from the origin, 2.4e-7 at 3 and at 8, and 8.3e-7 at 12.

It is now computed as cosh d - 1 directly. On the sheet (a0 b0)^2 - (a.b + 1)^2 = |a - b|^2 + sum over i < j of (a_i b_j - a_j b_i)^2, Lagrange's identity as in `tangent_norm_sq`, so cosh d - 1 is that sum of squares over a0 b0 + a.b + 1, a denominator of at least 1 since a0 b0 >= |a||b| >= |a.b|, and d = 2 asinh sqrt((cosh d - 1)/2). Equal points give exactly 0 at every distance tried (0.5 to 12 on H^2, 9 on H^4), and two points on one geodesic through the origin, 7 and 7.5, 10 and -4, 1e-4 and 2e-4 from it, agree with |s - t| to 1e-9 relative. A NaN coordinate still gives NaN, with no comparison or max() in the path, which is what the NaN test over Real50 holds.

The snapshot of laws, discovered again: 834 laws became 637, the 226 that went being inequalities among terms equal to d(q,q) that are now the term 0 and merge with it, 29 new ones are the same inequalities against d(p,exp(p,w)) and the like. No law that held on a space stopped holding on it, and the open tasks on H1, H2, H4 and S2xH2 closed.
