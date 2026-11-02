---
date: 2026-10-07
title: The first space with no finite basis, L^2 of an interval
stage: 2.6
files: [include/spatium/spaces/l2.hpp, tests/test_l2.cpp]
open:
  - the Chebyshev basis, orthogonal for a weight, so a different inner product
  - a point as a lazy Field expression rather than a truncated series
  - the Fourier basis for periodic functions
  - ProductSpace with a function space (its dimension is infinite, the points are pairs)
  - from_function costs n integrals by the door (a few hundred evaluations each): no fast transform
  - the tail estimate is the L2 mass of the last three coefficients: a heuristic, not a bound (Certified<T>, stage 4)
---
Stage 2.6, the consumer of the infinite dimension. `spaces/l2.hpp`: `L2Interval<T>{a, b}` is the space of square-integrable functions on an interval, `dimension = Dimension::infinite()`; a point is a `LegendreSeries<T>`, a truncated series sum c_n P_n(t) with t = (2x - a - b)/(b - a). The Legendre polynomials are orthogonal for the plain L^2 weight, so the inner product is Parseval's identity on the coefficients, (b - a)/2 sum 2/(2n+1) c_n d_n, with no quadrature -- and the norm of f^2 by the door of `integrate` is an independent path to the same number: they agree to 5e-16 relative for exp, cos 3x and 1/(1+x^2) on [-1, 2] with 40 terms. It is a `HilbertSpace` (inner product, complete) and not a `EuclideanSpace` (no finite basis), asserted at compile time; `verify_metric` passes on four sampled functions, the polynomials P_0..P_4 are orthogonal with the known norms (b - a)/(2n + 1), Pythagoras and Cauchy-Schwarz hold. A series reproduces an analytic function (24 terms, max error 2e-12 over [-1, 2] for exp(sin x) + 0.3 x). A series carries an estimate of what it neglected -- the L^2 mass of its last three coefficients -- and the estimate says when the projection cannot be trusted: 1e-15 for exp at 24 terms, 3.2e-3 for |x|, whose coefficients decay like 1/n^2 and not geometrically. It is an estimate, not a bound; the bound is what `Certified<T>` is for (stage 4). Built from a function by `from_function(f, n)`, the coefficients (2k + 1)/2 times the integrals of f P_k by the door. Not done: the Chebyshev basis (orthogonal for a weight, so a different inner product), a point as a lazy `Field` expression rather than a series, the Fourier basis for periodic functions, and `ProductSpace` with a function space (its dimension is infinite; the points are pairs).
