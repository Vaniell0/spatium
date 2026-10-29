---
date: 2026-10-29
title: the chance a walk on Z^3 comes back, for walks whose steps are one orbit of the cube's symmetries
files: [tests/test_return_probability.cpp, tests/CMakeLists.txt]
open:
  - the library has no summation of a slowly falling series with a bound on its tail; the series route here knows its asymptotic by hand and is good to 3e-9 at a thousand terms
  - the integral's error estimate is not a bound: at full precision it is below the true error by a factor of two to five (1.5e-16 against 2.9e-16, 1.2e-15 against 5.7e-15)
  - near the corner singularity the innermost integrals stop at the depth cap (about 1300 of them, all within 5e-21 of the corner); the outer integral converges because their weight is negligible, which nothing in the result states
  - walks in other dimensions and on orbits with more than one nonzero length (the knight's (1,2,0)) have no closed form here to be held to
---

A walk whose steps are the orbit of a vector v under the signed permutations of three axes comes back to its start with probability p = 1 - 1/W, where W is Watson's integral of 1/(1 - phi(k)) over [0, pi]^3 / pi^3 and phi(k) the mean of cos(k.s) over the steps. The steps come from `algebra::orbit` over `SignedPermutationAction`; three orbits are the three cubic lattices: (1,0,0) the simple cubic with six steps, (1,1,1) the body-centred with eight, (1,1,0) the face-centred with twelve. Each W has a closed form in Gamma functions, and three nested `integrate` calls, the doubly exponential rule meeting the singularity at the corner of each interval, give it: at the rule's own tolerance to 2.9e-16, 5.7e-15 and 5.0e-15 relative (about 55 s at -O2 in all, under `[.deep]`), at tolerance 1e-3 to 4e-7..9e-7 within estimates of 2e-5..5e-5 (0.28 s, in the default run). The return probabilities are 0.340537329551 (Polya's), 0.282229988954 and 0.256318236505.

Written as 1 - (cos a + cos b + cos c)/3 the denominator cancels to zero below k of about 1e-8 and the integrand is infinite inside the cube; the quadrature returns NaN with a status that is not trusted, and the test holds it to that. Written as the mean of 2 sin^2(k.s / 2), which is the same function without the cancellation, it converges.

For the simple cubic walk the series of return probabilities after 2n steps is the second route. Its terms, C(2n,n) S(n) / 36^n with S(n) the sum of squared trinomial coefficients (two forms of S agree to 1.4e-14), fall like n^(-3/2): a thousand terms are 1.5e-2 short. The leading asymptotic 2 (3/(4 pi n))^(3/2), its tail summed exactly through zeta(3/2), leaves 1.8e-6, and Richardson in N leaves 2.6e-9; doubling N divides the two errors by 2.9 and 5.7, the 2^1.5 and 2^2.5 the asymptotic predicts, and the test holds both ratios. The integral reaches in 3 s what the series, knowing one term of its asymptotic, does not reach at all.
