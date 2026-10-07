---
date: 2026-10-06
title: Quadrature over a domain that is a type
open:
  - integrate_surface
  - a Monte Carlo estimator
  - Certified<T>: a proven bound instead of an estimate (stage 4)
  - a NaN in the integrand inside the domain is read as the end of that side: the value is a wrong number and the status DepthCapped, where Failed with NaN is honest
  - the doubly exponential rule has a level cap of 9 and a t range capped at 10, and no subdivision for an integrand with a feature away from where its nodes crowd
  - Gauss-Kronrod bisects to depth 20 and halves the tolerance per half instead of keeping one error budget as QUADPACK does, and its estimate |K15 - G7| is unscaled
  - a singularity at an end b != 0 is truncated at an ulp of b
  - Gauss-Kronrod beyond long double: its constants stop at 33 digits
resolves: [2026-10-06-04#1, 2026-10-06-04#2, 2026-10-06-04#3]
stage: 2.7
---
`algebra/quadrature.hpp`: `quadrature(f, Finite<T>{a, b})` (tanh-sinh), `quadrature(f, HalfLine<T>{from, toward_infinity})` (exp-sinh) and `quadrature<T>(f, WholeLine{})` (sinh-sinh), each returning the `IntegralResult` of `integrate_with_error`, and `gauss_kronrod(f, a, b)` (adaptive G7K15, scalars up to long double; a `static_assert` refuses Real50 rather than answer to 19 digits) as the second witness. The nodes are formulas in the scalar, so the same call is right on float, double, Real50 and `Dual`, with the tolerance following the scalar (eps^(3/4)); the doubly exponential rule's estimate is the difference of its last two levels and its status reads how that difference fell. What the audit measured as NaN is answered: 1/sqrt(x) and log x on [0, 1] exactly (56 evaluations, estimate 3e-15 and 2e-13), exp(-x) on [0, inf), exp(-x^2) on the line and exp(-x)/sqrt(x) on [0, inf) with the singular end and the infinite end together, to the last digit of a double; on Real50 log x over [0, 1] to 7e-51, exp(-x) over [0, inf) to 1e-52, exp(-x^2) over the line to 3e-52. The derivative passes through the domain: d/dp of the integral of exp(-p x) over [0, inf) at p = 2 is -1/4 on `Dual`, and d/dp of the integral of x^p over [0, 1] is -1/(p+1)^2. Held by: closed forms, the rule against Simpson against Gauss-Kronrod (the difference of two witnesses no larger than the sum of their estimates), and the cases that must NOT be called converged -- 1/x on (0, 1], exp(x) on the line, NaN bounds, a function that is NaN everywhere (which first came back as a converged 0: a difference of two zeros). Limits stated, not hidden: a node near an end b != 0 is a point of the scalar, so a singularity AT such an end is truncated at an ulp of b (the near-1 end of 1/sqrt(x(1-x)) gives pi to 1e-6, not 1e-12); Gauss-Kronrod's bisection does not conquer an endpoint singularity (sqrt x and x^0.7 on [0, 1] come back `DepthCapped`, with a true error of 1e-14 and 1e-16, 630 evaluations) -- which is the case tanh-sinh is the other witness for. Not done: `integrate_surface`, a Monte Carlo estimator, `Certified<T>` (a proven bound, stage 4).
