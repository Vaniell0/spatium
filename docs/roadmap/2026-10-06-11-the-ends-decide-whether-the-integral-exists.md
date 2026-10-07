---
date: 2026-10-06
title: The ends decide whether the integral exists
open:
  - integrate_surface
  - the kappa family and Minkowski<N> in geometric units
  - only the ends are read: a singularity in the interior is not seen
  - a conditionally convergent integral (sin x / x) and a decay that is not a power (exp(-x)) are Undetermined, not decided
  - f must be generic over the scalar, and no concept gives a readable diagnostic when it is not (testing invocability would hard-error inside a generic lambda's body)
stage: 2.9
---
`improper_ends(f, domain)` and `quadrature_checked(f, domain)` in `algebra/quadrature.hpp`, the first consumer of `Series`. Near a finite end at distance d, f ~ c d^p is integrable iff p > -1; at infinity f ~ c x^-p iff p > 1 (the p-test): the leading term alone, which dominates a one-sided neighbourhood, so where the series gives a verdict it is a proof and not an estimate. `quadrature_checked` runs it first and returns the new status `IntegralStatus::Divergent` with NaN and no evaluations spent where the ends prove the integral does not exist (1/x on [0, 1], 1/x and 1/sqrt(x) on [1, inf), x on any half line, x/(1 + x^2) on the line) -- where `quadrature` alone would have run its levels out and said `DepthCapped`; the verdicts also carry the order (1/sqrt(x) on [0, 1] has order -1/2, 1/x^2 at infinity 2), and x^-1.5 at infinity is integrable by ramification 2. Where the series cannot say it says `Undetermined` and the numerical rule goes on alone with its own answer: log x at 0 (-1), exp(-x) at infinity (decays faster than any power, not a Laurent series), sin(x)/x (conditionally convergent). It reads the ends and only the ends -- a singularity in the interior is not seen -- and needs f callable on a Series as well as on a T, a generic function with math through ADL; a plain `double(double)` function still goes to `quadrature`. Not done: `integrate_surface`, the kappa family and `Minkowski<N>` in geometric units.
