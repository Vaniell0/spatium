---
date: 2026-11-21
title: Hyperbolic log with no cancellation, a point's log at itself is 0
files:
  - include/spatium/spaces/hyperbolic.hpp
  - tests/test_hyperbolic.cpp
resolves:
  - spaces/hyperbolic-distance-of-a-point-from-itself-is-zero#1
---

Chains of laws that each hold were replayed as rewrites (eight premises, a start term, every hop an equality of the
snapshot) and every hop checked again on every space. Ten chains reached a term the oracle refuses, and eight of
them through one hop on Hyperbolic<1>, <2> and <4> and their products with S2 and E3: d(q, exp_q(log_q q)) = d(q, q),
on the breaker's tuples whose q is 8 from p. Real50 on the same inputs gives 0 on both sides; double gave 3e-5 to
4e-4, 1e11 times more than an ulp of the inputs moves the exact answer.

`log_map` read d off acosh(-<p,q>), which `distance` stopped doing for the same reason: once a point is a few units
out, -<q,q> is 1 plus rounding, so log_q q was a vector 8e-5 long 9 units from the origin, and not tangent at q.
The same cancellation spoils every short way far out: log_q(exp_q v) against v, |v| = 1e-4, was 3.5e-2 wrong 8 from
the origin and 8.8 (no digit) at 12.

`cosh_minus_one` is now the one place both read: Lagrange's sum of squares over a0 b0 + a.b + 1. `distance` is
2 asinh sqrt(c / 2) of it as before, and `log_map` is proj (d / sinh d) with proj = q + p <p,q> = (q - p) - p c and
sinh d = sqrt(c (c + 2)), so nothing cancels in either. Equal points give the zero tangent at any distance; the
error of log_q(exp_q v) at |v| = 1e-4 is 1e-13 at 0.5 from the origin, 2e-12 at 3, 1e-9 at 8 and 2e-9 at 12, the
rounding of the inputs themselves. All eight chains hold. The law snapshot is unchanged (637 laws on H1, H2, H4,
S2xH2 and H2xE3).
