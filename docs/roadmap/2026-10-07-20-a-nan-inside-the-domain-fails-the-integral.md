---
date: 2026-10-07
title: A NaN inside the domain fails the integral
stage: 2.7
files: [include/spatium/algebra/quadrature.hpp, tests/test_quadrature.cpp]
resolves: [2026-10-06-09#4]
open:
  - a value that is not finite is told from an overflow at an end by one comparison, the last term of the side against the tolerance times what the integral had summed; an integrand that decays slowly and then goes NaN is read by that comparison, not by a proof
  - a removable singularity exactly at the middle of a domain (the 0/0 of sin(x)/x over [-1, 1]) is now Failed, and the caller removes it; the rule could instead take a node beside the middle, and that is not done
---

The doubly exponential rule walks each side of its nodes until the values stop being finite, which is right where the integrand has overflowed at an end that weighs nothing, and wrong where the integrand is NaN in the middle of the domain: the side was cut there, the sum of what came before was reported, and `integrate` of a function that is NaN on half of [0, 1] answered 0.5008 with `DepthCapped`. The rule now looks at what the side had added just before the value that is not finite. If that was already negligible against the tolerance times the integral summed so far, the integrand had decayed and the overflow is at an end, as before (`x^2 exp(-x)` far out is `inf * 0`, and still converges to 2). If the side was still contributing, or had not contributed at all, the integrand is not finite where it matters and the result is `Failed` with NaN, as `sqrt(x)` over [-1, 1] and a function undefined past a point on a half line are now. The middle node is the same case: it was skipped when it was not finite, which left its weight out of every sum (`sin(x)/x` over [-1, 1] came back 1.8891 for 1.8922, untrusted), so it fails the integral too.

The finding was recorded when the door was written, as the fourth open item of the entry on quadrature over a domain that is a type; this closes it. It has a consumer on the day it lands: `integrate`, and through it `integrate_volume`, which weighs the status of each inner integral and had been reading a wrong number as a materially converged one.
