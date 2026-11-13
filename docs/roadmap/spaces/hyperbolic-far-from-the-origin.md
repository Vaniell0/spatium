---
date: 2026-11-13
title: Hyperbolic far from the origin, exp and log without cancellation
files:
  - include/spatium/spaces/hyperbolic.hpp
  - tests/test_hyperbolic.cpp
open:
  - the hyperboloid model still rebuilds a point near the origin from one d away as q cosh d plus a tangent of the same size, so about e^{2d} of rounding stays (1.9e-10 at d = 8, 9e-7 at d = 12); a model whose coordinates do not grow with d (the Poincare ball of ConstantCurvature) is where a long way would be computed instead
---

Enumerating laws over the operations of `core/access.hpp` on the spaces of the connectivity matrix -- terms in log, exp, the metric and the distance, kept where they agree on every sample, with inputs moved to break what holds by chance -- found exp_q(log_q p) = p true on every exact space at the sampler's scale and false on every hyperbolic one once q was 8 away. Checked against Real50, both halves lost digits as e^d grew.

A tangent at q has coordinates of size d e^d and length d. `exp_map` read its length as -v0^2 + |v_s|^2, where terms of size e^{2d} cancel to d^2; from the exact tangent the error of exp was 4e-11 at d = 4, 2.6e-4 at d = 8, 1.4e3 at d = 12 and 5e10 at d = 16, with no refusal. It now reads the same number for a tangent at p as (|v_s|^2 + sum over i < j of (v_i p_j - v_j p_i)^2) / p0^2, Lagrange's identity turning the difference into a sum of squares: 2.5e-15, 9.6e-14, 1.9e-10, 9e-7, 3e-3 at the same distances, which is the e^{2d} the model cannot avoid. `log_map` divided by |proj|_L, sinh d, read off a projection that cancels the same way (5.5e-11 at d = 8, 2.1e-3 at d = 16); past d = 1 it now takes sinh d from the inner product as sqrt((x - 1)(x + 1)), and keeps its relative error at 1e-16 to 2e-16 at every distance tried. Below d = 1 the projection's own length is kept, which is the better one where the points are close. Reading sinh d off the inner product holds only for points on the sheet, and an iteration of exp steps had let them drift: the Frechet mean at H^8 ended 1e-6 off the sheet after fifty steps and stalled 5e-4 from the mean, which the old log, normalising by the projection itself, had absorbed. `exp_map` now puts its result's time coordinate back as `project()` does, sqrt(1 + |x_s|^2), where nothing cancels, so no step leaves the sheet.

The test holds a point 2, 4, 8 and 12 away from one near the origin in H^3: the distance and the length of log at the far point to 1e-12 relative, and exp after log to 64 eps cosh^2 d, which the old code missed at 8 by four orders.
