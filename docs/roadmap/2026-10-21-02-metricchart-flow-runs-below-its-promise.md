---
date: 2026-10-21
title: MetricChart's flow runs at a hundredth of the accuracy it promises
files: [include/spatium/spaces/metric_chart.hpp, tests/test_metric_chart.cpp]
open:
  - the cost rises by about a third (T6 measured 102 -> 137 metric evaluations per exp on the sphere at 4e-13); a cheaper route needs a change in geodesic_rhs (several Dual seeds at once, diagonal or conformal metrics), not in the tolerance
  - the check covers the sphere against its closed form only; hyperbolic, torus and the skewed metric were checked against Real50 in a scratch program, not in the test suite
  - the five disagreements of rsc/spatium/find-report.md are the same after this change (only four error figures of TorusTyped@Dual lines moved); the TorusChart exp one is therefore not explained by the tolerance alone and stays open
---

`MetricChart` promised eps^(2/3) (4e-11 for double) and its flow, integrated at that same tolerance, missed it on 7-11% of inputs (worst 7.9x on the sphere). The flow now runs at tol/100 while the shooting keeps tol. Against Real50 on 150 held-out inputs per space (hyperbolic N=2 and 3, sphere, torus, skewed): 0 above 4e-11, worst 3.4e-12; before: 2-17 of 150 above, worst 3.14e-10. tol/10 was tried first and left one sphere input at 1.3x. `test_metric_chart.cpp` gets a test of 200 sphere inputs against the closed form that fails before the change (18 above) and passes after.
