---
date: 2026-10-30
title: a field rewritten without changing a value, and a bound name folded to a constant
files: [include/spatium/io/field.hpp, include/spatium/io/build.hpp, tests/test_field_simplify.cpp]
open:
  - a coefficient that is a small fraction (0.16666666666666663 for 1/6) is a different number from the fraction, so writing it as one is a step of the synthesis with its own certificate, not a rewrite
  - lower() does not simplify: a constant folded on the host in double would reach a device that computed it in float, and the two differ
  - no rewrite that reads better but costs an op (x + x as 2 * x)
---

`Field::simplified()` makes only the rewrites that give the same bits on every input: an op whose children are all constants becomes the constant the evaluator itself computes for it, x * 1, 1 * x and x / 1 become x, and x - 0 and x + (-0) become x. One pass in order, so 2 * 2 * u folds to 4 * u. Not made, because each changes a value somewhere: x * 0 (NaN, infinity), x - x and x / x (infinity), x + 0 (-0 + 0 is +0), reassociation. Held on +0, -0, the infinities and NaN as well as ordinary inputs, on hand-made fields and on 300 random ones, which go from 1799 ops to 1402.

The consumer is the scene: a field handed to a builder has its names bound to numbers and is then simplified, so a thickness written 2 * r * r with r = 0.5 is one constant in the node, while the node keeps it as written for the description.
