---
date: 2026-10-30
title: the quadratic solver keeps the small root, in the form that does not cancel
files: [include/spatium/algebra/polynomial.hpp, tests/test_solution_steps.cpp, docs/consumers-baseline.json]
resolves: [2026-10-30-11#1]
---

`solve_quadratic` computed both real roots as (-b ± sqrt(D)) / 2a, and the one whose numerator subtracts two numbers of nearly one size lost its digits: for x^2 + 1e8 x + 1 the small root came out -7.45e-9 against -1e-8, which Vieta's formulas showed (the product of the roots 0.745 instead of 1). It now uses q = -(b + sign(b) sqrt(D)) / 2, which adds two numbers of one sign, and takes the roots as q / a and c / q, in the same order as before; a double root is still -b / 2a for both, so the two are one number to the bit. The small root of x^2 + 1e8 x + 1 is -1e-8 to 1e-22, and Vieta's formulas hold to 1e-12 on it and on 20000 integer quadratics; roots off by the old loss are the case the check must fail. Every test that solves a quadratic -- the ray-quadric hits and misses among them -- passes unchanged.

The consumer baseline is brought down to the names unused today (163, from 168), so the five that gained a caller this week cannot lose it again unnoticed.
