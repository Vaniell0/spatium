---
date: 2026-10-30
title: an equation solved in closed form with its steps, and the gallery written by the library
files: [include/spatium/algebra/polynomial.hpp, include/spatium/io/steps_text.hpp, rsc/spatium/closed_form.hpp, examples/gallery.cpp, site/gallery/gallery.json, tests/test_solution_steps.cpp]
resolves: [2026-10-30-05#2]
open:
  - steps only for the quadratic; the cubic and the quartic solvers keep no record yet
  - a question with a parameter shows its answer and certificate, not the route it was found by as steps
  - the bounds of a bracket are not yet labelled in the steps
---

`solve_quadratic_steps(a, b, c)` is the quadratic solver with what it found the roots through: the coefficients, the discriminant and which case the equation is in (every x, no x, linear, two real roots, a double root, a complex pair). `solve_quadratic` returns its roots, so the steps and the answer cannot differ in a bit. `to_steps` writes the record as a chain a page reads, each step a label and a line of LaTeX, the formula and then the formula with the numbers in it: `D = b^{2} - 4ac = \left(-3\right)^{2} - 4 \cdot 2 \cdot \left(-5\right) = 49`. Held on 20000 integer quadratics, every case reached, against Vieta's formulas in complex arithmetic.

`rsc::solve_closed_form("2x^2 - 3x - 5 = 0")` reads an equation in one unknown and takes its coefficients from the field's values at 0, 1 and -1; that it is of degree at most two is checked at five more points, so `(y - 1)(y + 2) = 4 - y` is y^2 + 2y - 6 and `x^3 - 1 = 0` is refused by name. No expression is rewritten. The parser reads a product after a number or a closing parenthesis without its '*' (2x^2, (x - 1)(x + 2)); after a name it does not, so u v is still refused.

The consumer is the gallery: `examples/gallery.cpp` writes `site/gallery/gallery.json` -- two equations solved in closed form with their steps, four questions with a parameter answered by synthesis or refused, the prism's volume, a scene in LaTeX -- and its `--check` is a ctest, so the file cannot drift from the library. An error is written as the power of ten above it, so the file does not move with the last digits of a libm.
