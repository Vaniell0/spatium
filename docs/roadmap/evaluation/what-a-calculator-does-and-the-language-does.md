---
date: 2026-11-13
title: what a scientific calculator does, against what the language does, asked of the door
files: [scripts/calculator_inventory.py, site/wasm/core.cpp]
open:
  - complex numbers in the language (the library's Complex, its conjugate, abs and from_polar); the editor's ImaginaryUnit is refused by name
  - matrices and vectors as values (Matrix determinant and inverse, Vec dot and cross), once the matrices-and-systems fork is decided
  - factorial and binomial coefficients as functions of the language (discrete/combinatorics.hpp has both, overflow-checked)
  - inverse trigonometric and hyperbolic functions in the field and its table (Dual has asin and sinh; the field has neither)
  - angles in degrees as the editor writes them (Degrees); the library has the _deg literal only
  - physical constants and units as values (physics/mechanics/units.hpp, compile-time only)
  - regressions of data, statistics of a sample, distributions; a sum of a series; prime factorisation, gcd and lcm; random numbers; the vertex of a parabola; a table of values
---

`scripts/calculator_inventory.py` asks the library's own door -- the host build of the page's core, `site/wasm/core.cpp` -- each feature of a Casio fx-991EX as a person would type it: a line of the language, or the MathJSON row the page's formula editor writes for its notation (`5!` is `["Factorial", 5]`). A feature the door answers is in the language; one it refuses that the library declares, named by its header and the name, checked to be there, is in the library only; the rest are absent.

Of 45 features, 22 are in the language: fractions and decimals exactly, mixed numbers and percent through the editor, roots in simplest form, pi in exact results, variables, a value substituted, every function of the table, exact sines, a derivative at a point, a definite integral, every kind of equation the calculator solves (a root numerically, the quadratic with its complex roots, cubic and quartic, linear and non-linear systems, inequalities as sets, proportions). 12 are in the library only, refused at the door by name: degrees, inverse and hyperbolic functions, factorial and binomials, polar coordinates, complex arithmetic, matrices, vectors, regressions, physical constants and units. 11 are absent: display formats, bases and logic, degrees-minutes-seconds, random numbers, prime factorisation, a sum of a series, the vertex of a parabola, statistics of a sample, distributions, a spreadsheet, a table of values.

The order to close them is the open list, the library's first: each is plumbing what the library has through the language's table and the door, with no page work.
