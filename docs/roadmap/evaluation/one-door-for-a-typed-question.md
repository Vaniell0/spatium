---
date: 2026-10-30
title: one door for a question typed as text, answered as JSON, with what its certificate is worth
files: [rsc/spatium/respond.hpp, rsc/spatium/synthesis.hpp, examples/gallery.cpp, site/gallery/gallery.json, tests/test_rsc_respond.cpp]
open:
  - the libm a certificate was computed with is not recorded
  - a power is printed as the product it is read into: 2(3 + 4)^2 reads back as 2 * (3 + 4) * (3 + 4)
  - a value is shown in double only; 1/3 is 0.3333333333333333, not also the fraction
---

`rsc::respond(text)` takes what a person types and sends it to the route it belongs to: an expression with no unknown is evaluated, as a calculator does (`2 + 2` is 4, `2 + 2 =` too); an equation in one unknown is solved in closed form with its steps; a question with a parameter (`solve ... for x in [a, b], p in [c, d]`, `integrate ...`) is answered by synthesis. A text that is none of these is refused with the library's own reason -- a parse error with its column, a cubic by name, a name with no value with the equation to write instead. The answer is the JSON a page reads, so the page composes nothing; the gallery is written by the same door, so a card and a live answer cannot differ.

What a synthesised answer's certificate is worth goes with it (`Certificate` on every `Solution`): the strength (sampled -- the error measured at the points of a Chebyshev grid, not bounded between them), the norm (relative or absolute), the number of points, the scalar, and whether it was computed on the host or in WebAssembly. The error is written as the power of ten above it; the test holds that bound above every piece's error and within a factor of ten of the largest.
