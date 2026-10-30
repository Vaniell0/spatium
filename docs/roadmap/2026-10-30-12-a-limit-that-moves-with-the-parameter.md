---
date: 2026-10-30
title: a limit that moves with the parameter, and the special functions it writes as questions
files: [rsc/spatium/query.hpp, include/spatium/io/field_parse.hpp, include/spatium/io/field_text.hpp, examples/gallery.cpp, site/gallery/gallery.json, tests/test_rsc_query.cpp]
resolves: [2026-10-30-05#1]
open:
  - a function with no closed form in the operations comes back in many pieces of a few short terms (Si on [0.5, 6] in 59, erf on [0.25, 3] in 46, the Fresnel C on [0.25, 2] in 93); a polynomial of higher degree per piece would take a few
  - an integral over an infinite interval is not a question yet (Ei, the exponential integral)
---

A bound of the unknown may now read the parameter, so a question can have its parameter as a limit of the integral: `integrate sin(t)/t for t in [0, x], x in [0.5, 6]` is the sine integral Si(x), which the operations do not have. The interval must hold something wherever the answer is wanted, checked at 17 points; a bound may read only the parameter, and the parameter's own bounds are numbers. `pi` is a number in the language, the double nearest pi, and is printed as its name, which reads back to the same bits.

The reference of each such question is held against a second route: Si against its power series to 1e-13 on [0.5, 6], erf (written 2 / sqrt(pi) times the integral of exp(-t^2)) against `std::erf` to 1e-15, the Fresnel C against its series to 1e-13. The gallery has the three as cards.
