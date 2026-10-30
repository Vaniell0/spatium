---
date: 2026-10-30
title: where no short expression answers a piece, a Chebyshev polynomial does, and says so
files: [rsc/spatium/synthesis.hpp, rsc/spatium/query.hpp, rsc/spatium/respond.hpp, examples/gallery.cpp, site/gallery/gallery.json, tests/test_rsc_synthesis.cpp, tests/test_rsc_query.cpp]
resolves: [2026-10-30-12#1]
---

A function with no closed form in the operations came back in many pieces of a few short terms: the sine integral on [0.5, 6] in 59, the error function on [0.25, 3] in 46, the Fresnel C on [0.25, 2] in 93. `chebyshev(reference, spec)` is the fallback: the polynomial that interpolates the reference at the Chebyshev points of the piece, of degree 4, 8, 16, up to `spec.polynomial_degree`, the first the check grid accepts. It is a field like any other -- t = (x - centre) / half summed by Clenshaw's recurrence, about four ops a degree -- so it prints, lowers to plain data with the interpreter's bits, and is certified on the same 2001 points as an expression; it is not a closed form, and the piece says so (`Form::Chebyshev`, its coefficients beside it).

`synthesize_piecewise` tries it on a piece no short expression answers, before halving, when the degree is set; a question asks for it (`Query::polynomial_degree`, 32), the measurements of the search itself do not. Each of the three special functions, and the integral of exp(-p x^2), is now one piece (Si: degree 16, error below 1e-12), and the gallery is written in 0.76 s instead of 17. The case it must fail: |x - 1| on [0, 2] is no polynomial of degree 32 to 1e-10. The gallery card keeps a Chebyshev piece's form and degree but not its coefficients, which carry the last digits of the libm that computed the reference; a live answer has them.
