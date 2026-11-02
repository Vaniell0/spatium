---
date: 2026-11-03
title: a power of any number, not only a whole one from 1 to 64
files: [include/spatium/io/field.hpp, include/spatium/io/field_pod.hpp, include/spatium/io/field_glsl.hpp, include/spatium/io/field_text.hpp, include/spatium/io/field_parse.hpp, site/data/reference.json, site/wasm/questions.txt, tests/test_field_parse.cpp, tests/test_field_pod.cpp]
open:
  - synthesis does not search over Pow, so the root of x^p + p = 0 is answered as a Chebyshev polynomial and not as (-p)^(1/p)
  - a question whose bracket ends where the equation is infinite (0^-1 at x = 0) is refused as having no value; the refusal does not say why
  - a metric field refuses Pow
---

`x^p` was refused at the reader: "expected a whole power from 1 to 64 (x^a for another a is exp(a * log(x)))". A whole power typed as one, x^2 to x^64, is still a product, which is exact and keeps every answer that used it. Any other power is now `Op::Pow`: a name, a fraction, a negative, zero or a parenthesis (`x^p`, `x^0.5`, `x^-1`, `x^0`, `x^(p + 1)`).

`Op::Pow` is `std::pow` on the host and `pow` in GLSL. Where the exponent is not a whole number, it has a value for a > 0 only; the language's table says so, and its example `2^0.5` is answered by the library. In text a Pow whose exponent is a whole number from 1 to 64 keeps its parentheses, `x^(2)`, since a bare x^2 would be read back as the product. The round trip of 50,000 random fields through text, now including powers, still holds.

The every-op interpreter test includes Pow, and the tests that pinned the old refusals of `u^2.5`, `u^0` and `u^-1` now hold their powers. The WebAssembly core and the host agree on 40 of 40 questions, a power and a graph of `x^p` among them.
