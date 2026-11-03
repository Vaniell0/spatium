---
date: 2026-11-03
title: the rows of a LaTeX editor are a document and its question
files: [include/spatium/io/field_mathjson.hpp, rsc/spatium/respond.hpp, site/wasm/core.cpp, tests/test_field_mathjson.cpp, tests/test_rsc_respond.cpp, site/wasm/questions.txt]
open:
  - the page does not send rows yet: the question as rows of MathLive fields is the page's step
  - an equation with no interval below a parameter, x^2 - p = 0 under p in [0.5, 4], is refused by synthesis (x has no bracket); its roots as functions of p, ±√p, are a closed form the solver does not give yet
  - a root on its bracket written as a line, x: x^2 - 2 = 0 on [0, 2], is answered in closed form with both roots, the bracket not narrowing them
  - an integral with one end open, or with no ends, is refused; the language's integral has both
---

A page asks a question as rows, one formula a row, the way a person writes a problem down: the parameter, then the equation and where its unknown is.

```
p \in [0.5, 4]
x^2 - p = 0,\ x \in [0, 10]
```

The rows come from a LaTeX editor whose parser (the Compute Engine's `latex-syntax`) writes MathJSON before any rewriting, and that form is read too: a product written without a sign (`InvisibleOperator`, `2x`), `Lb` and `Lg` for the logarithms in base 2 and 10, and `e` as the number e, so `e^{x}` and `\exp(x)` read alike as `exp(x)`. A row is one line of a document (`mathjson_row_text`): `p \in [0.5, 4]` the parameter `p in [0.5, 4]`; an equation and its unknown's interval, which the editor writes as a list, the root `x: x^2 - p = 0 on [0, 10]`; `I = \int_0^1 ... dt` the integral `I = integral ... dt on [0, 1]`; an equation or an expression as written.

What the lines mean is the language's. `respond` now takes several lines: the lines above the last are a document, and the last is the question -- the shorthand a person types on one line (`2 + 2`, an equation in its one unknown), read with the names above it, or else the document's last statement. A single line that is a statement and no shorthand, `x: x^2 - 2 = 0 on [0, 2]`, is answered as that document; a shorthand answered before is answered as before, and the gallery and the reference did not change. `respond_rows` (`spatium_rows` in the WebAssembly core) reads the rows into those lines and answers them, and a row that does not read is refused with its number: `row 2: at [2]: 'D' is not in the language`.

Two routes are held equal: the rows above and the same two lines typed give the same answer, to the byte (√p, synthesised and certified), and `2x^2-3x-5=0` from the editor and typed give the same steps. Each row form is taken from the parser's output, measured in a browser, and the rows as a document are the document written by hand, field by field. The WebAssembly core and the host agree on 51 of 51 questions, five of them rows: √p, the roots 5/2 and -1, 3a = 12 with a = 2 + 2, an integral of a parameter synthesised as a Chebyshev series, and a refusal.
