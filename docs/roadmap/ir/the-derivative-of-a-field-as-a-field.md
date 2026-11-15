---
date: 2026-11-13
title: the derivative of a field as a field, checked against a dual number
files: [include/spatium/io/field_derivative.hpp, include/spatium/io/field.hpp, include/spatium/io/field_parse.hpp, include/spatium/io/document.hpp, include/spatium/answer/query.hpp, include/spatium/answer/evaluate.hpp, include/spatium/answer/respond.hpp, examples/site_reference.cpp, site/wasm/questions.txt, web/src/registry.ts, tests/test_field_derivative.cpp, tests/test_io_roundtrip.cpp]
open:
  - the derivative is written as its rules leave it, d/dx x^3 as (x + x) * x + x^2; no rewrite brings it to its shortest form yet
  - noise and a net that read x are refused; a net's derivative is a sum of tanh' terms and could be written
  - a derivative of a derivative is two lines, h = d/dx g; there is no d^2/dx^2
  - the variable of a derivative that is no parameter is checked and drawn on [-2, 2], whatever the function's domain
---

`derivative(f, "x")` (`io/field_derivative.hpp`) differentiates a field in a symbol and writes the result as another field: an op at a time, in order, by the chain rule. The rules of + - * /, of a power (b a^(b-1) a' + a^b log(a) b') and of the comparisons are in code; a function of the language carries its own in its row of the table, `Function::derivative`, written in the language over its arguments a and b -- sin's `cos(a)`, min's `less(a, b)` and `1 - less(a, b)` -- so the reference page shows the same text the rule is read from. An op that reads nothing that reads x adds no term, and a factor of exactly 1 is left out, so d/dx (2x) is 2 and not 0 * x + 2 * 1. Constants are left as written: x^3 / 3 at 1 is exactly 1 in fractions, where folding 1/3 first gave 9999999999999999/10^16.

In a document it is a statement, `g = d/dx f`, a function of x; `g = d/dx f at 2`, a number; and, when x is a parameter, the number at the value chosen for it. The calculator takes `d/dx f` alone, and a formula editor's `\frac{d}{dx} f` and `\frac{d}{dx} f, x = 2` as MathJSON (`D`, and `D` in a `Sequence` with the point). The derivative is taken after the values the function reads are put in and before a chosen value is, which would make its variable a number. A row of the editor that holds to one pattern and departs from it no longer stops the search: a root's row and a derivative's both start as a `Sequence`.

The answer goes by its own route, `derivative`, and is checked by a second one: the derivative a dual number computes from the function by its own rules (`algebra/dual.hpp`), both in intervals, at sixteen points of the variable's interval, or at the point. Where the two do not overlap the route refuses, since one of the two sets of rules is wrong. Every row of the table, ten composite functions and the check refusing two wrong derivatives are held in `tests/test_field_derivative.cpp`; the function, the point and the chosen parameter give one number; the derivative's LaTeX row read back by the page's own parser answers as the typed line.
