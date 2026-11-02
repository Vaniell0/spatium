---
date: 2026-10-30
title: a field read from text, and a question typed as one
files: [include/spatium/io/field_parse.hpp, include/spatium/io/field_text.hpp, rsc/spatium/query.hpp, tests/test_field_parse.cpp, tests/test_rsc_query.cpp]
open:
  - the language is the field's own operations; no tan, no pi, no power with an exponent that is not a whole number from 1 to 64 (x^a is written exp(a * log(x)))
  - no implicit product: 2x is refused, 2 * x is read
  - the inputs a field evaluated per instance or per point reads (id, origin, point, hash, the coordinates of spacetime) and the noise and gather ops are printed but not read
  - a question names exactly one unknown and one parameter
---

`parse(text, u, v)` reads a field written in the usual notation into the ops `to_text` writes it from: sums and products to the left, a unary minus tighter than a product and looser than a power, x^n as n factors, the names given for the inputs as U and V and every other name as a symbol, `inf` and `nan` as numbers. A text that is not a field is refused at its column with what was expected there (`parse: at column 4, expected a number, a name, a function or '(', found the end`).

The two are held against each other on 400 random fields over the inputs, a symbol, special constants and every operation of the language: the text reads back to the same bits on twelve inputs each way (zeros of both signs, infinity and NaN among them) and prints to the same text again; on the fields where the printer shortens nothing (no factor of exactly 1, no negative constant added) it is the same ops by canonical key. This found that the plain text was not yet a program: u + (v + 1) was printed u + v + 1, which reads as (u + v) + 1, and floating-point addition does not reassociate (at u = 1e300, v = -1e300 the two give 0 and 1). The text now keeps the parentheses of a right operand of a sum or a product; LaTeX, read as mathematics, still drops them. And a sum whose right operand began with a minus was printed as a subtraction whatever it was (u + (-2 - v) as u - 2 - v); only a negative constant is now.

The consumer is the question: `rsc::ask("solve x^2 - p = 0 for x in [0, 10], p in [0.5, 4]")` is the `Query` with those names and brackets, whose answer is the answer to the question built in code, to the last bit; `integrate f for x in [a, b], p in [c, d]` likewise. A bound may be an expression in no name (`sqrt(2) * 10`). A question that is not one is refused where it stops being one, an empty bracket and a name used twice by name.
