---
date: 2026-10-29
title: named symbols in the field, a question asked in them, and a field as LaTeX
files: [include/spatium/io/field.hpp, include/spatium/io/field_pod.hpp, include/spatium/io/field_glsl.hpp, include/spatium/io/field_text.hpp, rsc/spatium/query.hpp, tests/test_field_symbols.cpp, tests/test_rsc_query.cpp]
resolves: [2026-10-29-11#2]
open:
  - an unbound symbol does not lower; a slot the device reads (a uniform, a slider on a page) is the next step, with the scene that needs it
  - the scene DSL does not take symbols yet: a parameter of a scene given by a number or by an equation solved at build
  - no simplification: a symbol bound to a number leaves 2 * 2 in the field, which a shader compiler folds and the host interpreter does not
  - the vector field has no symbols; a motion reads one only through its scalar factors
---

`Field::symbol("r")` is a named number the field does not fix, the `Symbol` op appended after `Exp` so saved op codes do not move. A field with one is a family of fields: `free_symbols()` names them (derived from the ops the root reads, like `is_structural()`), and `bind(name, value)` picks a member, a number or another field, rebuilt by appending so the order holds and equal ops merge. Bound, it is the field written out: the sphere x^2 + u^2 - r^2 with r bound to 2 and x to v has the canonical key, the content id and the values, bit for bit, of v^2 + u^2 - 2 * 2. Unbound, it evaluates to NaN and does not lower, and the error names the symbol.

A question can be asked in names: `x * x - p` with `unknown = "x"` and `parameter = "p"` is bound to U and V before anything is evaluated, and its answer is the answer of the question written in U and V, to the last bit. A symbol the question reads but does not name is refused by name.

`to_text` now walks the field itself rather than its lowered form, since the lowered form keeps no names, and puts parentheses where precedence needs them: `u / (v * 2)` used to print as `u / v * 2`, which reads as (u / v) * 2. `to_latex` is the same walk spelled for a page: `\frac`, `\sqrt`, `\cdot`, an Iverson bracket for `less`, `x^{k}` for a coordinate of spacetime, Greek names by their command and `r_0` with a subscript, so the answer to the root of x^2 - p reads `\sqrt{p}` in code and on a page from one source.
