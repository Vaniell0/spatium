---
date: 2026-11-16
title: a function of one argument, f(t) = 726t + 2815, called on the lines below it as f(-5)
files: [include/spatium/io/field_parse.hpp, include/spatium/io/document.hpp, tests/test_document.cpp, tests/test_io_roundtrip.cpp]
open:
  - one argument only; f(x, y) is refused as a call with two arguments
  - a formula editor's row for a definition, f(t) = ..., has no MathJSON pattern yet, so the page reaches it as text only
  - a name an equation above solves (x: 2x - 4 = 0) is still not readable below as its value
---

A school mathematics benchmark (DeepMind's Mathematics Dataset, 54 modules) asks "Let c(t) = 726*t + 2815. Calculate c(-5)." and builds most of its composed questions the same way; the language had no form for it, and those questions could not be written at all. A line `f(t) = body` is now a value with its argument recorded (`Statement::variable`), written, printed and carried through JSON by one more row of `statement_forms`; the argument is a new name, as an integral's variable is. A call `f(a)` on a line below is read by `FieldReader` as the body with the argument put in, through a name no text can hold, so `f(t + 1)` with f's argument t is not read twice; `read_document` hands each line the functions above it. The line keeps its text as written, `answer = c(-5)`. On the benchmark the module that evaluates a defined function went from none to all 30 questions, and the total from 244 to 302 of 1620.
