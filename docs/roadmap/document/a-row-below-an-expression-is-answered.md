---
date: 2026-11-04
title: a row below an expression is answered
files: [include/spatium/answer/respond.hpp, tests/test_answer_respond.cpp, site/wasm/questions.txt]
open:
  - a root's value is not read by a row below it: x^2 = 2 and then 2x + 1 is refused (the expression reads x, which evaluate takes for neither a parameter nor a value)
---

On a page of several rows every row after an unnamed expression was refused, "line 1: parse: at column 1, expected a name": the rows above the last were read as lines of a document, and an expression typed alone is no line of one. Each row above the last is now a line as written, or else the same shorthand the last row is -- an expression a value named answer_k by its row, which a row below may read (answer_1 + 1), an equation a root -- and a row that is neither is refused by its row. The last row is read against the document the rows above it made. The pairs: five pages of rows, two refused by the row they are on (`test_answer_respond`); a page of two rows is a question of the WebAssembly check, which agrees with the host on 52 of 52.
