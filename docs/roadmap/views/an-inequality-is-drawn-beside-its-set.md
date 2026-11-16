---
date: 2026-11-16
title: an inequality is drawn beside its set
files: [include/spatium/answer/respond.hpp, tests/test_answer_evaluate.cpp]
open:
  - an inequality in two unknowns has no region drawn; the language answers one unknown
---

An inequality's answer was its set on the axis and no graph: `graph_of` drew only an equation's roots. It now draws an inequality the same way, its left side less its right along the unknown, over the roots and breaks the interval method proved on its axis with a margin on each side; with none, over its bracket, or [-5, 5] on the whole line. `x^2 - 4 >= 0` is drawn on [-4, 4], `1/x > 0` on [-1, 1] about its break, `x^2 + 1 > 0` on [-5, 5]; each request draws (`test_answer_evaluate`).
