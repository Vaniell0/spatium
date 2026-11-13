---
date: 2026-11-04
title: a number shows how it was computed, and is exact where it is a fraction
files: [include/spatium/io/steps_text.hpp, include/spatium/io/field_mathjson.hpp, include/spatium/core/rational.hpp, include/spatium/answer/evaluate.hpp, include/spatium/answer/respond.hpp, web/src/answers.ts, web/src/page.ts, web/src/view/answer.ts, web/test/notation.test.ts, examples/gallery.cpp, site/gallery/gallery.json, site/wasm/questions.txt, tests/test_solution_steps.cpp, tests/test_answer_evaluate.cpp, tests/test_answer_respond.cpp]
open:
  - only an evaluated number has steps from its computation; a closed form's are still written by hand (`to_steps`), and the proof of roots has none
  - a program with one operation fractions do not have (a root, a sine) is in intervals throughout: 1/3 beside sqrt(2) is written 0.333333333333333, not 1/3
  - a power a whole number long is a chain of products in the IR; its step is shown as one, the chain itself is not a power
  - the answer as a page shows it is still the route's own block; one sheet for every answer is the next step
---

An evaluated number now says how it was computed, an operation a step: `2(3 + 4)^2 - 1/4` is 3 + 4 = 7, 7^2 = 49, 2 · 49 = 98, 98 − 1/4 = 391/4 = 97.75. The steps are read off the lowered program as it ran (`evaluation_steps`, `io/steps_text.hpp`): each operation written on the values of what it reads, a number or a name as it was written, a power the reader made into a chain of products as one step, a fraction as written no step of its own. The labels are the language's own, from the table the reference shows (square root, sine, absolute value). Nothing is worked out again from the answer.

The values are fractions where every operation is one on fractions, in 64 bits, and the digits an interval of the same program proves elsewhere, written with ≈. Where no digit is proved the step says how far from zero the value may be: 10^20 + 1 − 10^20 is 1, past 64 bits, and a double's interval holds it only within 10^5 of zero, which the step writes as 0 ± 10^5.

The answer is exact where the computation is: the evaluation runs the program in fractions before intervals, and a value that is one is the value (`Evaluation::exact`). 0.1 + 0.2 is 3/10 = 0.3, (3 + 7/12 − (2 + 11/18)) / (19/84) is 245/57, and the exercise below is 12, not 12.000000000000002 with fifteen proved digits. `Rational::decimal` writes a fraction whose decimal ends, in integers.

The exercise the gallery now carries is a page of problems' own, ((3,2 − 1,7) : 0,003 / ((29/35 − 3/7) · 4 : 0,2) − (1 13/20 − 1,5) · 1,5 / ((2,44 + 1 14/25) · 1/8)) : 62 1/20 + 1,364 : 0,124. Typed into the formula editor it found three gaps, all in reading: the editor's parser reads a colon as a relation looser than a sum, so the sum was divided; it reads 3{,}2 as 32; and it writes a mixed number 1 13/20 as the bare juxtaposition a product is. The page now writes a colon as \div and a decimal comma as a point before parsing (`notation`), and the core reads a whole number before a fraction of whole numbers as their sum. Its 18 steps are exact, 1 13/20 = 33/20 among them, and it is 12.

The pairs: a chain's last value is the fraction the rational engine computes, on four expressions; the exercise as the editor's parser writes it is a test, a question of the WebAssembly check (51 of 51 agree with the host) and a card that reads back from its rows (18 of 18).
