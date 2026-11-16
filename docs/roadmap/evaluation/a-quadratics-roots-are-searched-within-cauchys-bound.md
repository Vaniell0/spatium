---
date: 2026-11-16
title: a quadratic's roots are searched within Cauchy's bound, so its check proves there is no other and says nothing about overflow
files: [include/spatium/answer/evaluate.hpp, tests/test_answer_respond.cpp, tests/test_answer_evaluate.cpp, site/gallery/gallery.json, site/data/reference.json, docs/language.md]
open:
  - an equation no closed form here reads (a cubic, x^3 - x = 0) is still searched to the ends of the doubles and leaves those ends undecided; a polynomial of any degree has the same bound once its coefficients are read
---

Every quadratic's Check carried a line "Undecided [4.494e+307, 1.798e+308], where the equation overflows a double": the second route, every real root proved, searched the doubles' whole range, and near its ends the equation's value overflows. The closed form has already read the equation as a polynomial of degree at most two, so every root lies within Cauchy's bound, 1 + max |a_i / a_n|. The search is now held to twice that and the bound is a step of its own. The check is complete, so the verdict gains its strongest note: "Exact, closed form, checked by Vieta's formulas, no other real root". The gallery and the reference carry the new verdict.
