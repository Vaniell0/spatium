---
date: 2026-11-13
title: an answer in its shortest form, kept only where its values hold
files: [include/spatium/io/field_shortest.hpp, include/spatium/answer/evaluate.hpp, include/spatium/answer/respond.hpp, tests/test_field_shortest.cpp]
open:
  - the text form of a leading negative coefficient keeps its parentheses, (-2) * x; the page's LaTeX does not
  - no factoring: x^2 - 1 is not written (x - 1)(x + 1), and a common factor is not taken out of a sum
  - a coefficient that was a decimal may come back a fraction, 0.1 * x as x / 10, when the two show as much
  - log, exp and powers of a name are not rewritten among themselves (log(exp(x)), exp(a) exp(b), x^a x^b)
---

`shortest(f)` (`io/field_shortest.hpp`) writes a field as briefly as it can and keeps the rewrite only where a page shows no more of it and its values hold. It reads f as a sum of terms, each an exact fraction times a product of atoms to whole powers -- a name, pi, the root of a squarefree whole number, a function or a power with its arguments written shortest in turn, a denominator that is no single term -- so like terms are collected, terms share one denominator, a root comes out of its square factors with its square a number, a number to a fraction is its root where that is a fraction, and a sine or cosine of a multiple of pi over 1, 2, 3, 4 or 6 is its exact value. What a page shows is counted from the LaTeX, a command as the marks it draws (`shown_length`): (x + 1)^2 stays, since x^2 + 2x + 1 shows more.

A rewrite is safe by its check, not its algebra: `same_values` compares the two at fourteen points, every name given a value, where both have one and where neither has, then on a small box about each point in intervals. So x / x is not 1, which has a value at 0.

What an answer shows goes through it: a derivative, (x + x) * x + x^2 as 3x^2 and (1 / x - log(x) / x) / x as (1 - log(x)) / x^2; a function of a parameter; and a number that is no fraction beside its proved digits, sqrt(24) + sqrt(150) = 7 sqrt(6) ~ 17.146. The values are still computed from the expression as written. In passing, a value that reads pi is no longer called exact: the double nearest pi read as a fraction made 2 pi "exactly" 3141592653589793/500000000000000; it is now 2 pi, proved to 15 digits.
