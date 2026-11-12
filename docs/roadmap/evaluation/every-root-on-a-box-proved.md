---
date: 2026-11-04
title: every root of an equation on its box, each proved alone in its interval
files: [include/spatium/core/interval.hpp, include/spatium/algebra/roots.hpp, include/spatium/algebra/dual.hpp, include/spatium/io/field_pod.hpp, modules/algebra/roots.cppm, include/spatium/answer/evaluate.hpp, include/spatium/answer/respond.hpp, web/src/view/answer.ts, examples/gallery.cpp, site/gallery/gallery.json, site/wasm/questions.txt, tests/test_interval.cpp, tests/test_answer_evaluate.cpp, tests/test_answer_respond.cpp]
resolves: [evaluation/a-number-proved-to-its-digits#2]
open:
  - the ends of the doubles: where x^2 overflows (past a quarter of the largest double for 2x^2 - 3x - 5) a piece stays undecided; a chart at infinity, x = 1/t with `Series`, would close it
  - infinitely many roots (sin x = 0.5 with no bracket) exhaust the budget and prove none; breadth first spreads it over the whole line
  - a double root is undecided, never counted; its multiplicity is not sought
  - an expression read twice, x*x + y*y - 2*x*y, leaves a strip undecided that narrows like the square root of a piece's width; the mean-value form through the same dual number would narrow it like the width
  - roots are in intervals of doubles only; under Real50 they are not proved to fifty digits
  - the outcome is a refusal's text where no root is proved; a closed set of outcomes (roots, none, undecided where) is not a type yet
  - inequalities on a box, the other consumer of `prove_on`, are not asked from the question line
---

An equation's roots are now proved, with or without a bracket. `2^x = 8` has no closed form here; it is answered x = 3.000000000000000, to 16 digits, and the whole range of the doubles is proved to hold no other root, in 56 evaluations. `sin(x) = 0.5 on [0, 10]` gives its four roots to 14 and 15 digits, `x^3 - 2x - 5 = 0` gives 2.09455148154233, `x^2 + 1 = 0` is proved to have no real root in 4 evaluations. A question line takes `on [a, b]` as a document does.

`prove_on(box, statement)` (`core/interval.hpp`) asks a statement of a box of intervals and cuts the box across its widest side where the answer is undecided, breadth first, to a budget. A side is cut at zero, at the middle of its magnitudes when it spans many (the doubles' range reaches unit size in ten cuts), else at its middle. An answer whose computation took a comparison by guess is undecided whatever it said. The result is the partition, each piece with yes, no or undecided, and the verdict. `Undecided` now counts two things apart: a guessed comparison, which proves nothing (`guessed`), and a root or logarithm reaching below zero, whose interval still holds the values of the rest (`outside`).

`roots_on(f, box)` (`algebra/roots.hpp`) asks of each piece whether every root on it is accounted for. Yes where f's interval leaves out zero. Yes where Newton's interval step, m - f(m) / f'(X) with f' enclosed by `Dual<Interval>`, misses X or lands inside it: then f' keeps one sign, and X holds exactly one root (Moore). A root on a cut sits at the end of two pieces and inside neither, so a step that contracts but sticks out is retried on itself widened (Rump's epsilon-inflation). A proved root is narrowed by the same step until it stops shrinking. The lowered field runs in `Dual<Interval<double>>` unchanged: the interpreter takes a dual number of enclosures as an enclosure, and `pow(Dual, Dual)` was added for 2^x.

`evaluate` keeps it as a route, "every root, proved", after the closed form and before bisection. Where the closed form answers, the proof is the second path, and the page shows its digits beside the formula. `respond` writes `enclosed: {box, roots: [{digits, text, lo, hi}], undecided, complete}` and the roots as LaTeX. The page says "1 root proved, and no other on the box", or what is still undecided. The gallery has two new cards.

The cases it must get right, and did. A jump is no root: `less(s, 0.5) - 0.5 = 0` changes sign at 0.5, and bisection returns 0.5, while the proof leaves a piece 1e-8 wide undecided. A corner is no root: abs's derivative guesses a branch at its kink, so a piece is first asked of f alone, and `abs(x - 2) = 1` is complete with two roots. `x^4 = 0` is undecided within 2.3e-10 of zero. `sqrt(x) = 1` on [-1, 4] is one root: the part of a piece below zero has no value and no root. A piece of one magnitude where f overflows a double at both ends is left uncut; before that rule, a quadratic with no bracket spent 8156 evaluations at the end of the doubles, after it 175.

The pairs: 106 quadratics, each root held against the formula's double and against sixty digits, at least 14 digits proved. The four roots of sin x = 1/2 against π/6 + 2πk and 5π/6 + 2πk. (x - y)^2 >= 0 written out is undecided only on pieces the diagonal crosses: undecided area 1.47, 1.06, 0.73 of 4 at 250, 1000, 4000 evaluations. WebAssembly and host agree on 56 of 56 questions.
