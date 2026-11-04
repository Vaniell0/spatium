---
date: 2026-11-04
title: a number's digits proved over the reals, by intervals
files: [include/spatium/core/interval.hpp, modules/core/interval.cppm, include/spatium/io/field_pod.hpp, rsc/spatium/evaluate.hpp, rsc/spatium/respond.hpp, web/src/view/answer.ts, tests/test_interval.cpp, tests/test_rsc_evaluate.cpp, site/gallery/gallery.json]
open:
  - roots, integrals and synthesis carry no proof yet: a root over the reals without a bracket (interval Newton on subintervals), a synthesis proved on its pieces
  - a box of inputs is not split when a comparison is undecided; the value is then only refused a proof
  - the connectivity matrix has no interval column (the level "true on the whole box")
  - `Certified<T>` (kind, strength, digits) is not a type yet; the proof rides on the evaluation
  - an end of an interval at infinity is a number of T, not a limit (`Series` could carry it)
  - the libm bounds are taken from glibc's and musl's documentation, Boost's from a test against twice the digits
---

Every number the calculator evaluates now says how many of its digits are proved. `(3 + 7/12 - (2 + 11/18)) / (19/84)` is 4.2982456140351 to 14 digits, `0.1 + 0.2` is 0.300000000000000 to 15 (the double says 0.30000000000000004), `2 + 2` is exactly 4.

`Interval<T>` (`core/interval.hpp`) is a Scalar over a double or `Real<N>`. Each operation rounds its lower end down and its upper end up, unless an fma or a two-sum shows it was exact, so the value of what was written, over the reals, stays inside. A constant enters by its text: 0.1 is between the two doubles beside it, 0.5 is itself. An integer power of one number is taken by squaring, so 10^20 stays one number, and an even power keeps its sign ([-1, 2]^2 is [0, 4]). The digits proved are the most for which both ends round to the same decimal.

Comparisons are three-valued (`less`, `equal`: yes, no, undecided). The operators generic code calls take the midpoints' order on an undecided pair and raise a flag, as IEEE raises inexact, and `Undecided` reads it. So every Scalar algorithm runs on intervals unchanged, and a result that followed a guessed branch is not called proved. A square root or logarithm reaching below zero raises it too: `sqrt(0.1*3 - 0.3)` is 7.45e-9 in a double, and no proof, since over the reals it is the root of zero.

`evaluate` runs the same lowered program in `Interval<double>` (`Interval<Real50>` under the Real50 engine), and `respond` writes `proved: {digits, text, exact}`, or `{digits: 0, bound}` when the interval holds zero. The page names it: "Proved to 14 digits", "Exact".

The case a double must fail: `sin(3.141592653589793)` is 1.2246e-16 in a double, while the decimal written is short of pi by 2.3846e-16. The interval proves no digit and bounds it within 1e-15 of zero; under Real50 it proves 2.384626433832795028841971693...e-16. That test found a wrong bound on the way: Boost's sine near its zero is off by epsilons of the argument, not of the result, so fifty digits of the argument leave about thirty-three of the sine. Functions over `Real<N>` now widen by the argument's size.

The pairs: 1143 decimals in [-4, 4] through a formula with every operation, `Interval<double>` against `Real<60>`, all held, and the midpoint never equal to the truth. `Interval<Real<50>>` is held against `Real<120>` on 199 cases, the cuts of pi among them. Proved digits agree with Real50's to half a unit of the last. WebAssembly and host agree on 51 of 51 questions.
