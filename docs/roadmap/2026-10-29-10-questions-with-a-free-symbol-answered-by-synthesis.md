---
date: 2026-10-29
title: questions with a free symbol -- a root, an integral -- answered by synthesis, not by rewriting
files: [rsc/spatium/query.hpp, rsc/spatium/synthesis.hpp, tests/test_rsc_query.cpp]
open:
  - two kinds of question, a root on a bracket and an integral over an interval; a limit, a series and an equation of motion have their references in the library (Series, integrate_extrapolated) and no question yet
  - one free symbol; a question in several parameters needs the synthesis in several variables
  - the scene description does not take questions yet: an equation in a scene, solved when it is built, is the use this is for beside a person's question
  - the root's reference is bisection on the bracket the question gives; finding a bracket is the asker's
---

`rsc/spatium/query.hpp` answers a question with a free symbol. The question is a field of the scene IR in two symbols, U the unknown or the variable of integration and V the parameter, and its kind says what the library computes for it: a root on a bracket by bisection for `solve`, the integration door for `integrate`. That number at each value of the parameter is the reference, and the answer is a field of the parameter synthesised on its interval by `synthesis.hpp`, in pieces where one expression does not reach the tolerance. Nothing is rewritten: where the operations hold a closed form it is found, where they do not the pieces are the answer and say so. A question that has no value somewhere on the interval is refused before anything is fitted.

The root of x^2 - p on [0, 10] for p from 0.5 to 4 is one piece of one term, sqrt(p), to 1e-10 at 200 parameters it was not fitted on. The integral of x^p over [0, 1], written exp(p log x) with the IR's new logarithm and exponential, is one piece, 1 / (p + 1), to 1e-9. The integral of exp(-p x^2) over [0, 1] is sqrt(pi / p) erf(sqrt p) / 2; erf is not among the operations, and the answer is twelve pieces, held to 1e-6 against erf. A root of x^2 + p on [0, 10] is refused: there is no change of sign.

The synthesis now counts every error below rounding as the same error, so among sets that fit to the noise the smaller one wins -- the rule of description length; before it, sqrt(p) came back with three more terms that only refitted the noise.
