---
date: 2026-10-30
title: synthesis searches a size whole, and drops a term that only refits the noise
files: [rsc/spatium/synthesis.hpp, tests/test_rsc_synthesis.cpp]
open:
  - a sequence known only at the integers (a formula from a list of (n, a_n)) is not a question synthesis can ask: it fits and certifies on a continuous grid, so a table read through an interpolating lambda gives 32 pieces of the interpolant, none accepted
  - a step function (the number of ways g(n) to write n as a sum of two primes) is 256 pieces, none accepted
  - 2n - 1 is found exactly but written n - (1 - n), one node fewer than 2 * n - 1
---

Run on the growth 2n / (log n)^2 of Goldbach pairs on [10, 1e6], with the logarithm among the operations and terms of up to seven nodes, synthesis returned 22 pieces instead of the one term n / (log(n) * log(n)). Two causes.

The budget of candidates (4000) stopped the enumeration inside size seven: it held 2483 of that size's 4837 expressions, whichever came first in the loop, and the one wanted was not among them; the whole size costs 43 ms. The budget is now spent by whole sizes -- no larger size is begun once it is reached, a size begun is finished -- so the terms are every expression up to the last size they hold.

With the term found, the answer was 2 * n / (log(n) * log(n)) + 0 * n / log(log(n)): the pair the search starts from is drawn from the shortest terms, a swap then replaced one of them by the long term and the error fell to rounding, and the other stayed with a coefficient of 1e-17, which the fraction step wrote as 0. A chosen term without which the error does not grow is now dropped, the description-length rule applied after the swaps. The answer is 2 * n / (log(n) * log(n)), one term, error 0 on the check grid. With six nodes it is not found, the case the test must fail.
