---
date: 2026-11-13
title: a law from noisy data, chosen by its description length, its constants named only where chance cannot have put them
files: [include/spatium/answer/law.hpp, include/spatium/answer/synthesis.hpp, examples/zeta_laws.cpp, tests/test_answer_law.cpp]
open:
  - the fourth moment to 5e6 (zeta_moments 2) does not give Ingham's law: fitted directly by x log(x)^k, its data resolve the top three powers and no more (the residual is the same with four and five, about 200 sqrt(T)), and the leading coefficient comes out 0.0523 against 1/(2 pi^2) = 0.0507; the sixth moment to 2e8 resolves 5 of 9 (docs/zeta/sixth-moment.md). Terms x log(x)^k past k = 2 are also beyond enumeration to size 5
  - the named constants are a fixed basis (1, gamma, log 2pi, pi, 1/pi, log(2pi)/pi, gamma/pi, log 2): a prior, from which the search chooses and which it cannot extend
  - at noise 0.3 on [1e3, 5e6], 1 and 0.096 log(x) with the slope of x moved are told apart in some draws and not in others
  - the points of a cumulative count or integral are correlated, and the standard errors take them as independent
  - the constant 7/8 of the count of zeros is found to 1 standard error (0.837 +- 0.036) and not named: the jumps of the count leave it too uncertain
---

`synthesize` answers a question whose values are exact and stops at the first sum within its tolerance; on data that carries a remainder it found 0.1582 x log x for the count of zeros where the law has 1/(2 pi). `synthesize_law(x, y, sigma)` (`answer/law.hpp`) is the other job. Its terms are the same enumeration over x alone (`SynthesisSpec::scaled_variable` false); a sum of them is fitted by least squares weighted by the noise at each point, with standard errors from the factor R, and chosen by its description length: n ln(RSS / n), plus ln n a coefficient, plus a term's nodes times the log of the symbols it is written in. Pairs among the first 600 terms are tried in full, triples among the first 60, and the 20 best pairs each with any third.

A coefficient is then read as a named constant, one at a time and the surest first, the others fitted again once it is fixed: a combination of 1, gamma, log 2pi, pi, 1/pi, log(2pi)/pi, gamma/pi and log 2 with small whole or fractional coefficients, within three standard errors, taken only when the combinations no costlier than it, times the window, over their spread, give below 1/100 -- among the forty-odd thousand combinations one lies near any number.

From the library's own numbers, blind: the second moment of zeta to T = 5e6 (`zeta_moments 1`, the noise taken as T^(1/4)) gives T log T + (2 gamma - 1 - log 2pi) T, Hardy and Littlewood's law in Ingham's form, its coefficient -1.6834426 +- 5.7e-6 against -1.6834457, named with 0.0073 expected by chance. Odlyzko's first 100000 zeros give the count (1/2pi) T log T - ((1 + log 2pi)/2pi) T + 0.837 +- 0.036, the leading terms of Riemann and von Mangoldt's, its 7/8 left as digits. A coefficient that is no such constant (-1.7) is not named. `examples/zeta_laws` prints both.
