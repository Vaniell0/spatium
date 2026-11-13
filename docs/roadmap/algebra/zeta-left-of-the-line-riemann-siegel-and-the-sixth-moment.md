---
date: 2026-11-13
title: Zeta left of the line, Z by Riemann-Siegel, and the sixth moment to 2e8
files:
  - include/spatium/algebra/zeta.hpp
  - include/spatium/answer/zeta_request.hpp
  - examples/zeta_moments.cpp
  - docs/zeta/sixth-moment.md
  - docs/zeta/sixth-moment-blocks.tsv
  - scripts/external/zeta/analyse.py
  - scripts/external/zeta/odlyzko_check.cpp
  - tests/test_zeta.cpp
open:
  - Z by Riemann-Siegel is in doubles with long double phases, its cut bounded and its rounding measured; an interval version (phases and C_k enclosed) would make the moments proved, not computed with a budget
  - Gabcke's bounds are taken as Gourdon quotes them, to C_2 only; his constants for C_3 and C_4 would tighten the bound below t = 10^4
  - the file of blocks was computed with phases in doubles; recomputed with long double phases it moves by at most 1.2e-8 a block (24 blocks checked), but it is not recomputed whole
  - the zeros of zeta are still proved by the Euler-Maclaurin sum, |t| / pi terms; Riemann-Siegel in intervals would take them past 10^5
  - the remainder of the sixth moment grows with slope 0.90 over three decades; whether it is sqrt(T) log^6 T or a power is not told apart
  - zeta.hpp is header-only in the C++23 modules build: as a partition, instantiations its inline functions make (quadrature, series, std::min/max<double>) were taken by importers as emitted there and were not, so spatium.render failed to link under GCC 15
---
`zeta(s)` left of Re s = 0 by the functional equation zeta(s) = chi(s) zeta(1 - s), the sine's argument reduced by whole multiples of pi so the trivial zeros are tight (zeta(-20) within 2e-37, where Euler-Maclaurin gave 392.98 at -13 and a NaN past -58), held against zeta(1 - 2n) = -B_2n / 2n for n to 30 and against Euler-Maclaurin where both run; `log_gamma` is now one generic Stirling series for intervals, doubles and dual numbers. The zeta page takes windows below t = 0 by the mirror. `riemann_siegel_z(t)`: the formula to C_4, Psi's derivatives from `Series` around p = 1/4 in the sinc form (dividing Psi's series at 1/2 loses 2^n of the n-th coefficient), C_0..C_4 from a table of 4096 steps (7e-13 from the derivatives, 0.2 us against 8 us), phases reduced in long double; against the proved Z to 1e-11 on [1e3, 3e5] and against mpmath at 30 digits to 9e-10 at 2.6e8 (6e-6 with phases in doubles). `riemann_siegel_truncation(t)` bounds the cut: Gabcke's |R_2| <= 0.011 t^{-7/4} for t >= 200, plus the C_3 and C_4 terms. `zeta_moment` integrates |zeta|^{2k} on the line adaptively, each unit's tolerance relative to the units before (an absolute one at the rule's own never settled); `zeta_moment_fixed` by Gauss-Legendre on each unit with every node's phases turned by n^{-i} from one unit to the next, sixteen units a pass, 600-1300 times faster, its `error` the bound on the cut. `examples/zeta_moments` sets the sixth moment against the CFKRS conjecture block by block in threads: to 2e8, 4000 blocks, the ratio 0.999979 at 2e8; CFKRS's Table 5.1.2 reproduced, its conjecture column to 1e-12 and its computed column low by 2e-7 to 4e-6 in every block, which three routes here agree on to 5e-13. Odlyzko's tables of zeros: the first 100000 and 9998 near zero number 10^12, all found, none extra, within the tables' stated accuracy. Write-up and data: `docs/zeta/sixth-moment.md`. The request moved into the answer door, `include/spatium/answer/zeta_request.hpp`.
