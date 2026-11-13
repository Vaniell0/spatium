---
date: 2026-11-06
title: The zeros of zeta proved, and counted
files:
  - include/spatium/algebra/zeta.hpp
  - include/spatium/answer/zeta_request.hpp
  - web/src/zeta.ts
  - site/gallery/zeta/index.html
  - tests/test_zeta.cpp
open:
  - completeness above a few thousand rests on Backlund's bound at the window's ends, which needs Re zeta > 0 on a segment; where no height near an end has it, the count is not decided -- Turing's method would decide it from Z alone
  - the Backlund check cuts sigma into pieces evaluated naively in intervals; a centred form with zeta' would need far fewer
  - a zero whose bisection meets an undecided middle stops there (Lehmer's first zero holds 8 digits, its neighbour 11)
  - the zeros of zeta', the |zeta| surface, the moments, the two-dimensional cut of |zeta| > eps
  - the Euler-Maclaurin sum costs |t| / pi terms; Riemann-Siegel would cost sqrt(t / 2 pi)
---
`algebra/zeta.hpp`: `zeta(s)` for a complex interval s holds the true value -- the Euler-Maclaurin sum with N >= (|s| + 59) / pi, 29 Bernoulli terms (B_2..B_60 as fractions, read once), and the remainder bounded by |s + 59| / (Re s + 59) times the first term left out; the same template in doubles is the picture's sum, and the double lies in the proved interval (tested at four points from sigma = -0.5 to t = 700). `riemann_siegel_theta` is Im log Gamma(1/4 + it/2) - (t/2) log pi by Stirling's series after shifting |z| past 16, its remainder bounded by 2^{m+1} times the first term left out (`atan` of an interval is new in `core/interval.hpp` for the argument). `hardy_z` is Re(e^{i theta} zeta(1/2 + it)); its imaginary part holds 0 at t = 10, 100, 1000 and a theta moved by 1e-9 does not, which checks theta against zeta. `gram_points(t0, t1)` proves each g_n between two numbers where theta is below and above n pi, and whether Gram's law holds there. `critical_zeros(t0, t1)` brackets every sign change of Z, narrows each by halving while the middle's sign is proved, and counts: where Re zeta > 0 on [1/2, 2] + iT (proved by `prove_on`), |S(T)| < 1/2 and N(T) is the whole number nearest theta(T)/pi + 1 (Backlund); with N decided at heights below t0 and above t1, chosen where Re zeta(1/2 + iT) is largest, equal counts prove the zeros are all there are between them, each simple and on the line. Measured: the 29 zeros below 100 against mpmath's to 25 digits, every one inside its interval, widths below 1e-10, complete, 0.18 s; Gram's law first fails at g_126 (282.4547), the window [280, 284] complete with 9 of 9; Lehmer's pair at 7005.0628662 and 7005.1005647 complete with 4 of 4 in 0.59 s (5.8 s before the heights were ordered, 90x slower before the Bernoulli numbers were cached). The page `/gallery/zeta/` asks the core in its worker (`spatium_zeta`, `include/spatium/answer/zeta_request.hpp`: hardy, zeros, gram, grid, notes) and draws Z(t) with the zeros and the Gram points, zeta over sigma in [-1, 2] as colour (hue the argument, a band where |zeta| is a power of 2, the zeros as rings on the critical line), and the notes the library writes: the count that proves completeness, each Gram point where the law fails, zeros closer than a quarter of the mean spacing. Host and WASM agree on five zeta requests in `site/wasm/questions.txt`.
