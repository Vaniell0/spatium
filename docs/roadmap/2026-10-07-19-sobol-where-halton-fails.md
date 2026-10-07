---
date: 2026-10-07
title: Sobol where Halton fails
stage: 2.7
files: [include/spatium/algebra/monte_carlo.hpp, tests/test_monte_carlo.cpp]
resolves: [2026-10-07-12#4]
open:
  - the direction numbers go to dimension 40 (the Joe-Kuo file has 21201)
  - the copies are randomised by a digital shift only; Owen scrambling is not done
  - the net property needs a power of two points per copy, and nothing enforces it (the default 65536 / 16 shifts is 4096)
  - an integrand of high effective dimension gains nothing from any of the sequences, and nothing tells the caller which kind theirs is
---
The user asked whether a Sobol sequence would not be more efficient than Halton for higher dimensions, and it is. `quasi_monte_carlo` now defaults to Sobol (`Sequence::Sobol`; Halton stays as an option up to 32 dimensions): Gray-code order, the Joe-Kuo direction numbers of dimensions 2 to 40 (S. Joe and F. Y. Kuo, "Constructing Sobol sequences with better two-dimensional projections", SIAM J. Sci. Comput. 30 (2008), file new-joe-kuo-6.21201, fetched rather than recalled), and a random digital shift per copy -- an XOR of every point with one random word per coordinate, which keeps the net structure that Sobol's quality is made of, where Halton's uniform shift mod 1 does not -- the standard error being the spread of the copies' means. Held by what makes it a Sobol sequence: the first points of the first two coordinates are the known ones (0, 1/2, 3/4, 1/4, 3/8 and 0, 1/2, 1/4, 3/4, 3/8); those two coordinates are a (0, m, 2)-net, exactly one of the 64 points in every elementary box 2^-a by 2^-(6-a), for all seven splits; and every one of the 40 coordinates is stratified, the first 1024 points one in each interval of width 1/1024, which holds the whole direction-number table to account. Measured on a smooth integrand of low effective dimension (the product of exp(x_i/K)) at 65536 points, by the standard error each reports (and the estimate within four of the exact value): Sobol against Halton 8x smaller at K = 16 and 18x at K = 32 (Halton has no K = 40), and against random points 65x smaller at K = 40; the product of 2 x_i, whose effective dimension is its dimension, gains nothing at K >= 16 from any of the three and every estimator's own error estimate is consistent with its true error -- known behaviour of the method, here measured. A request for more dimensions than the sequence has is a `Failed`, not a quiet wrong answer.
