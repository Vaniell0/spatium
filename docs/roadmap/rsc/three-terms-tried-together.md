---
date: 2026-10-30
title: synthesis tries sets of three short terms together, held to a measure of recall
files: [rsc/spatium/synthesis.hpp, tests/test_rsc_synthesis.cpp]
open:
  - sets of four are still grown greedily from the best of three; on [0.1, 100] 52 of 60 sums of four terms come back
  - the Thomas-Fermi solution near the origin is found as pieces, not as its series 1 + B x + 4/3 x^(3/2) + 2B/5 x^(5/2) + ...: the fourth term is six nodes
---

Run on the Thomas-Fermi solution near the origin, 1 + B x + 4/3 x^(3/2), three terms of six nodes in all, synthesis accepted nothing on [0, 1e-3] or [0, 0.1]: sets of one and two terms were tried in full and a third was added greedily, and on a short interval the terms are nearly collinear, so the greedy pick took whichever leaned toward the residual. Sets of three among the shortest terms (about 60 of them, 34220 sets, 37 ms) are now tried in full when no pair is accepted on the check grid; a pair that is accepted is kept, since the smaller set wins.

The measure is a symmetry: synthesis is the inverse of evaluation, so a sum of terms of its own language with coefficients no fraction writes, evaluated on an interval, must come back accepted in no more terms than it has. On 60 random sums per count of terms, before and after:

| interval | 3 terms | 4 terms |
|---|---|---|
| [1, 1.01] | 60 -> 60 (49 -> 60 in at most three) | 59 -> 60 |
| [0.5, 3] | 46 -> 60 | 36 -> 56 |
| [0.1, 100] | 45 -> 54 | 33 -> 52 |

The test holds sums of up to three terms on a narrow and an ordinary interval. With the search wider the deflection of light without the logarithm takes five pieces of 106 nodes at 1e-4, against eight of 120 before.
