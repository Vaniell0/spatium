---
date: 2026-11-02
title: the tree search reads on past a value that is not a number
files: [rsc/core/tree_search.hpp, tests/test_tree_search.cpp]
open:
  - a state the net cannot read is counted, not named; a search over pieces of an interval will want to say which piece and why (a pole, an undefined region)
---

One NaN from the evaluator used to poison the whole tree. It was added into every edge on the path back to the root, and a NaN is never above or below anything, so from then on the search compared nothing. An infinite value took the search over instead: inf as a mean beats every finite one. A NaN prior meant the move it belonged to was never tried. In Nim from 9 stones, with the net unable to read one state two plies down, every standing at the root came out NaN; with +inf there, the losing take of two got 397 of 400 simulations.

Now a value outside [-1, 1], an infinity included, is the net's certainty and is clamped to it. A NaN value is a state the net could not read: its simulation backs nothing up, and since the state is already expanded, the next simulation through it goes on below. So the search reads what the net could not by playing it. A move's prior that is a NaN, infinite or below zero counts as none, and the rest are renormalised (uniform when nothing is left). `SearchResult` counts the states with an unread value and those with an unread prior, and a net that reads everything leaves both at zero.

The same position now gives the winning take of one 347 of 400 simulations with a NaN there, and 358 with +inf. Tests hold the search to the exact answer with a NaN value, with either infinity, and with the winning move's prior a NaN, below zero or infinite at every state.
