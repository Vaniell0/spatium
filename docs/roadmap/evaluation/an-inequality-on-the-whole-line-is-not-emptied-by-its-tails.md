---
date: 2026-11-16
title: an inequality on the whole line is not emptied by its tails -- a piece left undecided only by being wide is proved again, and one with no root and a bounded value is decided by a point
files: [include/spatium/answer/evaluate.hpp, tests/test_answer_respond.cpp]
open:
  - a pole at an irrational point (1/(x - sqrt 2)) stays an undecided piece: no fraction names it and its value is unbounded, so neither a break nor a sign is proved there
  - log_{x-3}(x) >= x leaves the pole at 4 undecided, where the value grows without bound on both sides; a break found by the sign of 1/f would close it
---

`(x - 1)/(x + 2) > 3` came back as the empty set, marked proved; on `[-10, 10]` it was `(-7/2, -2)`. The interval method cuts the whole line breadth first, and its 4096 evaluations went to the tails at magnitudes near 1e306 before the middle was cut fine enough: on `[-46, -16]` the quotient reads `[-2.6, 0.36]` for a value near -2, the piece stayed undecided, and a run of undecided pieces with no root was written down as a break, a place with no value. The root at -7/2 sat in one run with the pole at -2 and went with it.

Two changes. Each run left undecided is proved again on its own with a fresh budget, up to three rounds, while a round still decides something. And a run on which no root is proved is a break only where its value is unbounded or empty: a bounded value on the whole run means a value at every point and no pole, so the sign cannot change and one point decides it. The respond test holds the line against the bracket for this inequality; the full suite (1503 cases, the matrix aside), the gallery and the reference are unchanged.
