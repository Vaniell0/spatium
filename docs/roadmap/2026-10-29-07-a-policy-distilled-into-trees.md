---
date: 2026-10-29
title: a tree that scores a row, grown by CART, and a policy distilled into trees beside its own description length
files: [include/spatium/algebra/models/tree.hpp, include/spatium/statistics/tree_fit.hpp, tests/test_tree.cpp, examples/chess_agent.cpp, modules/algebra/models_tree.cppm, modules/statistics.cppm]
open:
  - a linear policy is cheaper than any tree that comes near it, so distilling one into a tree is the wrong direction by description length; the tree pays where the model is not linear in its columns (the tanh net of the search, 481 weights), which is not measured here
  - CART fits the scores by squared error, not the distributions by divergence: on held-out chess positions the divergence of depth 2 is above that of depth 1
  - the tree is not lowered to the field IR yet, and neither is the linear policy, which needs only products, sums and a comparison
  - a bit per number is the 32 of a float; a tree's thresholds and scores could be stored in fewer, which is not tried
---

`algebra/models/tree.hpp` is a binary tree over rows of K numbers: an inner node sends a row left when one column is at most a threshold, a leaf holds a score. It counts its description length in bits (a bit of shape a node, the column in log2 K bits, a number in 32). `softmax` of scores at a temperature is the same distribution `LinearSoftmax` gives from its own, held to 1e-15. `statistics/tree_fit.hpp` grows one by the greedy rule of CART on (row, target) pairs; it recovers a two-level step function exactly (squared error below 1e-24, three leaves, every held-out row away from the thresholds scored right) and one split short it cannot, and its squared error on the sample never rises from depth 0 to 8. A tree of one leaf scores every candidate alike, and its divergence from a policy is log n - H(p), held to 1e-12.

A linear softmax with four weights, distilled into trees of depth 2 to 8 on 300 lists of eight random rows, gives on 200 other lists a divergence of 0.29, 0.13, 0.085 and 0.050 nats and the same choice on 46.5%, 65%, 75.5% and 79.5% of them, at 237, 1053, 4317 and 16557 bits; the policy itself is 128 bits with nothing lost. On chess (`chess_agent --observe 256 --distill 40`: the agent's 13 weights, 416 bits, agreeing with the hand-written player on 72.3% of 7612 positions of random games), trees grown on 204341 moves of 7545 other positions reach at depth 4, 6 and 8 a divergence of 1.19, 0.68 and 0.12 nats, the agent's choice on 52.9%, 59.5% and 77.5% and the player's on 50.5%, 54.6% and 71.2%, at 1083, 4303 and 13403 bits. By the rule that the simplest rung passing the measure is the one kept, the linear rung wins here: the tree of depth 8 nearly matches its agreement with the player at 32 times its bits.
