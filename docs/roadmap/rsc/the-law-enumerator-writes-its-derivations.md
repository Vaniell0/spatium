---
date: 2026-11-16
title: the law enumerator writes each law it derives, with the chain that shows it, as training items
files: [rsc/spatium/rsc_laws.cpp]
open:
  - 590 items at size 7 (depth 1-6) is too few to train on; a generator that applies the premises as rewrite rules to random terms, checked on the same tuples, is the training set, and these items its held-out test of real laws
  - a congruence hop names the argument equality it rests on but not that equality's own chain
---

The enumerator drops a law the closure of smaller ones already holds; it now remembers why every two terms were joined -- a kept law, or congruence with the arguments already equal -- and `rsc_laws <size> <file>` writes each dropped equality as a question: the premises its chain uses and a few others true on the same spaces, the chain an equality a hop, its depth, what a one-hop rule would answer, and near misses a tuple of those spaces shows false. At size 7: 590 items, depth 1 to 6 (67, 100, 121, 202, 92, 8); the 637 laws it reports are the snapshot's, unchanged. A recurrent model that computes in its own latent steps learns depth past three hops only with a target at every step; these chains are that target, with an exact executor behind them.
