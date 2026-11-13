---
date: 2026-11-13
title: Muon, a step rule that moves a weight matrix as a matrix, and what it changes for chess_zero
files:
  - include/spatium/algebra/analysis/optimize.hpp
  - tests/test_network.cpp
  - examples/chess_zero.cpp
resolves: [rsc/a-tree-search-that-learns-chess-by-playing-itself#5]
open:
  - one matrix per rule; a net of several hidden matrices would want one Muon buffer for each
  - the orthogonalisation is five fixed Newton-Schulz steps in the working scalar; whether fewer do on matrices this small is not measured
  - on chess_zero's 13-32-1 net it does no better than Adam; a net where the matrix carries most of the weights is where it would have to show
---

`Muon<N, R, C, T>` in `algebra/analysis/optimize.hpp` is the rule of the reference muon.py of github.com/KellerJordan/Muon, read whole: the first R x C parameters are a matrix, the rest go to an Adam inside it with the reference's defaults (betas 0.9 and 0.95, epsilon 1e-10) at a rate scaled from the matrix's. The matrix keeps a momentum buffer (0.95, Nesterov), and its step is that update orthogonalised by five steps of the quintic Newton-Schulz iteration with coefficients (3.4445, -4.7750, 2.0315), normalised first by the Frobenius norm plus 1e-7 and worked in the wide orientation, then scaled by max(1, R / C)^0.5. The rate is in units of the spectral norm, 0.02 by default. The iteration does not land the singular values on 1 but between about 0.5 and 1.5, which the reference states and the test holds.

The tests: a 6 x 3 gradient with singular values 10, 1 and 0.1 comes back with all three in (0.5, 1.5) and nothing off its diagonal above 1e-12, and the same gradient scaled by a thousand gives the same step within 1e-7, the epsilon beside the norm. Fitted by Muon at 0.02, the tanh net learns tanh(2 x0 - x1) to a mean squared error below 0.01, as it does with Adam.

`examples/chess_zero --rule muon` trains the value net's 32 x 13 hidden matrix by Muon and its 33 other weights by the same Adam as `--rule adam`, so the two runs differ in that matrix alone. Over three seeds of 100 iterations (16 games, 32 simulations a move), measured every 25 iterations in 40 games against each fixed player, the pooled points are 0.683 for Adam and 0.672 for Muon against the greedy player and 0.285 and 0.278 against the two-ply material search; at iteration 100, 0.762 and 0.671, 0.283 and 0.229. One evaluation of 40 games moves by about 0.08 from seed to seed and iteration to iteration, so the rules are not told apart here. Both reach about 0.99 against a random player by iteration 25 and stay level after it, against the greedy player between 0.44 and 0.86: on this net the step rule is not what holds the player back.
