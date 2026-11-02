---
date: 2026-10-29
title: entropy, cross-entropy and divergence of distributions, and distilling a softmax at a temperature
files: [include/spatium/statistics/information.hpp, include/spatium/algebra/models/softmax.hpp, tests/test_information.cpp, modules/statistics.cppm]
open:
  - the first callers are the tests and a chess project outside this repository; the tree distilled from a policy, whose fidelity these measure, is the next change
  - the measures are over vectors of probabilities; a density, or a distribution given by samples, has no estimator here
  - nothing turns a cross-entropy in bits into an admission rule for a learned component yet, though it is the description length of what the component explains
---

`statistics/information.hpp` gives the entropy of a probability vector, the cross-entropy of one against another and the Kullback-Leibler divergence, in nats (`in_bits` for bits), templated on the scalar. A term with probability zero contributes nothing, and a q that gives zero to an outcome p can draw makes the cross-entropy and the divergence infinite, which the tests hold. The divergence equals the cross-entropy less the entropy on 50 random pairs, is never negative and is zero only at p = q, and on 41 or more of the 50 pairs it differs from its reverse by more than 1e-6. It adds over independent draws: for binomials of n = 1, 5, 20 and 60 draws, built from `discrete::binomial_coefficient`, it is n times the divergence of one draw to 1e-11.

`LinearSoftmax`'s cross-entropy and gradient now take a temperature (default 1, which leaves every earlier caller as it was): the target is a distribution at that temperature and the gradient is the difference of the two distributions over the temperature. Distilling a teacher is that with the teacher's distribution at the same temperature, scaled by T^2: T^2 times the gradient is held against Dual carried through `statistics::kl_divergence` of the two distributions at T = 1, 2.5 and 8, to 1e-11, and the cross-entropy at T is the divergence plus the teacher's entropy. A student distilled at T = 2 from a teacher's soft distributions on 1500 lists, by AdaGrad, chooses as the teacher does on more than 95% of 400 lists it never saw; at rate zero it stays below 30%, chance being 10%.
