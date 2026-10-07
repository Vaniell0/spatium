---
date: 2026-10-07
title: A Lie group as a space
stage: 2.1
open:
  - the exp and log of SO(3) and SE(3) read through core/access.hpp for riemannian_minimize
  - rotations as a row of the connectivity matrix
  - SE(3) has the Cartan-Schouten structure only: no geodesics of a left-invariant metric
resolves: [2026-10-07-05#2, 2026-10-07-06#3, 2026-10-07-07#1, 2026-10-07-08#1]
---
Stage 2.1, the last of its four. `spaces/lie_group.hpp`: `LieGroupManifold<G>` for any `LieGroup` -- points are group elements, tangent vectors are algebra elements carried by left translation, `exp_map(p, v, t) = p exp(t v)`, `log_map(p, q) = log(p^-1 q)`, the metric the algebra's inner product at every point -- so the algorithms written against exp, log and a metric (`frechet_mean`, `geodesic`, `midpoint`, `verify_exp_log`, `verify_metric`) run on rotations and rigid motions with no copy for groups. For SO(3) this is exactly the Riemannian geometry of the bi-invariant metric (geodesics are the one-parameter subgroups, the distance is the rotation angle), so the Frechet mean is the rotation average of SLAM and sensor fusion; SE(3) has no bi-invariant metric, so there exp and log are the group's own (the Cartan-Schouten connection) and invert each other, which is the mean everyone computes for poses but not the Levi-Civita geodesic of any metric -- said in the header, not left to be found. Held by paths that share nothing: the mean of two rotations about one axis against the half angle (1e-9) and against the customization points' midpoint (1e-12); equivariance, the mean of left-rotated data against the rotated mean (1e-8); the Karcher condition read off the logs, their sum at the mean zero (1e-8); `verify_exp_log` and `verify_metric` on rotations, and exp/log round trip on poses. Not done: the Lie-group exp/log of SO(3) and SE(3) read through `core/access.hpp` for `riemannian_minimize`, and rotations as a row of the connectivity matrix.
