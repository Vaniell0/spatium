---
date: 2026-10-07
title: Two paths through the torus, and the conditioning of a cubic
stage: 1
files: [tests/symmetry/pairs.hpp, tests/test_symmetry.cpp, include/spatium/spaces/spd.hpp]
open:
  - eigen_sym for a 3 by 3 matrix takes its eigenvalues from the characteristic cubic, so near a double eigenvalue it is good to about sqrt(epsilon) times the scale where a Jacobi or QR sweep is good to epsilon times the norm; the SPD log and exp inherit that, and the pair now states it instead of hiding it
  - the eigenvectors of eigen_sym are not in any pair, and a repeated eigenvalue's null space is the place they are most likely to disagree
  - the bound of the SPD pair is a model of the cubic's conditioning fitted to what was measured (about four times the worst), not a derived inequality
---

Two more pairs in the registry. `ray_torus`, the closed-form quartic, against the same torus as a `ParametricSurface` hit by `ray_parametric`, a Newton search over the chart that shares no code with it: over 3 000 rays no hit is gained or lost and the worst distance is 3.5e-7 against the route's own tolerance of 1e-6. And the eigenvalues of a symmetric 3 by 3 matrix, double against `Real50`, on matrices Q diag(1, 1 + gap, 3) Q^T with the gap log-uniform between 1e-10 and 1, so that the near-degenerate case is where the draws are, not an accident of them.

The second pair found that `eigen_sym<Real50>` did not compile: a unary minus on a Boost number is a lazy expression, and `solve_cubic<T>` cannot deduce `T` from four arguments of three types. Two casts. Then it measured something the code comment had only asserted. Away from a double eigenvalue the error is a few ulps; at a gap of 1e-3 it is 1.3e-12, at gaps around 1e-7 it reaches 2.2e-7, and below 1e-8, where the two eigenvalues are one for every purpose, it saturates near 5e-9. The bound of the pair is the conditioning of a cubic root, which makes the tolerance a function of the input as the plan asked of the registry: the coefficients carry a rounding of eps S^3, the root moves by that over its slope, which is about the gap, and the slope cannot fall below sqrt(eps) S.
