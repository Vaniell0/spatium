---
date: 2026-10-08
title: The 3 by 3 eigen solver is Jacobi now
stage: 1
files: [include/spatium/spaces/spd.hpp, tests/symmetry/pairs/spd3_near_identity.cpp]
open:
  - Jacobi iterates (a few sweeps, at most 32) where the 2 by 2 is a formula; through `Dual` the iteration is differentiated by the same code as any other, which is fine apart from ties, and the function of a matrix does not go that way (it uses the values' decomposition and the Daleckii-Krein formula), but a caller who takes `eigen_sym` itself over `Dual` at a tie still gets a derivative of the labelled eigenvalues, which is not defined there
  - general N is still not built; SPD(n) stops at 3
---

The entry before this one listed it as open: the 3 by 3 solver took its eigenvalues from the characteristic cubic and its vectors from cross products of rows of S - lambda I, and on I + hM the reconstruction V diag(l) V^T missed S by 1.7e-3 at h = 1e-2, by 1.7e-5 at 1e-4, and at h = 1e-5 by 0.95 with the vectors no longer orthonormal. Anything built on a function of an SPD(3) matrix near an isotropic one, a covariance of a nearly round cloud or the inertia tensor of a nearly round body, was wrong by an amount that depended on how round it was.

Cyclic Jacobi rotations replace it. A rotation zeroes one off-diagonal entry of a symmetric matrix exactly and keeps V orthogonal to rounding whatever the eigenvalues are, so a repeated or a nearly repeated eigenvalue costs nothing and the three special cases of the cross-product method (distinct, two equal, three equal) are one code path. Measured on the same matrices, the reconstruction error is 1e-15 and the orthonormality error 1e-16 at every h from 1e-2 to 1e-7. The pair `spd3_near_identity` holds it: sqrt, inverse sqrt and log of c I + hM for h log-uniform between 1e-9 and 1e-1, double against `Real50`, to 1e-13; it passes with a worst distance of 3% of that. The existing tests of `eigen_sym` (known eigenvalues, axisymmetric, scalar multiple of the identity) pass unchanged, since none of them depends on the order or the sign of the vectors; the values come back sorted descending, as the 2 by 2 does.
