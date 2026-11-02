---
date: 2026-10-22
title: gravity of many sources, free or on a torus, with a Barnes-Hut tree on the BVH
files: [include/spatium/physics/mechanics/gravity_field.hpp, include/spatium/physics/mechanics/gravity_tree.hpp, tests/test_gravity.cpp, examples/nbody_demo.cpp]
open:
  - the tree carries monopoles only; on a torus the damped short-range kernel makes it about ten times less accurate than in free space at the same opening angle
  - the reciprocal sum visits every wave vector for every point, with no mesh and no FFT
  - there is no tree, and no periodic sum, for a curved space
---

`PointGravity` is one source with no softening and a fixed law. `GravityField` takes many sources in N dimensions, an integer law exponent (N-1 is a point source's own field; 1 is the softer law), Plummer softening, an optional cap on the acceleration, a per-body ratio of gravitational charge to inertia (`GravityPull` feeds it to the existing single-body integrators), and a periodic box with a uniform neutralizing background. The box is summed by Ewald's method with the softening inside the split, so the answer is exact for every exponent from 1 to N-1 and does not depend on the split. `GravityTree` evaluates the short-range part on `spatial::BVH` with the monopoles in a side array indexed by node; the BVH is untouched. `examples/nbody_demo.cpp` is the caller: leapfrog on a torus, tree against direct, energy and momentum drift, an SVG of the result. Figures and the method are in `docs/ROADMAP.md` next to the BVH items.
