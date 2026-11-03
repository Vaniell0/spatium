---
date: 2026-10-23
title: gravity and flight on the closed 3-sphere, with the two stereographic charts
files: [include/spatium/physics/mechanics/sphere_gravity.hpp, include/spatium/physics/mechanics/manifold_body.hpp, include/spatium/spaces/sphere_charts.hpp, tests/test_sphere_gravity.cpp, tests/test_sphere_charts.cpp, examples/sphere_nbody_demo.cpp]
open:
  - there is no tree for S^3: the field is the direct sum over sources, O(sources) per point
  - the long free flight drifts with the speed that `geodesic_step` re-reads each step (4.6e-11 after 2e4 steps of 0.01)
  - the flight step and the charts are for the sphere; hyperbolic space has the same chart maps but no flight under a force yet
---

The potential of a unit mass on S^3 with a neutralising uniform background is `G(psi) = -(pi - psi) cot(psi) / (4 pi^2)`; it matches the zonal eigenfunction expansion (4e-6, with the constant -1/(8 pi^2) for the mean) and solves the equation with the right flux, so it stands as written. `SphereGravity` sums softened sources with a coupling per body and an optional clamp; the force is smooth and zero at the antipode. `flight_step` is kick - exact great circle - kick: reversible (1e-14), second order (energy error 6.4e-8, 1.6e-8, 4.0e-9 for dt 4e-3, 2e-3, 1e-3), and with no force the geodesic itself. The six components of angular momentum of SO(4) are conserved to rounding, and a clamp breaks them. A fan of equal-speed shots from P arrives at -P to 3.6e-14 with no masses and spreads in proportion to the mass with them. `SphereCharts` gives both stereographic charts of S^N, the inversion between them, and velocities through the conformal factor; a flight crossing the equator keeps the same coordinates at the rim. Figures are in `docs/ROADMAP.md` next to the gravity items.
