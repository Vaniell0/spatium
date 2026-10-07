---
date: 2026-10-07
title: Subdivision goes through the space's own midpoint
stage: 2.1
closes: [25]
open:
  - Schild's ladder at full length
  - riemannian_minimize without HasNormal
  - the geodesic and Lie-group algorithms on the customization points
---
Stage 2.1, findings 25. `subdivide_once` and the adaptive subdivision put a new vertex at `space.project` of the chord's midpoint, which is the point halfway along the geodesic only where the projection is radial -- a sphere. On `Hyperbolic` it projects vertically, off the geodesic and not halfway: on the triangle of the new test an edge of hyperbolic length 2.59 had its new vertex at 1.34 and 1.39 from its two ends, where half is 1.29 (the two sum to more than the edge, so the vertex is off the geodesic; on another edge 1.25 and 1.34 against 1.25), and the area of a refined mesh drifted from the true one (9.42 against 2.23 after five levels, as measured in the audit). Both now call `spaces::midpoint` of `core/access.hpp` -- a member if the space has one, an ADL function if not, else `geodesic(p, q, 1/2)` from its own exp and log -- so the rule is the space's, not a projection that happened to be right on the first space tried. Held by a test written to fail first (it did): after one subdivision of a hyperbolic triangle every new vertex is on the hyperboloid and at half the edge's hyperbolic length from each end, to 1e-10; the 68 existing mesh, adaptive, geodesic, transport, Voronoi, topology and heat-method tests are unchanged. Not done: Schild's ladder at full length (findings 16), `riemannian_minimize` without `HasNormal`, the geodesic and Lie-group algorithms on the customization points.
