---
date: 2026-10-21
title: distance is the geodesic distance everywhere, or it is absent
files: [include/spatium/core/access.hpp, include/spatium/spaces/parametric.hpp, include/spatium/spaces/implicit.hpp, include/spatium/mesh/geodesic.hpp, include/spatium/mesh/voronoi.hpp, rsc/spatium/gen_graph.py, rsc/spatium/graph.json, rsc/spatium/find-report.md, tests/test_minkowski.cpp, tests/test_parametric_typed.cpp, tests/test_implicit_typed.cpp, docs/connectivity.md, tests/connectivity/expected.hpp]
open:
  - the erased ParametricSurface and ImplicitSurface no longer have a distance, so their MetricAxioms cells fall to L1 (nothing left to check the axioms on); the chord is their `chord_distance`
  - mesh/geodesic.hpp and mesh/voronoi.hpp now call spaces::distance, which for a space with a log_map and metric is |log|_g, not the straight edge length they used on an erased surface
  - the TorusChart exp disagreement of rsc/spatium/find-report.md remains (the retraction through find_params differs between double and Dual)
---

`core/access.hpp` says `distance(p, q) = |log_p q|_g`. Two places broke that. The erased parametric and implicit surfaces had a member `distance` that returned the straight chord through the ambient space; it disagreed with the derived one by up to 3.7e3 of its tolerance. The member is removed and the chord stays under `chord_distance`. `Minkowski<4>` got the derived `|log|_g`, which is zero on light-like and imaginary on time-like separations; the derived form is now withheld from a space that has `interval` and `proper_time`, so Minkowski offers the interval and no distance. `gen_graph.py` follows the same rule, and the disagreements of the find report go from five to one.
