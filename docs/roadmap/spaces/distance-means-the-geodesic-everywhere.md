---
date: 2026-10-21
title: distance is the geodesic distance everywhere, or it is absent
files: [include/spatium/core/access.hpp, include/spatium/spaces/parametric.hpp, include/spatium/spaces/implicit.hpp, include/spatium/mesh/geodesic.hpp, include/spatium/mesh/voronoi.hpp, rsc/spatium/gen_graph.py, rsc/spatium/graph.json, rsc/spatium/find-report.md, tests/test_minkowski.cpp, tests/test_parametric_typed.cpp, tests/test_implicit_typed.cpp, docs/connectivity.md, tests/connectivity/expected.hpp]
open:
  - the erased ParametricSurface and ImplicitSurface have no distance at all: their log is a retraction, and |log|_g of it is neither symmetric nor subadditive (it failed the metric axioms over Real50 in a CI build); the chord is `chord_distance`, and the matrix probes and the mesh edge lengths use it for them
  - frechet_mean has nothing to minimise on an erased surface, so its cells there are L1
  - the TorusChart exp disagreement of rsc/spatium/find-report.md remains (the retraction through find_params differs between double and Dual)
---

`core/access.hpp` says `distance(p, q) = |log_p q|_g`. Two places broke that. The erased parametric and implicit surfaces had a member `distance` that returned the straight chord through the ambient space; it disagreed with the derived one by up to 3.7e3 of its tolerance. The member is removed and the chord stays under `chord_distance`. `Minkowski<4>` got the derived `|log|_g`, which is zero on light-like and imaginary on time-like separations; the derived form is now withheld from a space that has `interval` and `proper_time`, so Minkowski offers the interval and no distance. The same withholding applies to the erased surfaces (`exact == false`). `gen_graph.py` follows the same rule, and the disagreements of the find report go from five to one.
