---
date: 2026-10-20
title: Routes are found in the graph, run, and certified against a reference
files: [rsc/core/graph.hpp, rsc/spatium/gen_graph.py, rsc/spatium/graph.json, rsc/spatium/rsc_route.cpp, rsc/spatium/emit.py, rsc/spatium/rsc_find.py, rsc/spatium/test_rsc_find.py, rsc/spatium/find-report.md, tests/symmetry/rsc_pair_support.hpp, tests/test_rsc_graph.cpp, scripts/ci_changes.py]
open:
  - the independence rule is structural (no shared node of the graph); it cannot see what a node calls inside the library, so a product's member distance, which calls its factors, counts as independent of the product built from the same factors
  - five disagreements are listed in rsc/spatium/find-report.md and not investigated or fixed here (TorusChart exp_map over Dual, and the derived distance of TorusChart and SphereLevelSet against their member distance)
  - frechet_mean and verify_exp_log (algorithms over a set of points) have no emitter; project and normal have routes only as lifts of one member, so nothing independent
  - Real50 lifts are left out for the derived spaces and the erased surfaces (minutes per call); the typed surfaces have Dual lifts only
  - cheaper-than-reference is counted in calls and none was found; delivery of a cheaper route (the seam into the library) is not written
  - no model reads the graph and no store holds weights; RSC in CI on a schedule (rsc.yml) does not exist
---

The graph now has edges between types: a space over `Real50` or `Dual<double>` is a node of its own with an `entry` edge that computes there and converts back, a `ProductSpace` has `components` edges that need its factors' quantities, a `LieGroupManifold` has `group` edges that need the group's compose, inverse, exp and log. `rsc/core/graph.hpp` reads needs on other types (`Type#quantity`) and an `entry` flag; the matrix has no cell for these edges, so they sit at level 0.

`rsc/spatium/rsc_find.py` takes the routes `rsc_route --all --json` lists, turns each into C++ with `emit.py` (one function per edge, calling the member, ADL or derived form the edge names; nothing written per route), compiles it, and runs it against a reference on 300 inputs of the connectivity matrix's sampler, with a bound that is the sum of what the two paths promise. The first run, over every (type, operation) with two or more routes: 266 pairs, 1187 routes, 61 witness pairs written to `tests/symmetry/pairs/rsc_*.cpp` (53 files), 5 disagreements, 0 red cells, 0 routes cheaper than their reference, and an unmatched measure below 3/300 for the 1182 pairs that never failed. `rsc_find.py --all --check` is a CI job; `find-report.md` carries the per-route outcomes.
