---
date: 2026-10-20
title: The library as a typed graph
files: [rsc/core/graph.hpp, rsc/spatium/gen_graph.py, rsc/spatium/graph.json, rsc/spatium/rsc_route.cpp, tests/test_rsc_graph.cpp, scripts/ci_changes.py]
open:
  - `find` over the graph (a chain of operations to a goal, with a sample or a residual as the target) and `witness` are not written; `routes` only lists the ways to one quantity on one type
  - the model that reads the graph, and the store for its weights, do not exist
  - RSC in CI (the scheduled `rsc.yml` that opens draft PRs on `rsc/*`, and `rsc_run --check` in `ci.yml`) does not exist; the CI here is the isolation grep and `gen_graph.py --check`
  - operations are the ones of `core/access.hpp` and the generic algorithms the connectivity probes exercise; the geometry primitives have satisfies edges and no operations, and `project` and `normal` are declared at level 0 because no matrix cell exercises them
---

`rsc/spatium/graph.json` is the library as a typed graph over its concept hierarchy, the first building block of RSC: the old prototype worked on `vector<double>` and so could not see spaces at all. Nodes are the 75 concepts of `include/`, the types (every space of the connectivity matrix, `Minkowski<4>`, the geometry primitives, `SO3`, `SE3`, `LieGroupManifold` over them, `FiniteSet`) and four scalars. Refinements are parsed from the concept bodies (only the terms of a top-level conjunction: `A<S> || B<S>` refines neither); which type satisfies which concept is not parsed but compiled, one translation unit per type printing `Concept<T>` for every concept that takes one type, so no compile is expected to fail. An operation is one way to get a quantity out of a type: `distance` as a member, by ADL, or derived as `|log_p q|_g` are three edges, each with the level the connectivity matrix reached for the probe that exercises it.

`rsc/core/graph.hpp` reads the JSON and answers `routes(type, quantity, depth)`: every distinct set of operation edges that closes over its needs, found breadth first, shortest first. It includes nothing of the library (a CI grep holds that), so the same code will read the graph of another one. `rsc_route <type> <operation>` prints them; `tests/test_rsc_graph.cpp` holds the graph against the compiler on a sample of (type, concept) pairs, asks `Sphere<2>` for `distance` (the closed form and the derivation, both at L3) and a space given by its metric alone, which has only the derived one. `gen_graph.py --check` in CI regenerates the file and fails if it would change.
