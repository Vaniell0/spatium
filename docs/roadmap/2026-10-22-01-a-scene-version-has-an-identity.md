---
date: 2026-10-22
title: a scene version has an identity, so regeneration is told what changed
files: [include/spatium/io/field.hpp, include/spatium/io/build.hpp, include/spatium/render/cooked_scene.hpp, include/spatium/spaces/sample.hpp, tests/test_content_identity.cpp, tests/test_cooked_update.cpp, tests/test_build_dsl.cpp, benchmarks/bench_instance_update.cpp, examples/scatter_regen.cpp]
open:
  - the generator's time is part of its identity, so a scene regenerated at a new t is replaced whole even where the motion does not read t
  - a scatter with any zero-scale object is compared instance by instance, not by generator
  - lay_out still lays out the static triangles of the cooked scene on every call (about 17 ms for a 240x120 torus) before the cache is consulted
---

Fields have a Merkle content id and an exact `same_value`; trace nodes have a version that separates a scatter's generator from its count; cooked objects carry their scatter and site. `lay_out(previous, ...)` compares generators instead of hashing instances (at 1e5: 67 ms to 25 ms), and `cook(..., &previous)` places only the new sites of an unchanged scatter through a resumable sampler (103 ms to 8 ms). Results equal a fresh cook and lay-out for each kind of edit. Figures are in `docs/ROADMAP.md` next to the Field and BVH items.
