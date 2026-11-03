---
date: 2026-10-22
title: a scene version has an identity, so regeneration is told what changed
files: [include/spatium/io/field.hpp, include/spatium/io/build.hpp, include/spatium/render/cooked_scene.hpp, include/spatium/spaces/sample.hpp, tests/test_content_identity.cpp, tests/test_cooked_update.cpp, tests/test_build_dsl.cpp, benchmarks/bench_instance_update.cpp, examples/scatter_regen.cpp]
open:
  - a scatter with any zero-scale object is compared instance by instance, not by generator
---

Fields have a Merkle content id and an exact `same_value`; trace nodes have a version that separates a scatter's generator from its count; cooked objects carry their scatter and site. `lay_out(previous, ...)` compares generators instead of hashing instances (at 1e5: 67 ms to 25 ms), and `cook(..., &previous)` places only the new sites of an unchanged scatter through a resumable sampler (103 ms to 8 ms). Results equal a fresh cook and lay-out for each kind of edit. Figures are in `docs/ROADMAP.md` next to the Field and BVH items. A time that no field of a scatter reads is not part of its identity (`reads_time`, derived from the ops), so a new `t` reuses its placements; and the static triangles are found in the `StaticCache` by the description of the static part before any is emitted (23-26 ms to 7-9 ms per update step at 1e3 to 1e5).
