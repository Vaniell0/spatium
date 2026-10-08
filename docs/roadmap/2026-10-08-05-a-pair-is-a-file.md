---
date: 2026-10-08
title: A pair is a file, and the sweep is a target
stage: 1
files: [tests/symmetry/pair.hpp, tests/symmetry/common.hpp, tests/symmetry/sweep.cpp, tests/symmetry/pairs/ray_torus_precision.cpp, tests/symmetry/pairs/spd3_eigenvalues.cpp, tests/test_symmetry.cpp, tests/CMakeLists.txt]
open:
  - the sweep is run by hand; nothing schedules it, and it is the first step of the unattended loop the plan describes for RSC in CI
  - the pairs that live as single tests (member / ADL / derived, erased / typed, interpret / eval_into, host / device) are not files here yet
  - no pair declares the trusted base and the `uses` of its two paths, so the independence of a witness is not checked
  - a pair file that fails to compile stops the whole test binary, not only itself
---

The registry of pairs was one header of 540 lines that every new pair edited, which is the shape that makes any number of writers conflict at every merge. A pair is now a file in `tests/symmetry/pairs/`, nine of them, each with its own banner (what the two paths are, where the bound comes from, what the first run measured), and each registers itself at load with `Registrar`. `tests/CMakeLists.txt` globs the directory, so a new pair is a new file and nothing shared changes; each is its own translation unit, so one is compiled alone in seconds before anything is proposed. What two pairs share (the rays through the torus's box, a point from a draw of numbers) is `common.hpp`.

The sweep became a target. `ninja symmetry_sweep && ./build/tests/symmetry_sweep 300000 [fragment]` runs every registered pair on as many inputs as asked and exits 1 on any failure; the test takes each pair at its default size. A rare failure is found by volume (one ray in twenty thousand was off by 4e-4), and a scheduled run repeating this is the maintenance loop the RSC line describes for CI.
