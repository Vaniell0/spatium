---
date: 2026-10-29
title: the hot paths allocate nothing, held by a count of heap allocations
files: [tests/test_allocations.cpp, tests/CMakeLists.txt, include/spatium/algebra/groups/action.hpp, include/spatium/discrete/lattice.hpp]
open:
  - a field of more than 32 ops still allocates its scratch on every evaluation (the metric's ten fields already keep 256 on the stack); raising the bound costs every small field the zeroing of a larger array of Vec, so the bound waits for a measured caller
  - the counts are not yet a column of the connectivity matrix, and no time is measured: the benchmarks in benchmarks/ have no table across the project
  - the search and self-play of the chess series, and the synthesiser's loop, are not in the library yet and are not counted
---

`test_allocations.cpp` replaces the global operator new for the test binary with one that counts and calls malloc, and holds the number of heap allocations per call on the library's hot paths. A count is exact and does not move with the machine's load, so it can be held in CI where a time cannot. Zero for: the gradient through `Dual`, a 4x4 matrix product, `solve_direct` on 8x8, `integrate` with and without its witness, a field of up to 32 ops, `legal_moves` and `perft` on Kiwipete, `outcome`, `from_fen` and `to_san`. The counter is held to see an allocation it cannot be optimised out of.

`algebra::orbit` and `discrete::attacks` return vectors, and grew them one push at a time: 7 and 4 allocations a call. Each now reserves its bound (the group's order; one cell per leap and the longest axis per slide), so each is one allocation, held as at most one since an optimiser may drop one that does not escape.
