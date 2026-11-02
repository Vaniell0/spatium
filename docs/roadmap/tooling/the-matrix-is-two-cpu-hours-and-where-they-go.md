---
date: 2026-10-08
title: The matrix is two CPU hours, and where they go
stage: 1
files: [scripts/gen_connectivity.py]
open:
  - the 24 heaviest cells (three derived spaces over fifty digits, eight probes each) are about a quarter of the CPU and are the only reason one cell takes minutes; held to the nightly run alone they would cut a pull request's matrix by that share, but the double cells of those spaces are graded against the fifty-digit signatures, so a run that skips them needs the committed signatures to grade by
  - the tolerance the derived charts integrate to over fifty digits is far tighter than the L2 check needs (sqrt(eps) times 64); a looser one would shorten those cells, and it is not measured
  - the split of the CPU by scalar and by space is not measured, only the twelve slowest cells
  - the run is CPU bound, not bound by its longest cell (7302 s over 12 threads is 608 s, the run took 632 s), so sharding it over jobs would shorten the wall time and spend the same CPU
---

`gen_connectivity.py --check` compiles and runs 1168 cells, each its own translation unit, and a pull request waits for all of it. Measured now with `--times` on 12 threads: 7302 seconds of cell time, 632 seconds of wall, so the run is bound by the CPU and the order of the cells cannot help beyond starting the longest first (it does that now: the slowest cell is `Real50 S2Derived FrechetMean`, 489 s, then `ExpLog` 267 s, `VerifyExpLog` 150 s, `MetricAxioms` 140 s, `DerivedDistance` 135 s; the twelve slowest are 1630 s, a fifth of everything).

What changed is the part that is the same in every cell. Each cell compiled `probes.hpp`, which includes the spaces, the algebra and Boost.Multiprecision, from scratch; it is now compiled once to a precompiled header with the flags of the cells and put beside their sources, where the compiler looks first. A cell over `double` takes 0.56 s instead of 1.80 s and one over `Real50` 2.25 s instead of 3.35 s, with the same output, and the regenerated `docs/connectivity.md` and `expected.hpp` are byte for byte the committed ones (`--check` exits 0). The saving is a sixth of the CPU, not more, because most of what a cell costs is instantiating its templates and, for the heavy ones, running them.
