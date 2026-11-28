---
date: 2026-11-28
title: the law enumerator is a generic engine in rsc/core, and spatium is its first adapter
files: [rsc/core/law_engine.hpp, rsc/spatium/rsc_laws.cpp]
open:
  - the enumeration of terms by size stays with each adapter (spatium's keeps one layer of scaling, ordered pairs, base-matched tangents); a signature-driven enumerator in the engine is for when a second adapter shows what they share
  - a second project's adapter (SQLite: its functions and operators by sort, spaces = the context an expression is evaluated in and the optimizer's switches) is the consumer that will show whether the six questions the engine asks of a signature are enough
---

What finds laws no longer knows what a term means. `rsc/core/law_engine.hpp` holds the values of a term on a space and their comparison at a tolerance, the space's interface (a batch staged, a term accepted), the world that evaluates a batch on every space at once and merges a term equal to one already kept, and the report: every pair equal or ordered on some spaces, a law dropped when it is an instance of a smaller one or follows from smaller ones by equality, congruence or a constant factor, the chains of what the closure derived, the JSON, and `check` reading laws back. It asks a signature six things: which operations are variables, constants, the product and the sum of the ordered sort, which sort is ordered, which terms are worth a law. `rsc/spatium/rsc_laws.cpp` is what spatium gives it: the sorts P, V@x, R, the operations of `core/access.hpp`, the 25 spaces of the matrix and the breaker, the terms of each size, the prefix form read back (884 lines before, 466 now, the engine 497). At size 7 its output is byte for byte what it was: 787 terms, 360 merged, 637 laws, 590 derivations; `rsc_laws.py check` finds 0 changes, and a false law (d(p,q) = d(p,r)) is defined everywhere and holds nowhere.
