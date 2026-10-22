---
date: 2026-10-22
title: a pull-based reactive layer
files: [include/spatium/core/reactive.hpp, modules/core/reactive.cppm, tests/test_reactive.cpp]
open:
  - no in-tree consumer yet; only the tests exercise the header
---

`core/reactive.hpp` adds `Cell<T>` (an input whose version moves only when `set` changes the value) and `Memo<F>` (a lazy derived value). A memo records what its function read through a thread-local observer, caches the result, and recomputes on read only when a dependency's version moved. A memo whose new value equals the old one keeps its version, so its dependents do not rerun. Dependencies are rebuilt on every run, so conditional reads are re-tracked. There are no virtual calls and no allocation beyond each memo's dependency list.
