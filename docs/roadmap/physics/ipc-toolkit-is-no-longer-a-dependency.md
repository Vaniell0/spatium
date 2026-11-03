---
date: 2026-10-21
title: ipc-toolkit is no longer a dependency; the comparison with it runs out of tree
files: [CMakeLists.txt, tests/CMakeLists.txt, examples/CMakeLists.txt, scripts/external/ipc_compare/CMakeLists.txt, scripts/external/ipc_compare/cloth_sphere_probe.cpp, scripts/check_doc_file_refs.py, README.md, CLAUDE.md, include/spatium/physics/mechanics/contact.hpp]
open:
  - scripts/external/ipc_compare/ is not built by CI and was not configured here (it fetches and builds ipc-toolkit's dependency tree); it is carried over unchanged and unverified
  - docs/concept-driven-physics.md and the older ROADMAP entries still describe SPATIUM_IPC_TOOLKIT as live; they are left for the planned rework of the physics docs
  - the IPC replacement the comparison is meant to measure (no discreteness, no ipc-toolkit) is not written
---

`SPATIUM_IPC_TOOLKIT`, its `FetchContent` block, `SPATIUM_HAS_IPC_TOOLKIT`, `tests/test_ipc_toolkit.cpp` and the gated `cloth_sphere_probe` example are removed from the build. `cloth_sphere_probe.cpp` moved to `scripts/external/ipc_compare/` with its own `CMakeLists.txt`, which fetches ipc-toolkit v1.6.0 itself, so the comparison with IPC stays possible without Spatium depending on it. The last in-tree state is the git tag `ipc-toolkit-final`.
