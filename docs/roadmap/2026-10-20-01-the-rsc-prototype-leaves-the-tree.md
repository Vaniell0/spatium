---
date: 2026-10-20
title: The RSC prototype leaves the tree
files: [include/spatium/core/concepts.hpp, include/spatium/physics/mechanics/surface_ccd.hpp, tests/test_surface_ccd.cpp, scripts/check_doc_file_refs.py]
open:
  - the engine itself, `rsc/core` (generic: jobs, verification, the ladder of organisations behind one interface and one artifact format) and `rsc/spatium` (the binding to this library), is not written
  - the first synthesis job, a chain of components over a substrate object and checked on every dimension N=1..8, is not chosen by measurement yet
  - RSC in CI (a scheduled workflow that opens draft PRs on `rsc/*` branches, and a check step in `ci.yml`) does not exist
---

The 2026-09 RSC prototype is gone from the tree: `rsc/` entire (the registry of `Op`s on `std::function`, seven domains, the `Dispatcher`, REINFORCE training, checkpoints and the embedded base, twenty tools, the CCD chooser tree), the eighteen `tests/test_rsc_*.cpp`, `SPATIUM_BUILD_RSC_TOOLS`, and the commit SHA and base checkpoint that the top-level `CMakeLists.txt` baked into the build. The last commit that has it is tagged `rsc-prototype-final`; anything it did is read, or rerun, from that tag. Its numbers were never in the code: they are in `docs/findings-2026-09.md` section 4 (the prototype audited, including why the reproducibility checks protected nothing) and in `docs/ROADMAP.md` under "RSC as search" (the collapse of scalar chains by a factor of 706, the CCD chains, the transfer between spaces), and they stay there as the baseline the next engine is measured against. The files those documents name are listed in `scripts/check_doc_file_refs.py` as removed.

What stays is the library's side of the seam: `Chooser<C, F>` in `core/concepts.hpp` and the overload of `first_contact` in `surface_ccd.hpp` that takes one. It is the choice point a synthesised artifact will enter through, and its coverage moved from the prototype's test of the generated tree into `tests/test_surface_ccd.cpp`, held with a hand-written chooser that branches on the features it is given. The prototype's other CCD tests (the label gate, the chain primitives, and the learned tree against the default's cost) went with the code they tested. The comments in `include/` and `examples/` that named `rsc/` paths now name the tag. No behaviour of the library changed.

`rsc/` is to be rebuilt as a synthesis engine inside this repository, with a boundary drawn in the tree rather than a separate repository: `rsc/core` does not include `<spatium/...>`, `rsc/spatium` binds it to the library, and the library does not depend on `rsc/`.
