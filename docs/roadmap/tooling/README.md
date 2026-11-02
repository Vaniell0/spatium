# Tooling

## Where it stands

How the repository keeps itself true: CI with path-filtered jobs, the generators and their `--check` modes (connectivity, capabilities, consumers, dependency graph, roadmap), the C++23 modules build, and the conventions in `docs/conventions.md`.

## Now

- the roadmap by area (this change)

## Later

- a map of the architecture as one page: every layer and how they connect

The entries of this area, by date, and what they did not do: `python3 scripts/gen_roadmap.py --area tooling`.

## Reasoning and earlier history

### Performance foundation

- Expression templates for Vec: VecLike concept, lazy +/-/*/÷, zero intermediates
- SIMD infrastructure (SSE2/4.1/AVX2): `if consteval` dispatch, stubs for missing ISA
- BVH<Shape> spatial index: SAH build (binned, 12 bins), stack traversal. ray_cast 57-253x vs brute-force
- Google Benchmark suite: 4 files (vec, intersection, mesh, bvh)

### API stability, CI, docs truing-up, C++23 modules folded into main (2026-08-28)

**Redefined 2026-08-28.** The original "production release" scope (mesh simplification, remeshing, PLY) measured readiness as geometry-processing feature parity with CGAL/libigl — a different, older idea of readiness than what the project actually became. Those are real, valuable features, just not what should gate declaring the API stable; moved to Backlog → Mesh processing. What actually gated it: API stability going forward, honest documentation, and something real to show, not a feature checklist against a competitor.

- Concept hierarchy + core public API considered stable as of this point: breaking changes from here on are a deliberate decision, not a "we'll circle back" note.
- CI (GitHub Actions): nix-build, bare-build, dependency-graph-freshness jobs — see "Code and build conventions, dependency graph" below.
- Documentation trued up against the real code, not just written once and left to drift: `docs/conventions.md` (new), `docs/architecture.md`'s Module Map and dependency claims corrected against a real generated graph, stale directory-layout/namespace claims in `CLAUDE.md` fixed.
- `exp/modules` merge strategy decided: folded into `main`'s `CMakeLists.txt` behind the hybrid `SPATIUM_USE_MODULES` flag — see "C++23 modules migration" below.

*Still open from this point: a real hero visual (the production GR black-hole render) — see Backlog → GPU rendering.*

### C++23 modules migration (2026-04-24)

7 phases bringing 11 named modules online (`spatium.core`, `spatium.algebra`, `spatium.spaces`, `spatium.mesh`, `spatium.geometry`, `spatium.point`, `spatium.spatial`, `spatium.discrete`, `spatium.io`, `spatium.physics`, `spatium.std` umbrella). Dual-mode `SPATIUM_EXPORT` macro, `import std.compat;`, CMake 3.28 `FILE_SET`, `spatium_module()` helper. Hybrid build: modules optional via `-DSPATIUM_USE_MODULES=ON`, header-tree stays the default so legacy consumers are unchanged. Header-units deferred until gcc 16.

This later went stale (never wired into CI, missing partitions for everything added to the header tree afterward) and needed real catch-up work — see "C++23 modules: caught up to a working state" below.

### Architecture audit (2026-08-26)

Docs truing-up pass: `architecture.md` concept-hierarchy fix + a new "header-only spine, and three principled exceptions" section, `CHANGELOG.md` deleted (redundant with this file), the "Block A/B/C/D" physics-milestone naming stripped from source comments and — completed in a follow-up pass — from `ROADMAP.md`/`api-reference.md`/`concept-driven-physics.md` as well. `Result<T>`/camera-module dedup deferred.

**Corrected 2026-09-15, twice over.** This entry used to end by claiming `cloth_sphere_probe.cpp` had been removed as a superseded investigation tool. It had not: the file is in `examples/`, no commit ever deletes it, and two other places in these docs describe it as live (`concept-driven-physics.md` §6.3, and the intrinsic-cloth entry in Backlog). The audit wrote down an intention as an accomplished fact. The file stays; no removal is planned.

And the naming sweep the same entry reports as complete was not complete either — `concept-driven-physics.md` §4.2 still listed `test_block_b_finish.cpp`, `test_block_c_close.cpp` and `test_block_d_start.cpp` three weeks later, as *file names*, which read as checkable references and were not. Fixed 2026-09-15 along with that document's stale test and LOC counts. The pattern worth keeping: in every one of those cases the document's *ideas* held up unchanged, and only its checkable facts rotted — which is an argument for deriving such facts rather than typing them, not for proofreading harder.

### Code and build conventions, dependency graph (2026-08-28)

A second truing-up pass, prompted by the same "code drifts from doc, doc drifts from code" pattern the 2026-08-26 architecture audit found, this time going into the code itself rather than just the docs describing it: `docs/conventions.md` (new) names the namespace/subdivision/error-handling/`measure()`-alias/doc-comment rules that were previously followed inconsistently per-domain, not written down anywhere. Applied: `algebra/` moved to `inline namespace algebra` (matching every other domain's per-namespace convention while staying bare-`spatium::`-compatible for the rest of the tree's existing unqualified usage — hit and fixed a real GCC constraint along the way, that `inline` must appear on the *first* reachable declaration of a namespace in a translation unit); `physics/atomic/` split out from `physics/`'s flat top level; `geometry/triangle.hpp`/`circle.hpp`'s `measure()`/`area()` alias inversion fixed (the project's own documented rule, violated in 2 of 3 sampled files). Build: `cmake/SpatiumTarget.cmake`'s `spatium_add_example()` replaces 16 copy-pasted `examples/CMakeLists.txt` blocks; `CMakePresets.json` gives 8 named, reproducible presets for what previously existed only as 10 undocumented `build-*/` directories. `scripts/gen_dependency_graph.py` generates `docs/dependency-graph.dot` from real `#include` edges (aggregated to domain level), checked in CI (`.github/workflows/ci.yml`'s `dependency-graph` job) so `architecture.md`'s Module Map can't silently drift from the real dependency graph again — confirmed real cross-domain edges the prior narrative description missed (`spaces`↔`mesh` mutual dependency, `discrete`→`geometry`/`algebra`, `geometry`→`spaces`). All 793 tests green throughout, phase by phase.

### C++23 modules: caught up to a working state (2026-08-28)

The 2026-04-24 migration (see above) had gone stale — this closes the gap rather than deleting the subsystem, since the missing partitions were genuinely mechanical, not a design problem.

- **New algebra partitions**: `dual`, `ode`, `linear_solve`, `calculus` — added to `modules/algebra/`, registered in `algebra.cppm`'s export list and `CMakeLists.txt`'s `PARTITIONS`. `eigen_interop.hpp` stays deliberately header-only (SSE intrinsics from `Eigen/Core` conflict with `:vec_simd`'s BMI — this was already documented, not new).
- **`physics/relativity/*` folded into `spatium.physics`**: `schwarzschild.hpp`, `kerr.hpp`, `geodesic.hpp`, `accretion_disk.hpp` added to `physics.cppm`'s single-TU `#include` block (no cross-includes among the four, any order works; `geodesic.hpp`'s `Dual<T>`/`ode`/`linear_solve` needs are covered by the algebra partitions above). `physics/atomic/`'s own stale include paths (left over from its 2026-08-28 split out of `physics/`'s flat top level) were also fixed in passing — a real, separate bug: `modules/physics.cppm` still referenced the pre-split `physics/orbital.hpp`-style paths, breaking the modules build specifically since the header-tree build never exercises that file. Drive-by cleanup: the "Block A/B/C/D" milestone-naming comments the 2026-08-26 architecture audit was supposed to have stripped everywhere had survived inside this one file; removed.
- **`spatium.render`: new module**, `modules/render.cppm` — single partition covering `camera.hpp`, `color.hpp`, `spectral.hpp`, `parallel_for_rows.hpp`, `supersample.hpp`, `sky.hpp`. `write_image.hpp` deliberately stays header-only: its whole contract is "exactly one translation unit defines `STB_IMAGE_WRITE_IMPLEMENTATION` before including it," which is a link-time property, not a module one — folding it into the module would bake `stb_image_write`'s externally-linked function bodies into `spatium_render_module` itself, and any TU that also defines the macro (every example does, so does `tests/test_render.cpp`) and links against the module in the same binary would hit duplicate-symbol errors. Same category of documented exclusion as `eigen_interop.hpp` and `mesh/dec.hpp`/`mesh/primitives.hpp` — a real technical call, not neglect.
- **Real GCC modules-ts bug found and fixed, not routed around**: linking a consumer of `spatium.render` against `Sky`/`make_starfield()` failed with `undefined reference` to `std::vector<int, allocator<int>>`'s special members (`_Vector_impl`'s destructor, the default constructor, the move constructor) — reproduced from a clean isolation (removing `Sky` usage from the test made the link succeed; adding real `#include <vector>` in both the producing module and the consuming TU did not fix it). Root cause narrowed to `Sky::star_buckets`'s element type specifically: switching `std::vector<std::vector<int>>` to `std::vector<std::vector<std::uint32_t>>` (matching the index-type convention already used elsewhere in the codebase — `BVH::Hit::index`, mesh vertex/face ids — not a workaround invented for this bug) made it link cleanly in both the modules and header-tree builds. A genuine, narrow GCC/libstdc++-modules-ts limitation with plain `int` crossing a module boundary through a nested container inside an inlined function, not an application-level mistake.
- **CI**: `.github/workflows/ci.yml` gained a `modules-build` job (configure/build/test via the existing `modules` CMake preset) — the gap ROADMAP had flagged as the actual reason this was allowed to rot the first time.
- **Verified**: 825/825 tests green on the `modules` preset (`build-modules/`, fresh configure) and 825/825 on the existing `build/` tree (which already had `SPATIUM_USE_MODULES=ON` cached from earlier work). `scripts/gen_dependency_graph.py --check` and `scripts/check_claude_md_layout.py` both still pass (neither script inspects `modules/`, so unaffected either way, checked to be sure).

### Copy-constructibility, leaking upward from `std::function`

Named 2026-09-15 as one item, after finding the third instance. Each was
recorded as a code comment at its own site, which is how a family of
related constraints reads as three unrelated footnotes.

`std::function` requires its callable to be copy-constructible, and the
requirement propagates to everything built on top:

| the slot | requires copyable | consequence |
|---|---|---|
| `ParametricSurface::ParamFn` | the `(u,v)→R³` map | the root of the chain |
| `io::build::ScalarField` (thickness) | the thickness callable | a thickness flows into `offset_surface()` and so into a `ParamFn`; it cannot own move-only state |
| `TraceNode::exact` (`std::any`) | the exact shape | a user shape owning a device handle or a `unique_ptr` cannot be recorded |

`PointField` escaped it by moving to `std::move_only_function` (PR #26),
which is why the donut demo's 19 800 motion closures can share one
`PerlinNoise` through a `shared_ptr` instead of copying 9.7 MB of
identical tables. The others did not, and the reason is the same in each
case: the type underneath insists on copying what it holds.

**Reformulated 2026-09-15: the goal is sharing, not move-only.** The
original plan for this item was `ParamFn` -> `move_only_function`, making
`ParametricSurface` move-only. That is now withdrawn, for two reasons
that only became visible once fields had a shape.

It aimed at the wrong target. The measured problem was never "a hook
needs to own a `unique_ptr`" -- nothing in the tree has move-only state
today. It was 19 800 copies of a 512-byte table, and `move_only_function`
does not fix that: each closure would still carry its own table, moved
rather than copied. What fixes it is **sharing the state**, which
`donut_demo.cpp:639` already does by hand with one
`std::make_shared<const PerlinNoise>`. The "hand-rolled shared_ptr dance"
this item feared turns out to be two lines, already written, and working.

And it pulled against the rest of the series. A `Field` owns its pool and
is copyable, which is what makes `TraceNode` copyable, which is what lets
`cook()` copy a trace rather than consume it. Making `ParametricSurface`
move-only would take that back. Worth noticing that the choice is now
genuinely a choice: "rebuild means a new trace, the old one discarded" is
today forced by `PointField` being move-only, and after this item it
becomes a decision to make on its merits rather than a constraint to obey.

So what remains here is:

- **the full `EdgeRule` set** -- `RoundCap(radius)`, `FlatCap`,
  `ExtendTo(node)`, which need fields that carry data and so depend on
  the item above;
- **making heavy shared state easy** -- `PerlinNoise` holding a
  `shared_ptr<const Table>` so a copy of the noise is a copy of a
  pointer, 8 bytes rather than 512, shared by construction rather than by
  the user remembering to. The same pattern for any future heavy payload
  in a `Field` leaf. Measured already: the report's payload across 19 800
  leaves is 10 137 600 B by value against 158 400 B by reference.
- **`ParametricSurface` stays copyable.** Not a compromise: copyable
  surface and shared heavy state are both right, and they work together.

One thing to know rather than solve, so it is not rediscovered: a
`shared_ptr<const T>` copy is an atomic refcount. Copying fields in a
loop over tens of thousands of them would count atomics, and if that ever
lands in the frame path the answer is a non-atomic intrusive pointer for
single-threaded use -- which is exactly the hand-rolled work this item
originally flinched from. Not now; it is not in a hot path today.

Worth stating because it is the obvious wrong turn: **the escape is
`move_only_function`, not `copyable_function`.** The latter requires a
copy-constructible target exactly as `std::function` does. Reaching for
it here would be swapping one copyable erasure for another and changing
nothing.

Nothing in the tree hits any of these today — surfaces and shapes are
value types, thickness fields capture numbers. What makes it worth an
entry rather than three comments is that the failure is always a
template error from `<functional>` or `<any>` that says nothing about
this library's design, and that **the same question arrives a fourth
time with structural fields**: a `Field` whose leaves may own state will
meet exactly this wall.

- **[want]** Lift `ParametricSurface::ParamFn` off `std::function`. It is
  the root: fix it and `ScalarField` follows.

  **`std::copyable_function` is not the fix, despite the name.** It
  requires a copy-constructible target exactly as `std::function` does --
  that is what "copyable" means there. What it actually repairs is a
  different hole: `std::function::operator()` is const but happily
  invokes a non-const target, and it drops the `target()`/RTTI surface.
  Worth having for those reasons; irrelevant to this one.

  The only type that relaxes the requirement is
  `std::move_only_function`, which is what `PointField` already uses --
  so the fix is to make `ParametricSurface` **move-only**, not to swap
  one copyable erasure for another. That is a real change rather than an
  `#if`: `offset_surface()` captures its base into a new closure,
  `Placed::surface()` hands back an `optional<ParametricSurface<T>>`,
  `resolve_surface()` returns by value, and several tests copy surfaces.
  All of those work under move, but each is a line that has to change.

  Measured here, GCC 15.3 under `-std=c++26`:
  `__cpp_lib_move_only_function` is `202110`, while
  `__cpp_lib_copyable_function` and `__cpp_lib_function_ref` are absent.
  So the tool this needs has been available all along; what is missing is
  the work, not the vocabulary.

  Templating `ParametricSurface` on its callable instead is tempting and
  wrong at this boundary: `resolve_surface()` has to return **one** type
  for nodes whose maps differ, so erasure there is a requirement rather
  than laziness.

- **[want]** Decide whether `TraceNode::exact` should accept a move-only
  shape. `std::any` cannot; a `unique_ptr<void, deleter>` plus a
  `std::type_index` can, at the cost of hand-rolling what `std::any`
  gives for free. Worth doing when a shape that needs it exists, not
  before — but worth *knowing* before someone discovers it as a
  compiler error.

### Interop / ecosystem

- **[dare]** Heat-method log map, CGAL-grade exact polyhedral geodesics (geometry-central and CGAL each cover one half of this; Spatium currently ships neither on top of Dijkstra/heat-distance).
- **[want]** Own Vec/Matrix creates impedance mismatch with the Eigen ecosystem — interop adapters beyond the current `to_eigen`/`from_eigen`/`eigen_view`. Previously gated on native SVD/eigendecomposition, which landed 2026-09-06 (see Backlog → Native math), so nothing blocks this now.

---

## Competitive Position

| vs | Their strength | How we differ |
|----|---------------|---------------|
| CGAL (~500K LOC) | Exact arithmetic, robust predicates, maturity | Concepts as extension point; sdk header-only |
| libigl / geometry-central | Heat method, mesh processing in R³, Eigen ecosystem | Surface as concept, not hardcoded type; built-in axiom verification |
| Sophus / manif (~3-5K LOC) | Analytical Jacobians, Ceres/autodiff integration | Generalized algebra (Group/Ring/Field), not just Lie groups |
| Eigen (~150K LOC) | SIMD, autodiff, entire scientific C++ ecosystem | Spaces as first-class (not just matrices), morphism pipes |
| glm (~30K LOC) | GLSL-compatible API | Algebraic structures, manifolds |
| geogram (~200K LOC) | Robust Voronoi/remeshing in Rⁿ | Manifold-generic (once heat method lands) |

**Unique selling point:** the only C++ library where a mathematical space is an extensible C++23 concept, not a hardcoded type. A user-defined `Surface` gets mesh subdivision, geodesics, parallel transport, and Vulkan rendering from a single `struct` + `static_assert`.

**Caveats:**
- Geodesics: Dijkstra (O(h) error on mesh edges) ships unconditionally; heat method (Crane 2013) ships with `SPATIUM_EIGEN=ON` via pre-factored `HeatSolver<S>` (see "Eigen interop, heat method" above). geometry-central also exposes the heat log map; CGAL adds exact MMP. See Backlog → Interop / ecosystem.
- Own Vec/Matrix creates impedance mismatch with the Eigen ecosystem — see Backlog → Interop / ecosystem.
