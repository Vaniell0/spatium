# Roadmap

## Vision

Computational geometry on arbitrary Riemannian manifolds. The only C++ library that unifies abstract space theory and practical geometry through concepts. Target: geodesic Voronoi/Delaunay, shortest paths, remeshing — on any Surface, not just Euclidean.

This document tracks what's done and what's open, by content and by date — not by version number. Version-gated milestones invite the failure mode they're meant to prevent: negotiating what counts toward which release, uneven development where some gates get stuffed and others starve, and a "two parallel tracks" split that has to be reconciled later — the same fragmentation branches/forks cause when they diverge instead of staying one line. Completed work below is dated where the date is known, named by content where it isn't. Open work is grouped by topic in **Backlog**; an item moves up into **Completed** the day it actually ships, not on a schedule.

---

# Completed

## Concept hierarchy, core geometry, mesh, viewer foundations

- Concept hierarchy: Set → TopologicalSpace → MetricSpace → ... → Surface (12 concepts)
- Concrete spaces: Euclidean<N>, Sphere<N>, Hyperbolic<N>, ProductSpace
- Geometry: 10 primitives, 11 intersect pairs, distance, contains, BoundedRegion + clip
- Boolean ops: intersection_region (Sutherland-Hodgman), difference/symdiff area
- Algebra: Vec, Matrix, AffineTransform, SO(3), SE(3), Group/Ring/Field concepts
- Mesh: Mesh<Surface>, midpoint subdivision + project, LodChain, icosahedron/tetrahedron
- Viewer: Vulkan 1.3, orbit camera, per-face coloring, edge overlay, dynamic resize, nix run
- Discrete: FiniteSet<T> (∪ ∩ ∖ △ ⊆, power_set, cartesian), GeometricSet, SVG Venn
- IO: Table, Svg (polygon fill), std::format for all types
- Extras: Point<Space>, Morphism pipes, Real50/Real100, axiom verification
- 295 tests, 2 audits (0 CRIT/HIGH), ~10K LOC

## Performance foundation

- Expression templates for Vec: VecLike concept, lazy +/-/*/÷, zero intermediates
- SIMD infrastructure (SSE2/4.1/AVX2): `if consteval` dispatch, stubs for missing ISA
- BVH<Shape> spatial index: SAH build (binned, 12 bins), stack traversal. ray_cast 57-253x vs brute-force
- Google Benchmark suite: 4 files (vec, intersection, mesh, bvh)

## Geodesic algorithms

- MeshTopology: edges, vertex neighbors, face-edge map, boundary detection
- Geodesic distance field: Dijkstra with space.distance() — works on any Surface
- Shortest path: Dijkstra + predecessor trace
- Parallel transport: Schild's ladder (exp/log only)
- Geodesic Voronoi: multi-source Dijkstra + label propagation + face_labels

## Surfaces & usability

- ParametricSurface<T>: f(u,v)→R³ with auto metric. Torus, cylinder, cone, Mobius factories
- ImplicitSurface<T>: F(x,y,z)=0 with auto gradient. Sphere, torus, gyroid factories
- Marching cubes: uniform grid extraction → Mesh
- Mesh operations: merge, transform, flip_normals, compute_face/vertex_normals, centered
- Mesh primitives: grid_mesh, uv_sphere_mesh, box_mesh
- OBJ I/O: load_obj (polygon fan triangulation, negative indices), save_obj
- STL I/O: load_stl (auto-detect binary/ASCII), save_stl (binary)
- Quaternion<T>: from_axis_angle, from_matrix, to_matrix, rotate, slerp, inverse

## Audit fixes, multi-mesh viewer, atom demo

- Marching cubes: Paul Bourke 256-entry tri_table (replaces centroid-fan)
- BVH stack: dynamic vector with reserve(64) (replaces fixed array<64>)
- MeshTopology: shared_ptr<const Mesh> ownership (replaces dangling reference)
- Matrix::inverse() → Result<Matrix> (was silent identity on singular)
- Quaternion::inverse() → Result<Quaternion> (was division by zero)
- identity() morphism: self-inverse (was nullopt)
- OBJ loader: std::from_chars (was std::stoi with exceptions)
- Viewer: multi-mesh rendering, per-mesh color push constants, alpha blending
- Physics: AtomModel, orbital wavefunctions, element database (118 elements), atom SVG
- atom_demo: CLI + console + Vulkan viewer + SVG export + shell toggle
- 437 tests, ~15K LOC

## ImGui viewer, Bohr model

- ImGui viewer: element picker, shell toggles, view mode controls
- Bohr model electron orbital visualization
- Dirty buffers: lazy GPU upload on mesh change

## DSL, concepts, Complex<T>, ray-quadric

- Vector UDLs: `3.0_x + 2.0_y + 1.0_z` via expression templates
- Uniform distance(): 7 new overloads (Box-Box, Circle-Circle, Triangle-Triangle, Polygon-Polygon, etc.)
- Polygon boolean operators: `operator&` (intersection), `operator-` (difference), `operator+` (union), `union_of()`
- Generalized pipe-unwrap: `Result<Vec> | Transform`, `Result<Point> | Morphism`, chains
- Lazy transform chains: `lazy(translate) * lazy(rotate)` with `.apply()` / `.collapse()`
- Concept constraints: Shape/Measurable requires in boolean.hpp
- Generic algebra functions: `power()` (Group), `commutator()` (Group), `adjoint()` (LieGroup), `poly_eval()` (Ring)
- Complex<T>: full arithmetic, conjugate, magnitude, phase, from_polar, sqrt, cbrt, constexpr
- Polynomial solvers: solve_quadratic, solve_cubic (Cardano), solve_quartic (Ferrari) → Complex roots
- Analytical ray-quadric: Quadric<T> (sphere/cylinder/cone/ellipsoid), ray_quadric (hits), ray_quadric_proximity (miss via imaginary part)
- 530 tests, ~18K LOC

## Eigen interop, heat method, operator fix

- Eigen3 interop (optional, SPATIUM_EIGEN=ON): to_eigen/from_eigen, eigen_view (zero-copy Map)
- Heat method geodesics (Crane 2013): O(h^2) via Eigen::SimplicialLDLT. HeatSolver pre-factored.
- Discrete differential operators: cotangent Laplacian, lumped mass matrix, face gradients, integrated divergence
- GeodesicMethod enum: Dijkstra | Heat, backward-compatible default
- operator| disambiguation: Polygon union → operator+/union_of(), | reserved for intersection/pipe

## Analytical-render dispatcher, benchmark suite, table fix (2026-04-23)

- **`examples/primitives_demo.cpp`** — multi-scene Vulkan dispatcher: `--scene primitives` (unified sphere/box/torus/cylinder/cone/ellipsoid/triangle with BVH raycast hit-cloud + console report), `--scene torus` (Clifford torus S³→R³), `--scene klein` (animated Klein bottle R⁴→R³). Absorbs and removes the former `viewer_demo`, `torus_demo` and `klein_demo`.
- **`benchmarks/bench_raycast.cpp`** — BVH (60–71 ns) vs brute (8.7 µs – 519 µs) vs analytical `ray_quadric` (23–29 ns) across 320 / 5 120 / 20 480 triangles
- **Table renderer fix** (`include/spatium/io/table.hpp`) — separator width off-by-one (`w+1` → `w+2`); console tables now render with aligned borders. (Path corrected 2026-09-15: this said `src/io/table.cpp`, which has never existed in this tree — `src/` holds three files and none of them is the table renderer, which is header-only like the rest of `io/`. Found by `scripts/check_doc_file_refs.py`, not by hand.)
- **Eigen integration decision** — remain optional via `SPATIUM_EIGEN=ON`; heat method + differential + interop stay isolated. Vec/Matrix not migrated to Eigen; `Eigen::Ref<…>` planned for public APIs only once dense QR/SVD lands (see Backlog → Native math)
- ~~CPU raytracer: Quadric → viewer texture (analytical render without mesh)~~ done: `ray_quadric`/`ray_quadric_proximity` in `geometry/ray_surface.hpp`, exercised by `examples/parametric_analytical_demo.cpp`'s `glow_sphere.png`; ships PNG output rather than a live viewer texture, which covers the same "analytical render without mesh" goal for the gallery use case
- ~~Ray-ParametricSurface: Newton UV iteration for arbitrary f(u,v) surfaces~~ done: `geometry/ray_parametric.hpp` — no tessellation, no BVH, Newton-solves `S(u,v) = o + t·d`; exercised by `examples/parametric_analytical_demo.cpp` producing `klein_analytical.png`/`mobius_analytical.png`/`bumpy_analytical.png`/`parametric_gallery.png`
- ~~Ray-quartic: torus intersection via solve_quartic~~ done: `geometry/ray_surface.hpp`'s `Torus<T>` + `ray_torus()`/`ray_torus_proximity()`, exercised by `primitives_demo`'s analytic-torus raycast row and `torus_analytical.png` — this had drifted stale in Backlog → Analytical rendering (fixed 2026-08-31)

## API stability, CI, docs truing-up, C++23 modules folded into main (2026-08-28)

**Redefined 2026-08-28.** The original "production release" scope (mesh simplification, remeshing, PLY) measured readiness as geometry-processing feature parity with CGAL/libigl — a different, older idea of readiness than what the project actually became. Those are real, valuable features, just not what should gate declaring the API stable; moved to Backlog → Mesh processing. What actually gated it: API stability going forward, honest documentation, and something real to show, not a feature checklist against a competitor.

- Concept hierarchy + core public API considered stable as of this point: breaking changes from here on are a deliberate decision, not a "we'll circle back" note.
- CI (GitHub Actions): nix-build, bare-build, dependency-graph-freshness jobs — see "Code and build conventions, dependency graph" below.
- Documentation trued up against the real code, not just written once and left to drift: `docs/conventions.md` (new), `docs/architecture.md`'s Module Map and dependency claims corrected against a real generated graph, stale directory-layout/namespace claims in `CLAUDE.md` fixed.
- `exp/modules` merge strategy decided: folded into `main`'s `CMakeLists.txt` behind the hybrid `SPATIUM_USE_MODULES` flag — see "C++23 modules migration" below.

*Still open from this point: a real hero visual (the production GR black-hole render) — see Backlog → GPU rendering.*

## C++23 modules migration (2026-04-24)

7 phases bringing 11 named modules online (`spatium.core`, `spatium.algebra`, `spatium.spaces`, `spatium.mesh`, `spatium.geometry`, `spatium.point`, `spatium.spatial`, `spatium.discrete`, `spatium.io`, `spatium.physics`, `spatium.std` umbrella). Dual-mode `SPATIUM_EXPORT` macro, `import std.compat;`, CMake 3.28 `FILE_SET`, `spatium_module()` helper. Hybrid build: modules optional via `-DSPATIUM_USE_MODULES=ON`, header-tree stays the default so legacy consumers are unchanged. Header-units deferred until gcc 16.

This later went stale (never wired into CI, missing partitions for everything added to the header tree afterward) and needed real catch-up work — see "C++23 modules: caught up to a working state" below.

## Physics v2: Mechanics, Symplectic/DEC/Lie-group, Variational, Contact (2026-04-24)

Top-level abstraction: **discrete Lagrangian variational integrator on Riemannian manifold with DEC forms** — not the Bullet-style "RigidBody + Force + Constraint" tree.

- **Mechanics foundations.** Compile-time SI via `Quantity<M,L,T,I,K,N,J>` + `std::ratio`, `PointMass<N,T>`, `RigidBody<N,T>`, composable forces (`UniformGravity`, `PointGravity`, `Spring`, `Damper`), integrators (`euler_step`, `semi_implicit_euler_step`, `verlet_step`, `rk4_step`). Kepler orbit over 10 periods under Verlet: energy + L conserved to 1e-3.
- **Symplectic + DEC + Lie-group.** `SymplecticManifold<S>` concept + `CotangentBundle<M>` wrapper + `verify_symplecticity_drift`, `Yoshida4` composition integrator (KAM 1000-period proof), `lie_rkmk4_cf_step` (commutator-free 4-th order RKMK, Celledoni-Marthinsen-Owren 2003), `PointOnManifold<M>` + analytical `geodesic_step<Sphere<N>>`. DEC primitives in `mesh/dec.hpp`: `Form0/Form1/Form2`, `exterior_derivative_0/1`, `hodge_star_0/1/2`, `laplace_beltrami_dec`.
- **Variational + continuum.** `DiscreteLagrangian<LD,N,T>` concept, `SeparableMidpointLagrangian`, `variational_step_separable` (kick-drift-kick Verlet as closed-form DEL). LGVI on SO(3) via Cayley 1-cut (orientation + Noether + Casimir exact; energy drift ~3 % is a Cayley-1-cut limitation documented as follow-up). `DecHeatSolver` (backward-Euler `(M + dt·L) φ = M·φ_n` pre-factored via Eigen `SimplicialLDLT`). `DeformationMap<MFrom,MTo>`, `deformation_gradient` (FD), `right_cauchy_green`, `green_strain`, `StrainEnergy<W,T,N>` concept, `SaintVenantKirchhoff`.
- **Contact (partial).** `ipc_barrier(d, d̂)` = -(d-d̂)²·log(d/d̂) with gradient + Hessian + default stiffness, `XpbdParticle`/`XpbdDistanceConstraint` + Gauss-Seidel projection + `xpbd_step`, `build_distance_constraints`, `build_bending_distance_constraints`. **`ContactSurface<S,T>` concept + free `point_to(p, surf)` overload set** for `Sphere<2>` and `Torus` (closed-form) + generic `Surface` fallback that picks up `ParametricSurface`, `ImplicitSurface`, and any user type satisfying the Surface concept. `ipc_contact_energy`/`ipc_contact_force_on` concept-constrained helpers. Full implicit `cloth-on-obstacle` pipeline (and Klein-bottle USP demo) deferred pending ipc-toolkit integration — explicit position-based contact ceiling was hit with overhanging corners yanking cloth through the obstacle faster than the contact band reacts. (Continued: "Contact physics build dependency", "Implicit-contact Newton solver" below.)

### Tests + benchmarks

- **673 test cases, 6 549 Catch2 assertions**, all passing on modules ON build (GCC 15 + Eigen 3.4), at this point in the project's history.
- `benchmarks/bench_narrow_phase.cpp`: `point_to_sphere` 3.7 ns/op, `point_to_torus` 17.7 ns/op, `point_to` on a parametric surface (Newton-on-UV via `Surface::project`) 5 245 ns/op, full IPC pipeline (query + energy + force) 12.6–25.5 ns/op for analytical surfaces. Order-of-magnitude faster than GJK (~100-500 ns/op) on closed-form paths.

### Reference specification

- `docs/concept-driven-physics.md` — a rigorous specification of the concept hierarchy (Set → Manifold → Surface + `SymplecticManifold`, `DiscreteLagrangian`, `ContactSurface`) and the dispatch pattern. Fixes the current status matrix and limitations.

## RSC: a trained dispatcher on top of Spatium

`rsc/` — a trained dispatcher deciding *which* Spatium call to make and with *what* parameters, not a learned recurrence pretending to be a computer: Spatium stays the exact substrate. Minibatched policy-gradient (REINFORCE) training, hand-derived backward. Domain pipeline, ordered by how proven the underlying mechanism is (named by content, never by position):

- **Precision-critical dispatch** (`precision_ops.hpp`) — `solve_cubic` f64 vs. Real50 on casus-irreducibilis cubics; ~0.99 held-out accuracy.
- **Geodesic/mesh dispatch** (`geodesic_task.hpp`) — Dijkstra vs. Heat method; proven mechanism, partial accuracy (~0.536 on class-balanced sampling — real signal, not yet near-ceiling).
- **Root-finding dispatch** (`rootfind_ops.hpp`) — Newton vs. bisection around `f(x)=x^3-a`'s `f'(x)=0` inflection.
- **Cauchy/IVP dispatch** (`ode_ops.hpp`) — Euler vs. RK4 across Decay/Oscillator/CircularOrbit families.
- **Mesh-strategy dispatch** (`mesh_ops.hpp`) — uniform vs. anisotropy-adapted UV tessellation; plateaus ~0.77-0.79 over a ~0.26 baseline.
- **General linear-solve dispatch** (`linear_ops.hpp`) — Jacobi iteration vs. direct Gaussian elimination for `Ax=b`.
- **Rigid-body integrator dispatch** (`integrator_ops.hpp`) — five real `PointMass` steppers (`euler`/`semi_implicit_euler`/`verlet`/`rk4`/`yoshida4`), already present in `physics/mechanics/integrator.hpp` with no new algorithm work needed, across three closed-form test families (UniformGravity, Spring, PointGravity). Measured dispatch mix (not uniform): Verlet ~63%, RK4 ~24%, Euler ~8%, semi-implicit Euler ~4%, Yoshida4 ~1% — every candidate wins real problems. 4.9% → 84.5% held-out accuracy over 2500 REINFORCE updates.
- **Contact physics / soft bodies** — investigated and found structurally blocked: a 20-configuration compliance/substep sweep of explicit XPBD contact against a heavy cloth overhang failed on all 20 (stretch ≥23%, often exploding). See "Contact physics build dependency" below for the resolution in progress.
- **Real-time control of complex dynamics** (illustrative: an underwater drone) — see Backlog → Contact physics / RSC; sequenced last, depends on the base+custom deployment split actually working.

New Spatium primitives built specifically to support this pipeline: `algebra/dual.hpp` (forward-mode autodiff), `algebra/calculus.hpp` (`gradient`/`integrate`/`minimize`), `algebra/ode.hpp` (generic `euler_step`/`rk4_step`/`integrate_fixed` over `Vec<T,N>`), `algebra/linear_solve.hpp` (`solve_direct`/`solve_jacobi`/`diagonal_dominance_ratio`), `spaces/parametric.hpp`'s `parametrization_anisotropy`/`normal_at`. Full per-domain writeups, measured accuracy numbers, and postmortems (each domain surfaced at least one real bug during training, not just tuning) live in `rsc/README.md`.

## General relativity + analytical rendering

- `physics/relativity/` — metric-agnostic geodesic integration: `schwarzschild.hpp` and `kerr.hpp` (both metrics as templated, `Dual<T>`-substitutable callables), `geodesic.hpp` (exact Christoffel symbols via `Dual<T>` partials, Killing-vector conserved quantities), `accretion_disk.hpp` (Schwarzschild thin-disk redshift; Kerr's own BPT-1972 equatorial-orbit/ISCO/photon-orbit/redshift formulas).
- `render/` — `supersample_pixel()`, NxN jittered-grid antialiasing for any CPU raytracer built on Spatium's ray-surface primitives.
- `examples/tumbling_body_demo.cpp` — LGVI-integrated Dzhanibekov-effect (tumbling-body instability) demo.
- `examples/blackhole_gr_demo.cpp`, `wormhole_demo.cpp`, `geodesic_curvature_grid_demo.cpp` — gallery renders built on the above.

## GPU rendering (CUDA) — kernels built and cross-validated

`gpu/` exists and works, and is genuinely necessary, not optional polish: CPU-only was directly measured and extrapolated to ~20 days for the actual target render (1920×1080, ~750 frames, Kerr flyby) — not a guess, a real number from a real benchmark, after Mojo was investigated and rejected in favor of plain CUDA C++ (T4 fp64 throughput is 1/32 of fp32, which settles the language question regardless of Mojo-vs-CUDA specifics) and an fp32 precision gate was verified to pass cleanly across 22 cross-validated cases. `geodesic_kernel.cu` (batch Schwarzschild/Kerr null-geodesic tracing), `schwarzschild_render_kernel.cu`/`kerr_render_kernel.cu` (full per-ray volumetric disk emission-absorption, ~2K LOC total), `render_4k_frame.cpp` (the actual production render driver) — cross-checked against the CPU path per-ray at 160×90 (`gpu/verify_cuda_render.cpp`, abs/rel diff against the same real library calls the CPU demo uses), not just visually eyeballed.

Closed-form Christoffel symbols were hand-derived for exactly the two metrics the GPU path needs, rather than sharing the generic `Dual<T>`-templated CPU headers directly as `__host__ __device__` code — that was considered and deliberately rejected, since `Dual<T>`'s general autodiff machinery is branch/indirection-heavy in a way that's specifically bad for SIMT throughput. See Backlog → GPU rendering for what's still open (landing the production render, and the one remaining mechanization gap in how those closed-form symbols get into the header).

## Architecture audit (2026-08-26)

Docs truing-up pass: `architecture.md` concept-hierarchy fix + a new "header-only spine, and three principled exceptions" section, `CHANGELOG.md` deleted (redundant with this file), the "Block A/B/C/D" physics-milestone naming stripped from source comments and — completed in a follow-up pass — from `ROADMAP.md`/`api-reference.md`/`concept-driven-physics.md` as well. `Result<T>`/camera-module dedup deferred.

**Corrected 2026-09-15, twice over.** This entry used to end by claiming `cloth_sphere_probe.cpp` had been removed as a superseded investigation tool. It had not: the file is in `examples/`, no commit ever deletes it, and two other places in these docs describe it as live (`concept-driven-physics.md` §6.3, and the intrinsic-cloth entry in Backlog). The audit wrote down an intention as an accomplished fact. The file stays; no removal is planned.

And the naming sweep the same entry reports as complete was not complete either — `concept-driven-physics.md` §4.2 still listed `test_block_b_finish.cpp`, `test_block_c_close.cpp` and `test_block_d_start.cpp` three weeks later, as *file names*, which read as checkable references and were not. Fixed 2026-09-15 along with that document's stale test and LOC counts. The pattern worth keeping: in every one of those cases the document's *ideas* held up unchanged, and only its checkable facts rotted — which is an argument for deriving such facts rather than typing them, not for proofreading harder.

## Contact physics build dependency (2026-08-28)

`ipc-toolkit` (implicit IPC: Newton + log-barrier + line-search filter, Li-Kaufman 2020 §5) is now an optional CMake dependency: `SPATIUM_IPC_TOOLKIT` (default OFF), `FetchContent`-declared against ipc-toolkit v1.6.0, matching `SPATIUM_EIGEN`'s opt-in pattern. Linked directly onto consuming targets (`spatium_tests` today) rather than the exported `spatium_sdk` INTERFACE target — unlike `Eigen3::Eigen` (an IMPORTED `find_package` target), `ipc::toolkit` is built locally via FetchContent and isn't part of any install export set, so linking it into `spatium_sdk` breaks `install(EXPORT SpatiumTargets ...)`. Verified end-to-end, not just at configure time: builds ipc-toolkit's own dependency tree (libigl, TBB, spdlog, TinyAD, tight_inclusion, scalable_ccd, xsimd, robin-map, abseil, filib, all via its CPM recipes), links, and a `CollisionMesh` smoke test (`tests/test_ipc_toolkit.cpp`) runs and passes. Default build is unaffected — option OFF by default, `SPATIUM_HAS_IPC_TOOLKIT` compiles to 0, all 811 test cases / 80 800 assertions still pass with the option off.

## Implicit-contact Newton solver, mechanism proven (2026-08-28)

`examples/cloth_sphere_probe.cpp` — a hand-rolled implicit-Euler Newton solve (`E(x) = 1/2(x-x̂)ᵀM(x-x̂) + h²Ψ_elastic(x) + h²B(x)`) combining a quadratic-spring cloth energy with ipc-toolkit's `BarrierPotential`, replacing explicit XPBD's Gauss-Seidel projection for a cloth-on-sphere scene (the sphere triangulated and merged into the same `CollisionMesh` as a static obstacle — ipc-toolkit's own collision API has no entry point for Spatium's analytical `ContactSurface<Sphere<2>>` directly). Three real bugs found and fixed via diagnostic instrumentation, not guessed: Newton starting from the unchecked predictor `x̂` instead of the last known-safe state (could tunnel before any CCD check ran); an Armijo line search that silently accepted a non-decreasing-energy step whenever every backtracking halving failed instead of erroring out; a CCD min-distance floor set to a physically-meaningful gap (`1e-4`, comparable to `dhat`) instead of a numerical safety margin, permanently deadlocking Newton once contact settled near that value. After all three fixes: confirmed converging via the diagnostic trace (`gnorm` 3.2e-2 → 2.6e-6 over ~14 iterations, no explosion, no tunnelling) — the mechanism genuinely works, beating XPBD's 0/20 baseline on the one config actually traced. Not done: a full multi-config sweep matching the XPBD investigation's own 20-point discipline (see Backlog → Contact physics / RSC) — even one reduced-scale config (81 cloth vertices, coarse sphere, 90 frames) took multiple CPU-minutes under ipc-toolkit's TBB-based collision detection (thread-pool overhead dominates at this problem size). Deliberately not chased further at the time — mechanism correctness was the goal, not throughput.

## Code and build conventions, dependency graph (2026-08-28)

A second truing-up pass, prompted by the same "code drifts from doc, doc drifts from code" pattern the 2026-08-26 architecture audit found, this time going into the code itself rather than just the docs describing it: `docs/conventions.md` (new) names the namespace/subdivision/error-handling/`measure()`-alias/doc-comment rules that were previously followed inconsistently per-domain, not written down anywhere. Applied: `algebra/` moved to `inline namespace algebra` (matching every other domain's per-namespace convention while staying bare-`spatium::`-compatible for the rest of the tree's existing unqualified usage — hit and fixed a real GCC constraint along the way, that `inline` must appear on the *first* reachable declaration of a namespace in a translation unit); `physics/atomic/` split out from `physics/`'s flat top level; `geometry/triangle.hpp`/`circle.hpp`'s `measure()`/`area()` alias inversion fixed (the project's own documented rule, violated in 2 of 3 sampled files). Build: `cmake/SpatiumTarget.cmake`'s `spatium_add_example()` replaces 16 copy-pasted `examples/CMakeLists.txt` blocks; `CMakePresets.json` gives 8 named, reproducible presets for what previously existed only as 10 undocumented `build-*/` directories. `scripts/gen_dependency_graph.py` generates `docs/dependency-graph.dot` from real `#include` edges (aggregated to domain level), checked in CI (`.github/workflows/ci.yml`'s `dependency-graph` job) so `architecture.md`'s Module Map can't silently drift from the real dependency graph again — confirmed real cross-domain edges the prior narrative description missed (`spaces`↔`mesh` mutual dependency, `discrete`→`geometry`/`algebra`, `geometry`→`spaces`). All 793 tests green throughout, phase by phase.

## Geodesic procedural generation, render/ engine consolidation (2026-08-28)

- `examples/geodesic_procgen_demo.cpp` — the first two of three manifold-backlog ideas raised this session (the third, hyperbolic ray-marched rendering, is in Backlog → Manifold applications): `make_bumpy_sphere()` (new factory, `spaces/implicit.hpp`) → `marching_cubes()` → farthest-point sampling via repeated `geodesic_voronoi()` calls (each new site is the vertex the last call's own distance field says is farthest, so no separate single-source Dijkstra needed) → `face_labels()` → `BVH<Triangle3>::ray_cast()` against the real mesh, rendered through the same `render::Camera`/`parallel_for_rows`/`supersample_pixel`/`write_png_rgb` pipeline `tumbling_body_demo.cpp` established as the pattern. Zero new geometry algorithm — genuinely pure composition, as scoped.
- `render/color.hpp` (`hsv_to_rgb255`) and `render/sky.hpp` (`Sky`/`make_starfield`/`sample_sky_color`/`sample_sky`/`random_sky_dir`) promoted out of `examples/io_helpers.hpp` into the library proper — the same "found needed by another caller → promote to the engine" path `camera.hpp`/`parallel_for_rows.hpp`/`supersample.hpp` already went through, triggered this time by `geodesic_procgen_demo.cpp` needing the color/star-scatter pieces but not the whole-sky gradient. `examples/io_helpers.hpp` now holds only `confirm_overwrite()` (an example-CLI convention, not a rendering primitive).
- Real bug found and fixed, not just moved: `sample_sky_color()`'s background brightness gradient ("brighter near mid-latitude, dark at the poles") is correct for the wide-FOV whole-sky views the 3 GR raytracers use, but at `geodesic_procgen_demo`'s close, narrow-FOV single-object framing only a small patch of that gradient is ever visible — reading as an unexplained dark smear rather than a sky. Fixed with an additive `Sky::wide_sky` bool, default `true` (unchanged branch taken by every existing caller — spot-checked by re-rendering `blackhole_demo`, one frame of `blackhole_gr_demo`, and one frame of `tumbling_body_demo` after the move, all visually consistent with their established look); `false` gives a flat background at `tint` instead, used by `geodesic_procgen_demo`. New tests in `tests/test_render.cpp` cover both branches plus `hsv_to_rgb255` and per-star findability. 822/822 tests green.

## Riemannian optimization (2026-08-28)

Gradient descent on manifolds (Riemannian SGD): `riemannian_minimize()` in `algebra/calculus.hpp`, Armijo backtracking line search along the manifold. Retracts via each space's own `exp_map()` — no separate retraction abstraction needed. `raise_gradient()` + `project_tangent()` convert the ambient `gradient()` covector into the actual Riemannian gradient first — index-raising matters once the ambient metric isn't Euclidean, e.g. Hyperbolic's Minkowski form (a real bug caught by a failing test, not review: on Hyperbolic the optimizer silently never moved before this fix, because the uncorrected "gradient" pointed the wrong way and line search rejected every step).

## Hyperbolic ray marching (2026-08-29)

Third and smallest of the three manifold-backlog ideas raised 2026-08-28 (the other two, Riemannian optimization and geodesic procgen, are above). `examples/hyperbolic_tessellation_demo.cpp` — sphere-tracing through `Hyperbolic<3>`'s own hyperboloid-model metric, not the mesh+BVH stack (Euclidean ray-primitive intersection, no way to tessellate hyperbolic space into flat triangles) and not RK4 geodesic integration (unlike Schwarzschild/Kerr, `Hyperbolic<N>::exp_map()` is exact, no ODE to integrate). Camera fixed at `Hyperbolic<3>::origin()`; because its tangent space there is genuinely Euclidean (`metric_at` reduces to the ordinary dot product once the normal-orthogonality constraint zeroes the ambient time coordinate), the existing Euclidean `render::Camera`/`camera_ray_dir` pinhole formula gives the initial ray direction directly, embedded into the hyperboloid's `Vec<T,4>` with a leading zero — no new camera-basis math. Markers sit at three shells of hyperbolic distance along a once-subdivided icosahedron's 42 directions (`mesh::subdivide_once`, reused as unit directions, not as a mesh), sphere-traced via `space.distance()` to the nearest marker as the safe step size.

Two real, non-obvious things found while tuning it, not guessed: (1) a hyperbolic ball's apparent angular size follows `sin(alpha) = sinh(radius)/sinh(distance)`, not the Euclidean `radius/distance` — an initial `radius=0.25` at `distance=1.0` filled a third of the frame, corrected to `0.08`; (2) with only 12 base directions spread across the *entire* sphere of view (not clustered toward the camera's forward hemisphere the way objects in ordinary scenes are), a normal ~40-55° FOV only ever catches 1-3 of them by chance — needed a wide (130°) FOV plus the subdivision (12→42 directions) to show a genuinely scattered field in one still frame; the wide-FOV rectilinear projection itself stretches off-center markers into radiating streaks, an honest projection-formula side effect (not a hyperbolic-metric effect) left in because it reads as a striking, correctly-motivated part of the image rather than a defect. Palette also needed golden-ratio-conjugate hue decorrelation (`hue = fmod(i * 0.618034, 1)`) since consecutive icosahedron-vertex indices are often spatially adjacent, which otherwise clustered same-hued neighbors together on screen. Wired into `examples/CMakeLists.txt` and `flake.nix` (`nix run .#hyperbolic-tessellation`), `nix flake check --no-build` passes, 825/825 tests still green (no library code changed, only the new example + existing `subdivide_once`/`icosahedron`/render-engine reuse).

## SO(3)/SE(3) templated on Scalar (2026-08-31)

`SO3`/`SE3` (`algebra/groups/so3.hpp`, `se3.hpp`) moved from hardcoded `double` to `template<Scalar T = double>`, matching every other Scalar-templated type in the tree — `SO3<double>`/`SE3<double>` behave exactly as before (all call sites updated: `physics/mechanics/lgvi.hpp`, `examples/tumbling_body_demo.cpp`, test helpers), and `SO3<Dual<double>>`/`SE3<Dual<double>>` now satisfy `Group`/`LieGroup` (checked via `static_assert`, plus a real differentiation test — rotating a point about Z with a `Dual`-seeded angle recovers the exact closed-form derivative, no hand-derived Jacobian). This is what Sophus/manif already give C++ robotics/SLAM/pose-graph optimization; Spatium didn't until now.

Two real, non-hypothetical things found while templating, not assumed:
- `Dual<T>` gained `acos()`/`tan()` (`algebra/dual.hpp`) — `SO3::log()` needs the former, `SE3::log()` the latter; neither existed before since nothing had exercised `Dual<T>` through a Lie-group log map yet.
- A real latent bug, unmasked by the change itself: `SO3::log()`'s `sin(angle)` call had no `using std::sin;` in scope (unlike `exp()`, which did) and silently worked anyway because it resolved to the global `::sin` leaked in by `<cmath>` — invisible while `so3.hpp` didn't yet include `dual.hpp`. Once it did (needed for the `Dual<double>` `static_assert`), `spatium::algebra::sin<T>` (from `dual.hpp`) became visible via ordinary unqualified lookup inside the same namespace and shadowed the outer `::sin`, breaking the `T=double` case outright until `using std::sin;` was added explicitly. Caught by the compiler on the very first build, not by review.

C++23-modules partitions (`modules/algebra/groups_so3.cppm`, `groups_se3.cppm`) needed `import :dual;` added — the `Dual<double>` static_asserts reference a type those partitions hadn't previously needed to see. 825→827 tests, one pre-existing unrelated failure (`embedded_base_is_current`, a stale CMake-configure artifact from the 2026-08-31 release's git-history squash — `kSpatiumCommitSha` reads "unknown" in the cached `build/`, needs a `cmake` reconfigure to pick up the real HEAD SHA, untouched by this change).

**Follow-up same day: `SO3::exp()`/`log()` rewritten for real differentiability at the origin.** Building an actual rotation-averaging demo (`examples/primitives_demo.cpp --scene rotavg`, below) exposed that `gradient(f, Vec3{0,0,0})` came back a hard `[0,0,0]` on a genuine rotation-averaging objective — confirmed via a standalone probe against the same objective evaluated 0.01 away (a real, nonzero gradient there). Root cause: `SO3::exp()`'s old `if (angle < eps) return identity();` and `SO3::log()`'s old `if (angle < eps) return AlgebraType{};` were both exactly correct in *value* at the origin/identity but silently *v-independent* there — under `Dual<T>` this zeroed every derivative exactly at the single most common optimization starting point, even though `exp()`/`log()` are analytically smooth there (the singularity is entirely in the intermediate `axis = v/angle` division and `acos'(1) = -1/√0`, not in the functions themselves). Fixed by re-deriving both in terms of `θ² = v·v` (a smooth polynomial in `v`, unlike `angle = √θ²`) with Taylor-series coefficients near the origin instead of a constant-value shortcut — `SO3::exp()`'s `R = I + a(θ)K + b(θ)K²` now uses `K = skew(v)` directly (not `skew(v/angle)`) with `a,b` the (now removable-singularity-safe) `sin(θ)/θ`, `(1-cos θ)/θ²`; `SO3::log()`'s general formula was refactored around `vee(R-Rᵀ)` (always linear/smooth in `R`) scaled by a coefficient that's Taylor-expanded near `cos_angle=1` *before* ever calling `acos` there, since `acos`'s own derivative diverges at exactly that point. The `angle≈π` branch is untouched — a genuine, non-removable coordinate singularity of the axis-angle chart itself (two antipodal axes represent the same rotation there), not a Dual artifact, and out of scope for this fix. Two new regression tests (`tests/test_algebra.cpp`): `SO3<Dual<double>>` differentiates correctly *at* `v=0` (checked against the known first-order Taylor expansion `exp(v) ≈ I + skew(v)`) and through a `log(exp(v))` roundtrip at `v=0`. 827→829 tests, same one pre-existing unrelated failure as above.

Re-deriving `SO3::exp()` this carefully surfaced the same *shape* of problem in `SE3::exp()`'s translation part — and cross-checking it against an independent ground truth (a 40-term brute-force Taylor sum of the 4×4 se(3) generator's own matrix exponential, computed independently of `so3.hpp`/`se3.hpp`) turned up a **real, pre-existing, silent correctness bug**, not just a differentiability gap: the old `V` matrix used `skew(ω/angle)` (the unit axis) with coefficients `(1-cos θ)/θ²`, `(θ-sin θ)/θ³` that are only valid for `skew(ω)` (unnormalized) — wrong for any nonzero rotation combined with a translation. Confirmed numerically: the old formula gave `t=[0.897, 0.646, 0.218]` against the ground truth's `[0.975, 0.603, -0.062]` for a representative `(ω,v)` — not close, and the *z* component even has the wrong sign. No prior test caught it: `"SE3 exp/log roundtrip"` is self-consistent under either convention (both sides used the same wrong one), and `"SE3 exp known-good: 90deg Z rotation + translation"` builds `T` via `from_Rt()` directly, bypassing `exp()`'s `V` matrix entirely — a real, honest test-coverage gap, not carelessness caught late. Fixed by extracting a shared `translation_jacobian(ω)` (the corrected, ground-truth-verified `V(ω) = I + b(θ)K + c(θ)K²` with `K = skew(ω)`, Taylor-safe near `ω=0` the same way `SO3::exp()` is) used by both `SE3::exp()` (as `V` itself) and `SE3::log()` — which now inverts it via the library's own general `invert()` (`algebra/linear_solve.hpp`) rather than a second hand-derived `V⁻¹` closed form, deliberately trading a small runtime cost for not risking a *third* manually-derived formula after the first one turned out wrong. One more regression test (`tests/test_algebra.cpp`): `SE3::exp()` checked directly against the brute-force ground truth, closing the exact coverage gap that let the bug through. 829→830 tests, same one pre-existing unrelated failure as above.

**`examples/primitives_demo.cpp --scene rotavg`** — the demo that surfaced all of this: SO(3) rotation averaging (N noisy measurements of one rotation → the estimate that best explains them all, a real SLAM/sensor-fusion building block), solved via `calculus.hpp`'s existing `gradient()`/`minimize()` with zero hand-derived Jacobian, animated one Armijo-backtracking step per frame (the same per-iteration logic `minimize()` uses internally, exposed step-by-step) so the convergence itself is watchable — bright RGB gizmo (live estimate) converging toward a gray reference (ground truth) alongside dim RGB gizmos (the noisy inputs).

## C++23 modules: caught up to a working state (2026-08-28)

The 2026-04-24 migration (see above) had gone stale — this closes the gap rather than deleting the subsystem, since the missing partitions were genuinely mechanical, not a design problem.

- **New algebra partitions**: `dual`, `ode`, `linear_solve`, `calculus` — added to `modules/algebra/`, registered in `algebra.cppm`'s export list and `CMakeLists.txt`'s `PARTITIONS`. `eigen_interop.hpp` stays deliberately header-only (SSE intrinsics from `Eigen/Core` conflict with `:vec_simd`'s BMI — this was already documented, not new).
- **`physics/relativity/*` folded into `spatium.physics`**: `schwarzschild.hpp`, `kerr.hpp`, `geodesic.hpp`, `accretion_disk.hpp` added to `physics.cppm`'s single-TU `#include` block (no cross-includes among the four, any order works; `geodesic.hpp`'s `Dual<T>`/`ode`/`linear_solve` needs are covered by the algebra partitions above). `physics/atomic/`'s own stale include paths (left over from its 2026-08-28 split out of `physics/`'s flat top level) were also fixed in passing — a real, separate bug: `modules/physics.cppm` still referenced the pre-split `physics/orbital.hpp`-style paths, breaking the modules build specifically since the header-tree build never exercises that file. Drive-by cleanup: the "Block A/B/C/D" milestone-naming comments the 2026-08-26 architecture audit was supposed to have stripped everywhere had survived inside this one file; removed.
- **`spatium.render`: new module**, `modules/render.cppm` — single partition covering `camera.hpp`, `color.hpp`, `spectral.hpp`, `parallel_for_rows.hpp`, `supersample.hpp`, `sky.hpp`. `write_image.hpp` deliberately stays header-only: its whole contract is "exactly one translation unit defines `STB_IMAGE_WRITE_IMPLEMENTATION` before including it," which is a link-time property, not a module one — folding it into the module would bake `stb_image_write`'s externally-linked function bodies into `spatium_render_module` itself, and any TU that also defines the macro (every example does, so does `tests/test_render.cpp`) and links against the module in the same binary would hit duplicate-symbol errors. Same category of documented exclusion as `eigen_interop.hpp` and `mesh/dec.hpp`/`mesh/primitives.hpp` — a real technical call, not neglect.
- **Real GCC modules-ts bug found and fixed, not routed around**: linking a consumer of `spatium.render` against `Sky`/`make_starfield()` failed with `undefined reference` to `std::vector<int, allocator<int>>`'s special members (`_Vector_impl`'s destructor, the default constructor, the move constructor) — reproduced from a clean isolation (removing `Sky` usage from the test made the link succeed; adding real `#include <vector>` in both the producing module and the consuming TU did not fix it). Root cause narrowed to `Sky::star_buckets`'s element type specifically: switching `std::vector<std::vector<int>>` to `std::vector<std::vector<std::uint32_t>>` (matching the index-type convention already used elsewhere in the codebase — `BVH::Hit::index`, mesh vertex/face ids — not a workaround invented for this bug) made it link cleanly in both the modules and header-tree builds. A genuine, narrow GCC/libstdc++-modules-ts limitation with plain `int` crossing a module boundary through a nested container inside an inlined function, not an application-level mistake.
- **CI**: `.github/workflows/ci.yml` gained a `modules-build` job (configure/build/test via the existing `modules` CMake preset) — the gap ROADMAP had flagged as the actual reason this was allowed to rot the first time.
- **Verified**: 825/825 tests green on the `modules` preset (`build-modules/`, fresh configure) and 825/825 on the existing `build/` tree (which already had `SPATIUM_USE_MODULES=ON` cached from earlier work). `scripts/gen_dependency_graph.py --check` and `scripts/check_claude_md_layout.py` both still pass (neither script inspects `modules/`, so unaffected either way, checked to be sure).

## SPD(n): log-Euclidean manifold on symmetric positive-definite matrices (2026-09-01)

`spaces/spd.hpp`'s `SPD<N, T>` (`N=2,3`) — the natural sibling of SO(3)/SE(3) rotation averaging: covariance-matrix descriptors in computer vision, diffusion tensors in DTI/medical imaging, spatial-covariance features in BCI/EEG. Chosen over the affine-invariant metric deliberately, as the cheap first step (a fuller/affine-invariant treatment is a real, explicitly deferred follow-up, not dropped).

The defining move of log-Euclidean (Arsigny, Fillard, Pennec, Ayache 2006): matrix log is a global diffeomorphism SPD(n) ↔ Sym(n) (ordinary symmetric matrices, a flat vector space), so a "point" is represented directly as `vech(log(S))` rather than `vech(S)` itself. Under that parametrization every `Manifold`/`RiemannianManifold` operation (`distance`/`exp_map`/`log_map`/`metric_at`) is the ordinary flat Euclidean formula — the curvature of SPD(n) is entirely absorbed into the log/exp change of coordinates, not left for `exp_map`/`log_map` to handle the way `Sphere`/`Hyperbolic` must. `vech()` weights off-diagonal entries by `√2` so the plain dot product on `PointType` equals the Frobenius inner product on the matrix itself.

Sharp, checkable consequence of that flatness: the Fréchet mean of SPD matrices under log-Euclidean has a **closed form** (average the logs, exponentiate back) — `frechet_mean()` uses it directly, no `riemannian_minimize()` call, unlike `Sphere`/`Hyperbolic` where retraction genuinely has curvature to contend with. This is the textbook reason log-Euclidean is called "cheap," not just a performance claim.

Deliberately **not** a `Surface` (no `project()`/`normal()`): under log-Euclidean, SPD(n) is flat and full-dimensional in its own log-space coordinates, with no meaningful embedding as a hypersurface of some larger ambient space the way `Sphere`/`Hyperbolic` have — so there's no normal direction to define, and forcing one on (the way `Euclidean<N>`'s own placeholder `normal()` does, purely to satisfy the concept for `Mesh<Euclidean<N>>` compatibility) would silently break `riemannian_minimize()`'s `project_tangent()` for any caller who did try to use it that way, by removing one real degree of freedom from every gradient step. `Manifold + RiemannianManifold` is the honest concept-hierarchy fit here, not `Surface`.

`N` restricted to 2 and 3: log-Euclidean's real numerical engine is the matrix logarithm/exponential of the (small) SPD matrix, computed via closed-form eigendecomposition — for symmetric matrices, eigenvalues are the real roots of the characteristic polynomial, so `solve_quadratic()`/`solve_cubic()` (already in `algebra/polynomial.hpp`) give them directly with no iterative eigensolver. This reuses existing infrastructure and deliberately stays out of the separately-tracked, much bigger general SVD/eigendecomposition backlog item (Backlog → Native math) — general `N` needs that to land first. 3×3 eigenvectors come from the null-space-via-cross-product of the largest-norm row pair of `S - λI`, with a Gram-Schmidt pass against already-placed columns to stay robust for near-repeated eigenvalues (a real, if secondary, numerical edge case, not one this scope tries to handle exactly).

11 new tests (`tests/test_spd.cpp`): eigendecomposition reconstruction + orthonormality for hand-verifiable 2×2/3×3 cases (a known 45°-rotated eigenbasis, not just the trivial diagonal case), `to_spd(from_spd(S))` roundtrips, and `frechet_mean()` checked against an independently-derivable ground truth (the per-eigenvalue geometric mean, both axis-aligned and under a shared rotated eigenbasis) rather than only self-consistency. 825→836 tests, same one pre-existing unrelated failure noted throughout this file.

## SPD(n): affine-invariant metric (2026-09-01)

`spaces/spd.hpp`'s `SPDAffineInvariant<N, T>` (`N=2,3`) — the fuller treatment explicitly deferred when log-Euclidean landed above ("do it seriously later, not dropped"). Same matrices, different metric: `<U,V>_S = trace(S^-1 U S^-1 V)`, invariant under every congruence `S -> A S A^T`. Log-Euclidean's `<U,V> = trace(UV)` via `vech(log(S))` is *not* invariant under that action — the real mathematical gap this closes, not just an alternate style.

Where `SPD<N,T>` keeps a point as `vech(log(S))` specifically so every `Manifold` operation collapses to flat arithmetic, `SPDAffineInvariant` keeps a point as the literal SPD matrix and lets the curvature show: `exp_map`/`log_map` are the standard closed forms `Exp_S(V) = S^(1/2) exp(S^(-1/2) V S^(-1/2)) S^(1/2)` and `Log_S(Q) = S^(1/2) log(S^(-1/2) Q S^(-1/2)) S^(1/2)`. Both reduce to the same eigendecomposition-based symmetric matrix log/exp `SPD<N,T>` uses (`S^(-1/2) V S^(-1/2)` is symmetric whenever `V` is), so this reuses `detail::apply_eigen_sym` directly — no new numerical machinery, same `N=2,3` restriction for the same reason. `distance()` is defined as `sqrt(metric_at(p, log_map(p,q), log_map(p,q)))` rather than the textbook direct generalized-eigenvalue formula — the general Riemannian identity (distance = norm of the initial geodesic velocity), keeping the three consistent by construction instead of risking independent hand-derivations drifting apart (see the `SE3::exp()` translation-Jacobian bug above for what that risk costs in this codebase). No `from_spd()`/`to_spd()` needed here: `PointType` already *is* the matrix.

SPD(n) under this metric is a Hadamard manifold (complete, non-positively curved, uniquely geodesic) — geodesics genuinely bend away from the SPD cone's boundary (`Log_S(Q)` diverges as `Q` approaches singular) instead of being straight lines in disguise. That extra fidelity costs the closed form: `frechet_mean_affine_invariant()` is a real fixed-point (Karcher-mean) iteration — average `log_map()` to every sample, retract, repeat until the average tangent vector is ~zero — genuinely exercising the exp/log pair for the first time in this file, unlike log-Euclidean's `frechet_mean()`. No ambient projection step needed (unlike `algebra::riemannian_minimize()` on `Sphere`/`Hyperbolic`), since this space's tangent space is already the unconstrained space of symmetric matrices; hand-rolled rather than calling `riemannian_minimize()`, which requires `HasNormal`, deliberately absent here for the same reason `SPD<N,T>` has none.

**A real, pre-existing bug in `detail::eigen_sym` found and fixed along the way** (not introduced by this feature — latent since the log-Euclidean commit, just never exercised): the 3×3 null-space-via-cross-product construction silently produced a **zero eigenvector** whenever an eigenvalue had multiplicity ≥ 2 (e.g. `diag(3,3,10)`, an axisymmetric-inertia-tensor shape) — not just the exact `S = cI` case. `A = S - λI` then has rank ≤ 1, so the cross product of *any* two rows vanishes regardless of which pair gets picked, and worse, a genuinely *repeated* root is only located by `solve_cubic()` to about `sqrt(machine epsilon)` (a double root's location is a square-root-sensitive function of the polynomial's coefficients — standard numerical-analysis fact, not solver imprecision), so the "vanishing" cross product isn't even close to exactly zero — it's small-but-numerically-unstable noise (~2e-7 in the `diag(3,3,10)` case), large enough to sail past a naive near-zero check while still being garbage. First exposed by `SPDAffineInvariant`'s exp/log-at-the-identity test (`S=I`, the fully degenerate `rank(A)=0` case) — chasing it down surfaced the more general `rank(A)=1` case too. Fixed by detecting multiplicity directly on the eigenvalues themselves (compare `values[i]` against every other root, with a tolerance scaled to `sqrt(epsilon<T>())` — matching the solver's actual achievable precision for a repeated root, not machine epsilon) rather than inferring it from the cross product's shakiness, then falling back to a construction that only needs the single largest-norm row of `A` (crossed with a slot-index-varied seed so two occurrences of the same repeated eigenvalue don't land on the same raw candidate before Gram-Schmidt gets a chance to separate them). The 2×2 path had the analogous exact-tie bug (`S = cI` making the "closest diagonal entry" tie-break pick the same eigenvector for both slots) — no solver-precision subtlety there since 2×2 eigenvalues come from `S`'s own entries directly, not a root-finder, so a plain equality check on `a == d` is enough. 3 new regression tests target `eigen_sym` directly (2×2 and 3×3 scalar-multiple-of-identity, 3×3 axisymmetric).

12 new tests (`tests/test_spd.cpp`): the 3 `eigen_sym` regressions above, concept checks, `contains()` (symmetric + positive-eigenvalue membership, unlike log-Euclidean's unconstrained log-space), exp/log roundtrip at the identity cross-checked directly against `SPD<N,T>::matrix_log_sym` (not just internal self-consistency), exp/log roundtrip at a non-identity base point (2×2 and 3×3), `distance()` symmetry, `distance()` against the independent generalized-eigenvalue closed form for a shared eigenbasis, and `frechet_mean_affine_invariant()` checked both against the log-Euclidean closed form (commuting samples) and against its own first-order optimality condition (non-commuting samples, where no closed form exists to check against directly). 842→854 tests, same one pre-existing unrelated failure noted throughout this file.

## Declarative scene DSL: `io::build`, analytic-first, trace-based (2026-09-09)

`spatium::io::build::Trace<T>` — a getting-started demo (a donut, in the tradition of every 3D beginner's first Blender project) built out a real, general primitive along the way, not just example glue. Core decisions, all deliberate:

- **Analytic, not mesh-first.** `spaces/offset.hpp`'s `offset_surface(base, thickness)` composes a *new* `ParametricSurface` from `base(u,v) + thickness * base.normal_at(u,v)` — a function, not a second shape draped/projected onto the first (an earlier attempt built icing as a near-duplicate torus formula just to project it back onto the real one; redundant, and the mistake `offset_surface` exists to make impossible). `spaces/sample.hpp`'s `sample_surface_uniform` places points evenly *by actual surface area* via rejection sampling weighted by the first fundamental form's area element (`ParametricSurface::area_element()`, new) — no mesh, no geodesic-Voronoi graph, for the common case where the target is a `ParametricSurface`. Mesh only gets built at `materialize()`, for display, never as the scene's own representation.
- **Trace, not closures.** `Trace<T>` is a flat `vector<TraceNode<T>>`, each node a tagged operation (`Space`/`Offset`/`Scatter`/`Compose`/`Literal`) with explicit fields, addressed by index — not a tree of `std::function` closures (the first draft's mistake: opaque, un-inspectable, can't be walked or re-interpreted by a second pass later). `Handle<T>` wraps `{trace*, index}`; every factory method (`torus()`, `offset()`, `scatter()`, `cube()`, ...) appends one node and returns a handle. `Literal` is the deliberate escape hatch for shapes with no natural (u,v)→R³ formula (a cube has none) — holds a precomputed mesh directly.
- **One motion slot, not three.** Every node has `.moving(fn(point, t) -> point)` — a constant shift, a `Morphism` pipe, or genuine time-dependence all fit the same signature. The donut's opening act (delete the default cube, Blender-tutorial-style) uses it directly: displacement is `p*(1-smoothstep(t)) + turbulence(p,t)*smoothstep(t)`, a pure function evaluated wherever `materialize()` asks for it — no per-frame simulation state. (First attempt pushed vertices *outward* along their own direction for "explode" — on a mesh with fixed face connectivity that just inflates the solid into a blob, doesn't shatter it; collapsing inward reads as "deleted" correctly with zero per-fragment mesh work.)
- **Lazy, deliberately not reactive.** Nodes are immutable values; building one does no geometry work, only `materialize()` walking the trace does. No dependency-tracked incremental recompute (SolidJS/Vue-signals style) — that engine cost has no consumer yet in this scope. If a reactive layer gets built at all near-term, `gpu/derive_christoffel.py`'s declarative-derivation → auto-codegen pipeline (GPU rendering section below) is the more natural home for it than the scene graph.
- **New library primitive, reused nowhere else yet on purpose:** `algebra/noise.hpp`'s `PerlinNoise` (Ken Perlin's 2002 "improved noise", seeded — no shared mutable state, so each use gets an independent field) replaces the ad-hoc per-particle RNG velocity every hand-rolled particle effect in this tree would otherwise reach for.

**Real bug found and fixed while rendering it, not just tuning:** `examples/donut_demo.cpp`'s `--photo` path mixed two color scales — `render::sample_sky_color()` already returns [0,255] (matching every other offline demo's direct pixel write) but the new recursive raytracer's outer wrapper multiplied *everything* by 255 again, including sky contributions blended in through glossy reflection rays. Glossy icing reflecting the sky compounded the double-scale on top of itself, blowing the whole icing surface out to solid white. Fixed by keeping `trace_ray()` in a single [0,1] convention throughout (converting `sample_sky_color`'s result down once, at the point of use) and doing the one [0,255] conversion exactly once, at the outermost call.

`--photo` is now a real small raytracer, not a flat raycast: `Material<T>` gained a `roughness` field (1 = matte/Lambertian, the old implicit default; 0 = glossy) driving a Blinn-Phong specular lobe plus an actual recursive reflection ray for glossy surfaces (depth-limited to 3) — dough is rough, icing is glossy, and the difference is a real light-transport effect, not two different `base_color`s. Per-vertex smooth normals (area-weighted face-normal average, interpolated at the hit point via the BVH's own barycentric weights) replaced flat per-triangle shading, which was reading as faceted plastic on the demo's coarse-ish tessellation.

**Not done, explicitly deferred:** a live Vulkan window (CPU raytraces, Vulkan only presents the frame) — real new plumbing (`viewer::App` has no sampled-texture descriptor/pipeline path today, only mesh/point-cloud rasterization), too open-ended to land blind against a same-day deadline; console output + `--photo` PNG are the two supported modes for now, see Backlog below. Scattered items still sit slightly proud of the surface rather than lying fully flush (a small tangent-plane fix, not attempted this round).

1001 tests (was 989), 12 new (`tests/test_build_dsl.cpp`) covering `PerlinNoise`, `offset_surface`, `sample_surface_uniform`, and the `Trace` DSL's node-kind recording, materialization, and `.moving()` determinism. One pre-existing unrelated failure (`svd cross-validated against Eigen::JacobiSVD` — confirmed via `git stash` to already fail on a clean `main`, an Eigen `JacobiSVD` assertion about thin U/V needing a dynamic-column matrix, nothing to do with this work).

Docs: [`docs/getting-started-dsl.md`](getting-started-dsl.md) (new, zero-barrier-to-entry walkthrough), `docs/api-reference.md` and `CLAUDE.md` updated for every new primitive above — and, caught in the same pass, both were already stale against several *prior* PRs (`terminal_donut_demo`, `scene_demo`, `sound_synthesis_demo`, `spd_relationship_demo`, `cloth_sphere_probe`, and all of `discrete/combinatorics.hpp` were shipped but undocumented) — fixed together, not just today's additions.

## SPD(n) outreach demo: relationship-matrix training, naive vs. geodesic (2026-09-01)

`examples/spd_relationship_demo.cpp` — first concrete piece of a planned outreach/contributor-recruitment demo (see project vision doc), scoped deliberately small: not the full game/event-generation idea, just the one structural claim underneath it, isolated and made legible. Same argument as "Mario Plays on a Manifold" (arXiv:2206.00106), which uses Riemannian geometry in a VAE latent space to guarantee generated game content stays valid against a Euclidean baseline that visibly breaks — here the "content" is a relationship matrix between two entities instead of a game level.

Setup: a 2×2 relationship matrix starts at the identity (uncorrelated). A sequence of similar events keeps nudging it in a fixed "strengthening" direction (an off-diagonal-only symmetric perturbation) — modeling repeated causally-similar shocks, not independent random rolls, applied two ways in parallel: naive plain matrix addition, and `SPDAffineInvariant::exp_map`'s geodesic retraction. Naive addition breaks positive-definiteness at step 13 (nothing about plain addition prevents an eigenvalue from crossing zero); the geodesic version is still comfortably valid at the same step. Enforced as a real regression (`tests/test_spd.cpp`, "Naive matrix addition can break SPD validity; SPDAffineInvariant::exp_map cannot"), not just a printed narrative — 14 steps, matching what's verified numerically healthy. Also emits an SVG line chart of `lambda_min` vs. step for both trajectories, meant to become the visual for whatever writeup eventually surfaces this to the geometry-processing/geometric-ML audience identified as the actual target (see project vision doc) — not r/cpp/HN as the primary channel.

**Real numerical-precision finding, not swept under the rug**: pushed past ~18 steps with the same fixed raw nudge, the geodesic trajectory's `lambda_min` asymptotically approaches the SPD cone's boundary and ordinary double-precision underflow there eventually returns NaN from the eigendecomposition. This is a genuine floating-point-precision phenomenon near a manifold's boundary, not a counterexample to the structural guarantee (in exact arithmetic a geodesic never reaches the boundary of a Hadamard manifold at all) — but it's real and worth remembering if this demo is ever pushed further: repeatedly applying the *same raw ambient tangent vector* at each new point isn't the same as following one smooth geodesic ray, and doesn't inherit that ray's asymptotic slowdown. Both the demo and the regression test stay within the verified-healthy range (14 steps) rather than papering over this.

Deliberately not yet built: the game/event-generation layer itself (player actions, entity hierarchy, causally-coherent generation) — that's the next, separate piece, once this one lands. "По немногу" (little by little), not one large PR.

---

# Backlog

Grouped by topic, not by version — an item sits here until it's ready to become real work, then moves to Completed above with the date it landed.

Each item carries a tag for *intent*, not difficulty or sequencing — Spatium's goal is to become the de-facto standard for computational work on arbitrary mathematical spaces, and that happens by earning the interest of people who know this territory well, not by chasing a release checklist:

- **[course]** — genuinely useful, actively where we're steering the project right now.
- **[want]** — we'd just like this to exist, or it'd be fun/interesting to build; no strong pull yet, no promise of order.
- **[dare]** — big, uncertain, or hard enough that it won't get built by us any time soon — a standing, genuine invitation for whoever wants to take a swing at it.

A PR against any tag here is welcome. So is a PR against nothing here — if you build something we didn't list, that's a real step toward Spatium becoming a standard in its own right, not a detour from this list.

## What is actually live — the index

The rest of this Backlog is long and stays long: the reasoning is worth more than the list, and an item's entry is usually an argument rather than a ticket. But an argument is not an entry point, and until 2026-09-17 there was none — finding "what is being worked on, what is half-done, what was dropped and why" meant reading a thousand lines, and the answer was split across this file, a task list and a plan file, two of which are not in the repository at all.

So: the items we are actually steering by, with their state. Everything below this table is still real; it is simply not what anyone is doing right now.

| Item | State | What would unblock it | Where |
|---|---|---|---|
| Renderer consuming `cook()` | **in flight** | — | Declarative scene DSL |
| Per-instance parameters (`MotionEnv::origin`) | **next** | the row above | Declarative scene DSL |
| Scatter's fixed axis binding (local z to the normal) | **next** | nothing; a fix, and one of two preconditions for instancing scattered items | Object model as manifold substrate |
| Structural motions in the demo (`grow_scale` and friends) | **next** | nothing; the other precondition, and a matter of spelling | Object model as manifold substrate |
| A refusal that explains itself (and `is_placement()`'s name) | **next** | nothing; rewriting the demo fixes one user, not the trap | Object model as manifold substrate |
| Scatter's arbitrary in-plane directions | parked | a consumer — "follows the flow" is its own design, not this fix | Object model as manifold substrate |
| One object deforming another | parked | an index over the `(u,v)` domain — and a decision about the cycle it introduces, below | Object model as manifold substrate |
| `ball_pit_demo` cleanup | parked | nothing; it is just work | Object model as manifold substrate |
| Offset self-intersection, in the library | parked | nothing; a thickness check against the base's minimum radius of curvature | Analytical rendering |
| Full `EdgeRule` (`RoundCap`, `FlatCap`, `ExtendTo`) | parked | nothing | Declarative scene DSL |
| C ABI for the DSL (CUDA / WASM / JS) | next | **the field work, not RSC** — a scene exports only as far as its fields are structural; the donut measures 32 opaque leaves across 107 fields, and whether they fall into a closed vocabulary is the open question | GPU rendering |
| Rewrite `docs/gpu-abi-design.md` | **done 2026-09-18** | — | GPU rendering |
| RSC as search | parked | the ABI, by our own ordering | RSC as search |
| Ball tree on a manifold | parked | nothing. Note the dependency runs *backwards* from how it reads: the tree is what gives RSC a second method to choose between | Interop / ecosystem |
| Vulkan live display | parked | traversal getting cheaper — a four-second frame makes a live window a slideshow | Declarative scene DSL |
| Sound from geometry | parked | sparse linear algebra, of which there is not a line in the tree | Sound synthesis |
| Time as a property of the space | parked | a scene where the difference is visible | Declarative scene DSL |
| A rainbow, then interference | parked | the sound work; they share the machinery | Analytical rendering |

### Two rules this table exists to enforce

**An item shipped in half says which half in its first line.** Not its third paragraph. This is not style: "no renderer consumes any of it" sat stale for a full day inside an item marked done, and was then quoted back as the current state by the person who wrote it.

**An abandoned item is struck through *with the reason*.** The convention already exists here — see the scattered-items entry under Object model — and simply was not applied everywhere. An item that vanishes silently is indistinguishable from one nobody got to.

Both rules exist because a claim decays in three different places with three different mechanisms, and none of them is visible from the other two: a comment citing a line number rots when the line moves, a ROADMAP sentence rots when the code catches up with it, and a task in an external tracker rots by being invisible. Three scripts in CI catch the first kind. Nothing catches the second, which is what these rules are for.

## Mesh processing

- **[course]** Mesh simplification (edge collapse + QEM)
- **[course]** Isotropic remeshing on Surface
- **[want]** PLY import/export

## Analytical rendering

- **[dare]** Fragment shader raymarcher: GPU-native analytical render (GLSL quadric math)
- **[want]** Dual Quaternion
- **[want]** A real rainbow, and then interference — light given the same wave treatment the sound work gives sound. Proposed 2026-09-17, and the two halves are deliberately one entry because the second is what makes the first worth building. A rainbow alone is dispersion: one more `ray_hit` on a sphere with a wavelength-dependent index, two internal reflections, the 42° caustic falling out rather than being drawn. That is a demo. The step is the second half — treating light as the wave the `wave_string.hpp`/`wave_membrane.hpp` steppers already treat sound as, so that a two-slit pattern is *integrated* rather than plotted from the closed-form fringe spacing. Phase becomes a quantity a ray carries, and the same machinery that gives a membrane its modes gives a field its fringes. The honest caveat up front: this is classical wave optics, and the distance from there to quantum mechanics is a real one — amplitudes that interfere are not yet probability amplitudes, and saying otherwise would be exactly the overclaim this file keeps catching. What it genuinely is: the point where the library's wave equation stops belonging to acoustics and becomes a shared piece.

## Sound synthesis

- **[want]** Physically-based sound synthesis — resonant modes of a surface via Laplace-Beltrami eigenvalues, geometric room acoustics via ray-traced reflections on the existing BVH, and time-domain wave-equation integration (string/membrane/percussion) reusing `algebra/ode.hpp`.

  **Shipped:** the wave-equation half. `physics/mechanics/wave_string.hpp` and `wave_membrane.hpp` are explicit central-difference steppers checked against closed-form fundamental frequencies, `io/wav.hpp` is the zero-dependency PCM writer, and `examples/sound_synthesis_demo.cpp` renders with both. Audio output is done, not pending.

  **Blocked:** the modal half, and on a bigger thing than this entry used to claim. Two sentences here were false and are worth recording as such rather than quietly deleting, because both pointed the work at the wrong place.

  The first said modes were *"already computable through the existing heat-method infrastructure"*. Measured, they are not. `HeatSolver` factorises `A + t·L`, which is `t·(L − σ·A)` at `σ = −1/t` with `t = h²` — so the shift is pinned to the mesh spacing rather than to the spectrum. On a 162-vertex sphere that puts `σ ≈ −11.2` against a first non-trivial eigenvalue of `λ₁ ≈ 2.0`, and shift-invert's discrimination `μᵢ = 1/(λᵢ − σ)` comes out at 2.0× across the first ten modes where a spectrum-chosen shift gives 568×. It also gets worse under refinement, not better: `h → 0` drives `σ → −∞` while the physical frequencies stay put. The second cached factorisation is no better — `L_pinned` has no mass term at all and is clamped at vertex 0, which splits the exact threefold `ℓ=1` degeneracy into `(0.152192, 0.152192, 0.193486)`, so it would distort mode *shapes* and not merely shift frequencies.

  The second said *"the only missing piece is audio output"*. The real missing piece is a **sparse generalised eigensolver**, which is a dependency decision rather than an afternoon: Eigen's own eigensolvers are all dense, and its sparse one is a wrapper that links ARPACK, a Fortran library not in the tree. The three routes are ARPACK (a first non-header mandatory dependency), vendoring Spectra header-only under `include/spatium/vendor/` the way `stb_image_write.h` already is, or hand-rolling shift-invert Lanczos over the `SimplicialLDLT` the project already uses — the last being closest to how `json.hpp`, `wav.hpp` and the Jacobi eigendecomposition were done.

  `build_laplacian()` and `build_mass_matrix()` in `mesh/differential.hpp` return `Eigen::SparseMatrix` and need no conversion; what they need is a fresh factorisation of `K − σM` at a small negative shift, since `K` is exactly singular.

## Calculus

- **[course]** Region-aware integration over Spatium's own geometry (Box/Sphere/Polygon volumes, intersection overlap) — `integrate()` in `algebra/calculus.hpp` is deliberately scoped to 1-D; this is the next step the function's own docstring already flags.

## Manifold applications

- **[dare]** Fiber bundles (tangent/cotangent)
- **[dare]** Geodesic FEM (Laplace-Beltrami, heat equation)
- **[want]** Stochastic processes on Riemannian manifolds — Brownian motion / SDEs via Euler-Maruyama through `exp_map` retraction (drift + diffusion tangent vector, the same shape the existing Lie-group integrators already use). Reuses the `RiemannianManifold` concept directly, no new abstraction needed. Genuinely underserved outside bespoke research code — general manifold-ML libraries don't cover this.
- **[want]** Fisher-Rao information geometry as a `RiemannianManifold` — the manifold of probability distributions, metric = Fisher information, distance relates to KL divergence. Real ML utility on its own (natural gradient descent, distributional RL) for exactly the audience this project already names. Brainstormed 2026-09-09 as "a raytracer inside an ML model's mind" (light traveling through KL-divergence space) — the *space* is the genuine value; the raytracing-visualization framing is optional decoration on top and risks looking hollow if built before the underlying math means something concrete. Land the manifold first, decide separately whether a visualization earns its keep.
- Considered and explicitly parked, 2026-09-09: acoustic raytracing through the Ellis/Morris-Thorne wormhole (send a sound pulse down the throat, solve the wave equation against the curved metric, get a real .wav of "what curved spacetime does to a voice"). Spectacle, not utility — doesn't serve any of the audiences this project actually names, unlike the two items above. The underlying general capability (wave-PDE solving against a curved/non-Euclidean metric, not just the flat `wave_string`/`wave_membrane` domains) might be worth having on its own merits someday; the wormhole framing specifically isn't a roadmap item, just a joke that came up.

## Object model as manifold substrate

Unlocked by the object→exact-`Surface` bridge (extending `geometry/surface_adapter.hpp`'s `ShapeSurface`, in flight as `feat/scene-object-model`): any renderable scene object exposing a real Surface, not an approximate tangent-plane one, means every manifold operation in the library (geodesics, Voronoi, parallel transport, DEC, contact, Riemannian optimization) becomes available on it for free instead of needing per-shape hand integration. Brainstormed 2026-09-06, honestly triaged (not everything below is equally close):

- **[course]** Recursive procedural world generation — geodesic Voronoi/Delaunay cell on a surface, one object per cell oriented by the normal, recurse per-object for a fractal scene. Mostly composition of what already exists (`mesh/voronoi.hpp`'s geodesic Voronoi, `Point`/`Morphism` pipes) once the object bridge lands; not new math.
- **[course]** Manifold-native RL as a new RSC domain — agent state = point on a `Surface`, action space = tangent space, transitions via `exp_map`. Fits RSC's existing REINFORCE-dispatcher methodology directly; also the concrete way to eventually fill the "geometric/manifold-constrained control theory" and "real-time control" candidates raised in the 2026-08-31 industry-directions session.
- **[course]** `GeodesicPath<S>` as a first-class 1D `MetricSpace` — a geodesic curve carries its own arc-length parametrization, so points/distance/Voronoi work along a route. Real routing use case: waypoints along a path, optimal stop placement.
- **[want]** Geodesic-guided growth simulation (cracks/roots/vessels) via existing parallel transport + distance fields — doesn't need the object bridge at all, already buildable on `Sphere`/`Hyperbolic`/`ParametricSurface` directly. Gallery-tier demo, not new library code.
- **[want]** Manifold-native game movement (walk on an asteroid, "up" = surface normal, Mario-Galaxy-style) — same status as above, a composition/demo of existing primitives (`Surface`, `exp_map`, `normal`, `render::Camera`), not a new subsystem. Sharpest concrete version, 2026-09-09: a `Hyperbolic<3>` dungeon crawler — exponential volume growth with distance from center makes rooms feel unbounded and corridors branch fractally, and since collision/physics/traversal already go through `exp_map()`/`log_map()`, enemies walking geodesics and light bending through the Poincaré-disk geometry falls out for free rather than needing hand-written non-Euclidean-movement logic. Worth building specifically because it's the cheapest of the "moonshot" ideas raised this session and it substantiates a claim the HN launch copy already makes ("build a game where up is whatever surface you're standing on") rather than leaving it as an unproven line.
- **[want]** Multiscale geometry pipeline (scan → mesh → `Quadric`/parametric fit → geodesic segmentation → recurse per segment) — not a separate idea, a continuation of the already-noted point-cloud → mesh reconstruction industry candidate (2026-08-31 session): reuses `ImplicitSurface`/`marching_cubes`/geodesic Voronoi, the new piece is only the fit-then-recurse loop.
- **[dare]** Cloth/shell dynamics intrinsic to a curved surface — solving the material's own equations of motion in tangent space accounting for curvature, distinct from the already-shipped "drape cloth onto a fixed obstacle" work (`examples/cloth_sphere_probe.cpp`, Contact physics section above). Real research-level lift, not a quick extension of XPBD.
- **[dare]** Topology optimization on an evolving surface (solve a stress/heat PDE on the current surface, remove material, iterate on the new one) — needs real elasticity FEM, well beyond the current DEC/heat-method infrastructure (which only covers scalar diffusion, not stress). Exciting, not scoped.
- Explicitly NOT a new subsystem, just a consequence of the type system once the object bridge lands: composing nested local frames (a point on a robot → a point on one of its sensors → a point on that sensor's own chip surface) through the same `Point<Space>` vocabulary at every level. Worth a demo to show it falls out for free; not separate work.
- **[want]** Reconstruction pipeline for real-world objects — photos/video → geometry → `ImplicitSurface` → existing Spatium post-processing (smoothing, simplification, geodesics, self-intersection checks, OBJ export). Corrected 2026-09-15: this used to end "Not researched yet" and carry `[dare]`. It has been researched, and the research changed the shape of the item. The learned-latents-to-SDF-decoder route (a JEPA-style encoder feeding a small SDF decoder) was tried and abandoned — the decoder produced a blob, not a surface. A frozen VGGT-1B was then validated locally with no training at all and produced a recognizable mesh from a DTU scan. **None of that work is in this repository** — no directory, no reference, nothing to build on — so what is actually open here is small and specific rather than a research programme: a load path in `io/`, a conversion utility, and somewhere for it to live. The reconstruction itself is a solved external step.

## Declarative scene DSL

### Direction, settled 2026-09-14

Where this is going, so the individual items below read as one line of work rather than a pile of improvements.

**The authoring layer and the execution layer are different layers, and only the first one may be open.** Spatium's whole argument is that a space is an extension point, and the same should hold for scene operations — a user with a better geodesic solver or a new kind of node should not have to edit our headers. But openness is bought with indirection, and indirection is exactly what cannot be lowered to a GPU kernel or fused across operations. The resolution is not to pick one: type erasure, registries and anything else ergonomic belong to *describing* a scene, which happens once; *evaluating* it must be monomorphic, either at compile time or compiled from an IR. This is why the trace was built as a flat record of tagged operations rather than a tree of closures in the first place — closures cannot be lowered, data can. The trace is the IR.

This is not a novel architecture and should not be designed as though it were. Dr.Jit (Mitsuba 3's backend) traces C++ operations into an IR and JIT-compiles to CUDA/OptiX/LLVM for exactly this class of problem; Taichi captures an AST and compiles to CPU/GPU kernels; Halide separates the algorithm from its schedule; Kokkos gets one source onto both CPU and GPU through lambdas and execution spaces, and forbids `std::function` inside kernels for the reason above. What would be ours is the particular IR over Spatium's own operations, not the shape of the solution.

**Measured, not assumed** (`benchmarks/bench_trace.cpp`, on the shape `examples/donut_demo.cpp` actually builds: 19 800 single-mesh nodes, 12 vertices each, one motion hook per node). Three plausible performance hypotheses were tested and all three came back smaller than expected, which is the reason this section exists in this form:

| suspected cost | measured |
|---|---|
| the per-vertex indirect call | ~3.3 ns/vertex — about 4% of a realistic motion callee, invisible unless the callee is trivial |
| 9.7 MB of duplicated `PerlinNoise` tables captured into 19 800 closures | ~3% of frame time |
| materializing the nodes with no motion hook at all | **~20% of frame time** |

So the dominant structural cost is none of the things that looked like the problem: it is that the same 12-vertex mesh is stored and copied per instance, because the trace has no notion of instancing. `std::move_only_function` was still worth adopting, but for what it enables (hooks owning move-only or shared state, const-correct invocation through the `const Trace&` that `materialize()` holds) rather than for speed, and the commit says so.

The order that follows from this: instancing and the structural field representation together, then the full set of edge rules that fields make expressible, then opening the node set, then time. Revised 2026-09-15 — the field work moved ahead of the edge rules, because `EdgeRule`'s remaining values carry data (a cap radius; a reference to another node for "extend to this surface"), so they cannot be enum values and the rule set cannot be completed until fields exist. That argument is worth more than the original one: caching and GPU lowering will be true later, whereas a value with nowhere to put it is true now.

### Chart versus manifold, recorded 2026-09-15

Splitting `offset()` from `offset_shell()` needed a closedness test, and writing it surfaced something larger than the function.

`periodic_u()` / `periodic_v()` describe the **chart**, not the surface. A sphere's chart is the rectangle [0,2π]×[0,π], and as a subset of R² that rectangle has an edge — which is what `periodic_v() == false` honestly reports. But the chart's `v` edges map to the poles: two *points*, not two curves. The surface has no boundary; the chart does. So a closedness test cannot read the flags, because the flags answer a question about the parametrization while the question being asked is topological. `is_closed()` therefore looks at the geometry — a direction closes by being periodic, or by both its edge curves collapsing to a point.

The general form of this, which is what `SurfaceWithBoundary` will actually need, is not "does the chart have an edge" but **what the chart's edge maps to**, and there are three answers:

| the chart's edge maps to | meaning | example |
|---|---|---|
| nothing (the map is periodic) | closed, seam only | torus, in both directions |
| a point | closed, pole — the surface passes through, it does not stop | sphere, in `v` |
| a curve | a genuine boundary | a band cut from a torus; a tube's open ends |

Only the third case needs a rule at all, and `EdgeRule::ZeroThickness` is one answer to it, not the general one. The first two need no rule and must not be made to ask for one — which is exactly the trap the pre-split `offset()` fell into by having a single operation cover all three.

`is_closed()` today collapses this to a yes/no because that is all `offset()` needs. The three-way classifier is what the boundary concept should be built on, and it should be built before the concept, not after.

### Render levels, settled 2026-09-15 — and what the number actually says

A node now carries a `RenderLevel`: `Exact` (closed form), `Tessellated`
(triangles in a BVH), `Newton` (Newton's method on the `(u,v)` map).
Inferred in one place — an exact form when the node has one, tessellation
otherwise, and **Newton never inferred**, because at three orders of
magnitude nobody should pay it by accident. `.rendered_as()` overrides
and refuses on the spot a level the node cannot serve.

**The argument for it is build cost, not per-ray cost.** This was got
wrong twice on the way here and is worth stating flatly, because the
per-ray numbers point in whichever direction the leaf ratio happens to
push them:

| scene | exact leaf vs triangles | per ray | build |
|---|---|---|---|
| 19 800 specks | quadric vs 396k tris | exact **20% slower** | 11.8 ms vs 291 ms |
| one dough | torus vs 25 600 tris | exact **2× faster** | ~0 vs 10.8 ms |
| 64 tori, through the DSL | torus vs 147 456 tris | **a wash** (44 vs 49 ns) | 0.011 ms vs 78.8 ms |

The dough's 2× was read at the time as the speck result reversing. It
was not: one leaf against 25 600 is an extreme ratio, and adding 63 more
tori brings the per-ray cost back level. Nor is it true, as a related
guess had it, that the exact path's advantage grows with scene size in
*traversal* — the third row is the test of that and it came back flat.

What survives all three is that **the tessellation never exists**. Seven
thousand times the build cost in the DSL row, ~8× the memory in the
speck row. An animated scene pays the build on every frame while the
per-ray cost is paid once per pixel, and that is the whole case.

A finding from the demo, kept because it is the mechanism working:
`donut_demo` reports `0 exact, 19 809 tessellated, 0 newton`. The dough
is an `Offset` carrying a noise bump, so it stopped being a torus the
moment it got bread texture. An exact form is a promise about the shape,
and a bumped torus cannot keep it.

### `Chart`, settled 2026-09-15 — the requirement the operations already had

The DSL was open to new *operations* and closed to new *spaces*, and the
reason turned out not to be the one written down. It was not that `Kind`
is a closed enum, and not that the node type is fixed. It was that the
requirement every DSL operation shares had no name, and travelled through
the signatures disguised as the concrete type `ParametricSurface<T>`.

The claim is checkable, which is why it settled the design. Reading the
three operations for what they actually call:

| operation | what it calls |
|---|---|
| `offset_surface()` | `evaluate`, `normal_at`, `domain`, `periodic_u/v` |
| `sample_surface_uniform()` | `area_element`, `domain` |
| `parametric_mesh()` | `evaluate`, `domain` |

Not one of them calls `project()`, `normal()` or `exp_map()` — so not one
of them uses the `Surface` concept. Generalizing `resolve_surface` to
`Surface`, which is how this item was written for a week, would have
handed the operations methods nobody calls while losing every method they
need. The right name is `Chart`: a (u,v) parametrization, which is a
choice of coordinates *on* a space, not a kind of space. It is
deliberately not part of the `Set → Manifold → Surface` hierarchy — a
space can carry many charts or none.

`chart_of(space)` is the extension point, found by ADL, the same shape as
`point_to(p, surf)` in the contact concepts and for the same reason:
compile-time overload resolution, no runtime registry, no indirect call.
`ParametricSurface<T>` gets the identity overload and remains the erased
form every chart converts into, because a trace holds nodes of one type
and cannot be templated per node.

`Sphere<2, T>` gets its **first** chart here, not a generalization of an
existing one — as a `RiemannianManifold` it has no `evaluate(u, v)`
anywhere, which is the real reason it could never enter the DSL. Note
which sphere: `Sphere<N, T>` is the N-sphere in R^{N+1}, so the ordinary
sphere in R³ is `Sphere<2, T>`. `Sphere<3, T>` lives in R⁴ and gets no
overload, which is not an omission — there is no chart of that shape to
write, and `Chartable<Sphere<3, double>, double>` is `false` on purpose.

The chart's `periodic_v` is `false` and its `v` edges are the two poles:
the middle row of the chart-edge table above, and the reason `is_closed()`
reads geometry instead of these flags. The poles are also where the chart
degenerates — `area_element` is `r² sin v`, zero there, and the
finite-difference `normal_at` with it. That is a property of the chart,
not a defect of the sphere.

A sphere node deliberately leaves the `exact` slot empty. That slot means
"a renderer can hit this without triangles", nothing in the library can do
that for a `Sphere` today, and filling it would make `render_level()`
claim `Exact` for a shape no renderer can hit exactly — a picture right
about the shape and wrong about the scene. A `sphere()` factory recording
`BoundedQuadric` is available later and is a different change.

**`Kind` stays closed, and adding a space does not go near it.** Its five
values — Space, Offset, Scatter, Compose, Literal — are *operations*, not
shapes; a sphere is a new chart on the existing `Space` node.
`Kind::Literal` already covers shapes with no (u,v) map at all, so the
shape set was open through the chart or through `Literal` before any of
this. The single thing that would open `Kind` is **CSG boolean on
analytic surfaces** — `union(torus, sphere)` as an operation rather than a
mesh merge, whose result is a new analytic surface not expressible as any
of the five. It is not planned: `geometry/boolean.hpp` does this for
polygons, not surfaces, and there is no consumer. Recorded here so the
question is not reopened from scratch in a month.

### Open items

- **[course]** Classify what a chart's edge maps to — point, curve, or nothing — per direction, generalizing `is_closed()`'s yes/no. See "Chart versus manifold" above. This is the thing `SurfaceWithBoundary` (`boundary_distance`, `is_on_boundary`, what `exp_map` does past the edge) should be designed on top of, and it is cheap enough to have before there is a second consumer for the concept itself.
- ~~**[course]** A Debug job in CI.~~ — **done**, PR #34, 2026-09-15, and the entry stayed listed as open for the rest of that day. The reason it was worth doing, kept because it is the argument and not just the outcome: every job that ran tests built Release, so `NDEBUG` removed `assert` and a whole class of defect — misusing another library in a way that library catches with an assertion — passed green. Not hypothetical; it is how the `JacobiSVD` thin-U/V misuse survived, and how "the test passes now" got written down as "the old failure note went stale". Struck through rather than deleted, same as the `scatter()` orientation entry, because an item that outlives its own completion is the failure worth seeing.
- **[want]** `EdgeRule::ExtendTo(node)` as a structural node holding a trace index, exactly like `Offset::base` already does — independent of the field work, since a reference to another node is topology rather than a per-instance parameter. `RoundCap(radius)` is the one that genuinely needs a parameter and should not drag `ExtendTo` along with it.

- **Two boundaries on the scene work, written down 2026-09-16 because both were about to be treated as "while we are in there".**

  **Light, shadows and textures are edits; time is not.** A light is a field on `Material`, a shadow is a second ray in a recursion that already exists, a texture is a field plus a sampling call — three changes inside files that are already there. **Time is not in that list.** Animation already exists through `.moving()` and is enough for a cube that crumbles; "time as a property of space" is the series item about the metric — gravitational dilation, an object's own clock — and it needs a scene where the difference is *visible*, or there is nothing to check it against. Folded into a render pass it would become another multiplier on `t`, which is precisely the wrong shape.

  **Sound from geometry is gated on sparse linear algebra, not on the scene.** `io/wav.hpp` writes and `wave_string.hpp`/`wave_membrane.hpp` integrate, so the ends exist — but resonant modes of a surface are the generalised eigenproblem `L φ = λ M φ` on sparse matrices, and sparse is still `[course]` above with no native default at all. There is a narrower route that does not need a new sparse eigensolver: shift-invert Lanczos over the cotangent Laplacian and mass matrix in `mesh/differential.hpp`, reusing the `SimplicialLDLT` factorisation `HeatSolver` already performs. That is still its own numerical task, and it is `SPATIUM_EIGEN=ON`-gated, so a default build would not have it. Until then the honest options are a recorded sample — which defeats the point — or the wave equation on a 1D/2D domain, which is real physics but is not the surface, and calling it "the cube's sound" would be the same overclaim this file keeps catching.

- **[course]** **An IR evaluable in `Dual<T>` gives exact normals in `normal_at` and `offset_surface`.** A goal in its own right, not a by-product of structural fields, and recorded separately so it is not lost when someone later asks what fields were for.

  `ParametricSurface::normal_at` computes its normal by **finite differences** (`spaces/parametric.hpp:150-155`, step `1e-6 × extent`), and `offset_surface` is built directly on that normal — so every offset surface in the library, the donut's icing included, carries a differencing error that nothing currently removes. `Dual<T>` satisfies `Scalar`, so a field tree generic over its scalar evaluates in `Dual<T>` and yields an exact derivative with **no differentiation rule written per node**: the derivative comes from the arithmetic, not from a second implementation of the tree.

  What makes this worth its own line: it is the library used on itself. The IR exists to make scenes lowerable, and the same IR makes the library *more accurate* than it is without one. That is a different kind of argument from "20% of a frame", and a stronger one — an accuracy result does not depend on the scene's size.

- ~~**[course]** Instancing in the trace~~ and ~~**[course]** Structural fields~~ — **the library half landed 2026-09-16**, merged as #42 (fields), #45 (`cook()`), #44 (the instance leaf). What that means concretely, since "instancing" turned out to name three separable things:

  `Field`/`VecField` are flat pools of tagged ops with an opaque callable as a *leaf*; `cook()` turns a `Trace` into objects with deduplicated shapes; `Instanced<S>` is a reference to a shape plus a placement, so a BVH leaf is one test rather than N. Measured: 2 shapes for 19 801 objects and 138× fewer vertices held; 16.1× BVH build, 21.6× memory, 1.08× traversal.

  **What decides instanceability is not opacity but whether a motion reads the point.** `A(t) + p·s(t)` with both terms opaque is still a placement, because cost is per vertex and those terms are per object — so `PerlinNoise` never had to become an expression. The donut now reports `35209 placements, 21 deformations`.

  ~~`19804 placements, 9 deformations`, with the nine being the exploding cube's fragments, which genuinely deform.~~ — **both halves wrong, corrected 2026-09-17.** The counts rotted, which is routine. The explanation was wrong when it was written, which is not: there are no cube fragments. The cube is a single node whose motion is `p * (time < 0.12 ? 1 : 0)` — a uniform scale to nothing — and the explosion is the *dust*, which is separate `flake()` nodes with their own motions. Nothing about the cube deforms.

  What the 21 actually are, printed by the demo rather than asserted here, because a claim about this has already been wrong once:

  ```
  refused nodes: 1xLiteral 2xOffset 18xScatter
  ```

  The `Literal` is the cube, the two `Offset`s are the dough and the icing, the eighteen `Scatter`s are the sprinkle groups. The table is a `Literal` too and is *not* refused, because its motion is written structurally. Every refusal is over spelling, not over deformation.

  That also reconstructs the historical nine exactly: six `Scatter` (one per colour, before the bands were split three ways) plus two `Offset` plus one `Literal`. **The number was right and the explanation was wrong from the day it was written** — there were never any cube fragments; the cube is one node and the explosion is the dust.

  A number going stale is caught by re-running. A *reason* going stale is caught by nothing.

  ~~**Still open, and it is the part that moves a number:** no renderer consumes any of it.~~ — **wired 2026-09-16.** `donut_demo` builds a `BVH<Instanced<BoundedQuadric>>` beside its triangle tree and takes the nearer hit; `.moving()` stopped dropping the exact form for a *placement* (a deformation still drops it), which is what made the dust eligible at all. Measured on the same scene with only the render path toggled, at `t=1.5` where there is dust to see: **1.95 s instanced against 2.82 s tessellated, 1.44×, 31% of the frame.** The demo's own line now reads `render levels: 19800 exact, 9 tessellated`, where it had read `0 exact, 19809 tessellated` since it was written.

  Two measurement mistakes are recorded with it, because both were made here and both are easy to repeat. The first: `t=3.9` is a frame in which the dust has already dissolved, so a "56% of the frame" taken there measured traversal of geometry that is not in the picture. The second: comparing 2.02 s against 2.47 s changed the speck's *shape* and the render path at once, which measures neither.

  **Rotation in a placement, 2026-09-17.** A placement was a translation and a uniform scale, so every instance of a shared geometry sat not merely in the same pose but in the same *orientation*. Invisible while a speck was a ball; impossible to miss once it is a flake. `VecOp::Rotate` joins `Scale` as a second unary op — affine in the point, so a node that turns keeps its exact form and stays instanceable — and `Placement`/`Instanced` carry an orthonormal `R` beside the uniform `s`. The pair stays a pair rather than collapsing into one 3×3 because that is the contract: `Rᵀ` and a uniform scale are exactly what let a ray test enter local space without rescaling `t` or repairing a normal.

  **It also closed a disagreement nothing had noticed.** `materialize_mesh()` has oriented scattered items since the sprinkles complaint — each site gets its tangent frame plus a deterministic in-plane turn — while `cook()` kept only `site.position` and dropped the frame. Two renderings of one scene, differing on every scattered object. The reason it survived is plain: **`cook()` and `Instanced` had no tests at all**, so nothing ever compared the two paths. There is now a vertex-for-vertex agreement test, and it fails on 432 of its 443 assertions against the old code. The same fix corrects a second latent error in the same expression — the site position used to be added *outside* the node's scale.

  Raising the particle count is the next thing this unblocks; the frame is still linear in count (19 800 → 1.94 s, 59 400 → 4.27 s), but a flake is now one leaf rather than a dozen triangles.

- **The reasoning under those, kept because the decisions are still load-bearing.** The sub-sections below were written while the work was open; they record why each call was made and stay valid now that it has landed. A `Field` is an expression tree in which an opaque callable is a *leaf* rather than an alternative to the tree: under "tree or callable" a user picks a representation at the declaration site and silently forfeits caching, lowering and hoisting, whereas under "tree with opaque leaves" composing a structural field with an opaque one keeps the outer structure and loses only the one leaf.

  **Identity without names, measured 2026-09-15** (`benchmarks/bench_field_identity.cpp`). Sharing anything between fields requires deciding when two are the same function. The demo installs **39 611 field instances from 6 distinct closure types** — 19 800 dust motions and 19 800 dust colors each come from a single lambda written once inside a loop. By address they are 39 611 different functions; the discriminator that sees them correctly is `typeid` of the closure type, which is compiler-generated and so costs the user nothing and opens no door to registering or versioning function *names*.

  **The deduplication rule, and why it is narrower than identity.** Type identity alone is sound only for a *stateless* closure, where every instance really is interchangeable — `std::is_empty_v` decides this without asking the user, and the demo's 8 `grow_scale` instances are exactly that case. For a capturing closure the state that makes instances differ is hidden inside the object, so deduplicating on type would silently fuse 19 800 particles onto one trajectory: a correct answer under the assumption that the capture is empty, wrong otherwise, which is the same shape of defect as a NaN walking through a filter phrased to skip rather than to keep. So: deduplication is promised for stateless leaves, and for parameterised leaves on `(type, params)`; **a capturing leaf with empty params is opaque and is never deduplicated.** Wanting deduplication means putting the data in `params`.

  `params` is a flat POD map — numbers, strings, node references, arrays — and deliberately not nestable. Nesting would be a second IR living inside `params`, which is precisely how an IR becomes a language that has to be learned.

  **Where names stop being optional.** `type_index` does not survive serialization: once a scene reaches a browser as JSON there is no closure type anywhere in it, and none can be reconstructed. `io/build_json.hpp` is therefore the first real consumer of function *names*, not a hypothetical second contributor, and the date is already on the calendar even though it is not today. Recorded here so it does not get re-litigated from scratch: names are refused now for want of a consumer, not on principle.

  **Representation: a flat pool in topological order, not a linked tree. Settled 2026-09-15.** A `Field` is logically a tree; its representation is a `vector` of tagged ops addressed by index, children before parents, with a `Field` being a *range* in that pool. The reason is not the cheap structural hash and not the cheap deduplication — those are consequences. The reason is that a topologically ordered array evaluates in **one linear pass**: a loop over indices in which each node reads children that are already computed. No recursion, no explicit stack, no branching on depth. That is the shape a GPU can run, and it is the shape a linked tree cannot provide without a flattening walk first — so building the tree first would mean writing the flattening anyway.

  This buys the same question `Trace` already answered one level up: what is substitution, `f(g)`? In a linked tree it is a pointer swap, local and free. In a pool it is either a copy or shared structure. The answer is the one `compose()` already uses for nodes — **append, never mutate**: substitution pushes new ops onto the end of the pool and returns a new range, leaving every existing range valid. One principle applied at two levels, the same way type erasure and `chart_of` were.

  **The report must separate recognized from unknown.** A tally that says "39 600 fields, 1 distinct type" is indistinguishable between "they really are all one closure type" and "they all fell into one unknown bucket", and that is the class already named in `conventions.md` — a value wearing the costume of an answer. `field_report()` reports three numbers, not one: recognized, unknown, total. This is not hypothetical; a probe written while measuring the report's cost produced exactly that ambiguous line, counting `move_only_function` fields it could not classify alongside `std::function` fields it could.

  **Cost of the report, measured 2026-09-15.** A full walk classifying every field costs 2.8–5.1 ns/field — 112–204 µs at the demo's ~39 600 fields, against an 80.7 ms frame, i.e. 0.14–0.25%. So `field_report()` is called per frame and needs no incremental accumulation; the concern that an O(all fields) report would be too expensive to call, and would therefore die unused like the orphan mesh headers, does not survive the measurement.

  **Why `typeid` has to be captured by `Field` itself.** `PointField` is today a bare alias for `std::move_only_function`, which unlike `std::function` exposes no `target_type()` — so a motion field in the current tree cannot be counted by type, let alone deduplicated. That is a fact about an alias we do not own, not a property that blocks anything: the fix is `typeid(F)` captured in `Field`'s constructor, one line, and it works because a closure type is one type for every instance of one lambda (see the identity measurement above). It is recorded here rather than done because there is no `Field` constructor yet to put it in — this item is what creates one.

  **Ownership, settled 2026-09-15 and worth stating because the obvious phrase is wrong.** `Field` owns its pool: it is a value type, copies are deep, and `+=` appends into its own vector where nothing else can observe it. "Append, never mutate" is the right phrase for a *shared* pool — a `Field` that is a `(shared_ptr<Pool>, range)` and must keep other fields' ranges valid — and carrying it over to the owning version describes a guarantee the code neither needs nor gives. Here append-only protects exactly one thing, topological order, because a node can only name indices that already exist. If substitution across many fields later wants the shared form, the aliasing guarantee has to be added and proved, not inherited from the phrase.

  **Two cracks, known and not being solved now.**

  `typeid` gives one type per lambda *expression*, not per lambda body. Three `[noise]` lambdas written inside one loop are one type; three written in three places with identical bodies are three. That is correct for deduplication — different sites are different closures — but a report saying `types=3` cannot say *why* there are three, and the person asking "why did deduplication not fire on my scene" needs the report to answer rather than a person. Same seam as `unknown`.

  **Instancing: what is actually on the table, measured 2026-09-16 rather than estimated.** `cook()` on the donut's dust shape finds **2 shapes for 19 801 objects** — 159 552 vertices' worth of geometry if each object carried its own, against 1 160 actually stored, **138× fewer**. And it reports that **19 800 of those 19 801 are opaque-refused**: the sharing opportunity is total, and every instance of it is blocked by one thing, the motion being written as a deformation when it means a placement.

  The upper bound on what that is worth was then measured directly, by rendering the same frame with the dust removed: **2.02 s with, 0.895 s without — the dust is 56% of the frame.** That bound is what says the work is worth doing, and it is also what says *which* work: the ~45 ms of per-frame vertex building is the small part, and the large part is that 396 000 triangle leaves would become 19 800 instance leaves over one 12-triangle shape, a 20× reduction in what the BVH holds and traverses.

  Which means the three layers are not equal. Reading placements out of structural motions needs a scale op and a `t`-dependent leaf (`VecOp::Const` is a literal, so the cheap route is not available); extracting them in `cook()` is straightforward after that; but only the third — a per-leaf transform in the BVH — touches traversal, and `ray_hit(ray, shape)` takes no transform today, so that one is a change to load-bearing code rather than an addition. It is the one carrying most of the 56%.

  Layout. Once the threshold is payload against L3, the next question is arrangement: `FieldOp` carries a `std::function`, a `type_index` and a payload size inline, so walking a pool to evaluate it drags the opaque-leaf metadata through cache along with the arithmetic. Splitting the pool into arithmetic and leaf-metadata arrays (SoA) is the natural continuation of the same threshold, and it is not worth doing before there is a lowering pass that walks the arithmetic alone.

  **Scrubbing time backwards already works, and is a property rather than a feature.** `materialize(trace, idx, t)` is a pure function of `t`: a motion field takes `(p, t)` and the trace never stores a previous frame. So rewinding is passing a smaller `t`, with no history buffer and nothing to record — the thing games with time control build machinery for falls out of the trace being a description rather than a state. Worth writing down because it is easy to lose: it holds *exactly as long as motion stays a pure function of time*. The moment a solver writes a pose that depends on the previous one, rewind needs history again. That is the strongest practical argument for keeping fields pure and putting any statefulness somewhere it can be named, and it is an argument that costs nothing today.

  **BVH under mutation: build-once, no refit.** `spatial/bvh.hpp` builds by SAH and has no refit or incremental update, which is exactly the weak spot a destructible scene would hit — a voxel grid updates one cell, a BVH has to rebuild the tree. This is not a defect to fix speculatively; it is the shape of the thing, and the interesting angle is ours rather than generic: an exact analytic leaf *replaces* leaves, and rebuild cost scales with leaf count, so the same measurement that showed a 7 000× build-time difference in the DSL row is also the argument that exact leaves make a rebuilding scene cheaper, not just a static one. If a mutating scene is ever a real consumer, refit comes first and exact leaves make refit smaller.

  **Raised again 2026-09-17 as "optimise the BVH across updates", and the answer is still "not yet", now with the reason being a number rather than a preference.** The donut rebuilds both trees every frame, so an animation pays the build once per frame and `--video` pays it 750 times. Refit — walk the tree bottom-up, widen each box, leave the topology alone — is the standard answer and is `O(n)` with no allocation. What makes it the wrong thing to reach for *here* is that **instancing already took the build down to 19.8 ms against a ~2 s frame**, so refit is optimising about 1%. The 16× that instancing won on build was the large factor, and it has been collected.

  Two things would change that, and they are worth naming so the item is not re-raised from scratch a third time. **Traversal getting cheap** — shadows, textures and refraction all add rays, and every one of them amortises the build further, so build share *falls* rather than rises as the renderer grows. **Mutation within a frame** — a physics step that moves bodies many times per frame pays the build per step, and `ball_pit_demo`'s `step_physics` is exactly that consumer, already written and already not using a BVH.

  There is also a correctness-shaped caveat specific to this scene, and it is the interesting half. Refit never makes a tree *wrong* — the boxes still contain their geometry — it makes it *bad*, because the topology was chosen for where things were. The donut's dust is close to the worst case imaginable for that: particles burst outward from a single point, so a tree built at `t=0` has every leaf in one cluster and by `t=3` describes a cloud it no longer separates. So refit here cannot be unconditional; it needs a rebuild trigger, and the honest one is measured rather than a constant — track the SAH cost of the refitted tree against the cost it had when built, and rebuild when the ratio crosses a threshold. That is a number `BVH` can report about itself, which is a better thing to add than a refit that silently degrades.

  **`cook()` and `Cooked<T>`, settled 2026-09-15.** The phase boundary the
  DSL already had in practice and never expressed: a `Trace` is mutated
  while it is built, then read-only while it is rendered, and nothing in
  the type system said so. `trace.cook()` names it.

  **A type, not a flag.** A `frozen_` bool would be a value that looks
  like a guarantee and does not hold one — nothing stops `trace.torus()`
  after it is set, and no compiler notices. That is the class named in
  `conventions.md`, applied to our own design. `Cooked<T>` holds the
  guarantee structurally: it simply has no `torus()`, `offset()`,
  `scatter()`. Mutation after cooking is not discouraged, it does not
  compile. Same shape as `Trace` describing and `Placed` evaluating,
  one level up.

  By value, and no `thaw()`. `Cooked` owns what it holds rather than
  pointing into a `Trace` that has to outlive it — the lifetime coupling
  is exactly what the owning-pool decision in `Field` rejected one level
  down. Going back to an editable trace means editing the source and
  cooking again, and if that is ever made a method it costs an O(n) copy
  and should say so, because the alternative is a flag under a new name.
  Note that "rebuild means a new trace" is *today* forced rather than
  chosen: `TraceNode` holds a `move_only_function`, so a `Trace` cannot
  be copied at all. Once motion is a `Field` that owns its pool, the
  constraint lifts and this becomes a decision on its merits.

  **What cook() does is one thing, and that is the argument for it.**
  Five loose ends that were being tracked separately turn out to be the
  same work, all of which can only happen once the whole scene is known:
  expanding operations into objects; deduplicating shapes (which *is*
  instancing — "one shape, N transforms" cannot be formed before you know
  there are N); compacting the pool; the field report, computed once
  instead of per frame; and AoS→SoA. One well-placed boundary closing
  several problems is the signal that they were one problem.

  Deduplication is part of the expansion, not an optimisation on top of
  it — without it there is no instancing at all. SoA and compaction are
  optimisations; those two are not the same kind of thing and the
  distinction is worth keeping.

  **Objects need their own index space, and `Scatter` proves it.** A
  trace node is an *operation*; a scene object is a *result*. One
  `Scatter` node materialises 600 sprinkles, so no node index names the
  37th. `Compose` is grouping and not an object at all — it has no
  `Placed` of its own — which is why `.moving()` on one is refused rather
  than implemented.

  **One layout, three sources of value.** Instancing and a runtime share
  a shape: a table of shapes plus N per-instance slots. They are not one
  mechanism, because the slots are filled three different ways, and
  collapsing that difference would hide the reason authored motion cannot
  serve physics:

  | source | when the value is produced | storage |
  |---|---|---|
  | instancing | once, in `cook()` | the transform, computed and kept |
  | authored motion (`.moving()`) | freshly, every frame | none — a pure function of `(p, t)` |
  | physics | accumulated across steps | the integrator's state |

  A pure function of time and an integrator are different mathematical
  objects and neither reduces to the other. `.moving()` cannot be
  extended into physics; physics is a second branch.

  **The consumer for all of this already exists and is not connected.**
  `examples/ball_pit_demo.cpp`, `tumbling_body_demo.cpp` and
  `native_collision_demo.cpp` include no `io/build.hpp` at all — the
  physics in this repository lives entirely outside the DSL. And
  `ball_pit_demo.cpp:349` is `step_physics(std::vector<SphereBody>&
  bodies, ...)`: a flat array of objects carrying mass, velocity and
  position, mutated every substep. That is a runtime, already written,
  already working, answering questions (1), (2) and (5) above before they
  were asked.

  So the next piece of work is not "design a runtime" but **express
  ball_pit's scene through the DSL and record where it breaks**. The
  breakages are the requirements, found rather than invented, and the
  demo is a deliverable either way. Questions (3) and (4) will not be
  answered by it — ball_pit uses a fixed `dt`, no event timers and no
  metric — and that is correct: a consumer answers the questions it
  actually has.

  **Five architecture questions, open and deliberately unanswered.** `PointField`'s eventual shape depends on them, and none has a consumer yet, so answering them now would be designing blind — the same call as `SurfaceWithBoundary`. Recorded so the questions are not rediscovered as surprises: (1) where mutable state lives, if `Trace` and `Cooked` are both immutable; (2) what the unit of mutation is — rebuilding a trace, or writing into a parallel structure; (3) that "time" today names two different things, an object's local clock and the metric's dilation, which are items 5 and 7 and not the same mechanism; (4) whether a timer is a function of `t` or state, which decides whether an object can wait, pause or start on an event; (5) whether physics bodies live in the trace or beside it.

  What follows from them *now* is only this, and it is cheap: a field's input should be a **named environment**, not a bare `t`. Today it holds one thing and behaves exactly as today. Tomorrow a local clock or a solver-written pose is a field added to a struct rather than a change to every signature and every call site. Item 3 taught the same lesson from the other side — the slot was right and what filled it changed.

  **The next crack, unsolved on purpose.** Two different `particle_motion`-shaped functions in one trace with different parameter counts — one takes three, the other five. Either the shorter pads with empties, or two genuinely different functions end up sharing a name. No answer yet, and it does not block the work, but it is the same seam and worth expecting.

  Justification is **not** call overhead, and not hoisting either: the indirect call measures ~3.3 ns/vertex (~4% of a realistic callee), so even eliminating it entirely is within noise. Hoisting `t`-only subexpressions is worth well under 1% on today's demo — `grow_scale`'s one `smoothstep` per vertex — and only becomes material when the hoisted subexpression is expensive, which is the local-time item below. What fields actually buy is the parameter table and deduplication above. Selling them on hoisting would be a fourth refuted hypothesis in this series. The one thing a closure genuinely cannot do at all, rather than merely more slowly, is lower to branchless GPU code.

  Composition becomes substitution of one tree into another's point slot — which is why `.moving()` was made to compose first, so the callable and structural forms agree. Degradation has to be **visible**, not felt: `is_structural()` plus a per-trace report ("8 opaque leaves, deduplication covers N of M fields"), for the same reason the polynomial solvers' lost root is pinned as a `[!shouldfail]` test rather than a comment. Silent degradation that the user discovers through a frame-rate drop is the NaN-through-the-filter failure one storey up. The risk to watch is that an expression IR quietly becomes a language users must learn; a plain lambda has to stay a normal form rather than a permitted exception, and the donut demo is the place that failure would show.
- **[want]** Time as a property of the space rather than a global scalar. A scene's `t` is currently one number shared by every node. Deriving each object's local time from where it sits degenerates to exactly today's behaviour in flat space at no cost, and in a curved one gravitational time dilation falls out of a metric the library already has — `physics/relativity/` ships Schwarzschild and Kerr as substitutable callables. What is missing is small and specific: no proper-time or local-clock-rate helper exists there yet. Deliberately preferred over a conventional keyframe timeline, which was considered and dropped.
- **[want]** Open the node set — `Kind` is a closed `enum class`, so a user cannot add an operation without editing `build.hpp`. `io/scene.hpp` already solved the same problem in the same namespace with an open registry (`register_shape_kind()`), so this is a consistency gap, not an unknown. Sequenced after the field work because whatever opens `Kind` has to survive lowering.
- **[course]** Vulkan live display — CPU raytraces (`donut_demo.cpp`'s `--photo` engine, extended), a small new Vulkan path just presents the frame (swapchain + one texture, uploaded and blitted every frame), not `viewer::App`'s mesh/point-cloud rasterizer. Real new plumbing (no sampled-texture descriptor/pipeline exists today), deliberately deferred rather than landed blind against a deadline.
- **[want]** Texture/UV mapping — `io::Material` currently has only `base_color`/`roughness`, no texture at all. Flagged back on 2026-09-07 alongside conform-to-surface (now shipped as `offset_surface`) as one of two primitives needed before the DSL; still open.

- **[course]** **A motion written as a plain lambda cannot be seen into, so almost every node in the donut demo is refused as a deformation when it is nothing of the kind.** Found 2026-09-17 by asking what the demo's 22 tessellated objects actually were.

  ```
  scene: 35222 objects -> 35200 instances + 429848 triangles;
         10995 refused (motion deforms), 35200 shared
  ```

  The 22 are eighteen sprinkle `Scatter` nodes (six colours times three bands, one merged object each), plus the cube, the dough and the icing, plus the table — which *is* a placement and goes to triangles only because it carries no closed form. The 10 995 counts refused *instances* rather than objects: 10 992 sprinkles and three whole nodes.

  The cause is not that those motions deform. It is that they are written as plain `(p, t)` lambdas:

  ```cpp
  auto grow_scale = [](const Vec<double, 3>& p, double time) {
      return Vec<double, 3>{p * e};        // a textbook placement
  };
  ```

  An opaque leaf that touches `p` sets `reads_point`, and `is_placement()` then answers "deformation" — correctly, since it cannot ask a closure what it does. Written structurally as `scaled(point(), e)` the very same motion is a placement. The donut's dust was rewritten that way when instancing landed; **nothing else was**, so everything but the dust is refused for a reason that is a matter of spelling.

  So sprinkle instancing has **two** preconditions, not one: the axis binding above, and structural motions here. Both are fixes rather than features, and they land together.

  **And rewriting the demo's motions does not close this**, which is the part worth separating. The demo is one user; the trap is in the API. Someone who writes `p * e(t)` as a lambda — the shape `.moving()` itself advertises, `(point, t) -> point` — gets a silent refusal, no instancing, and no way to find out why. `docs/conventions.md` records the general form under "the variant where the value is honest and the *name* is not": a predicate that can only see a syntactic property must not carry a name promising a semantic one. Three fixes, not exclusive:

  - **The refusal explains itself.** `VecField` can distinguish "refused because an opaque leaf reads the point" from "refused because the expression genuinely is not affine", and the first should say so with the remedy attached. `Cooked::refused_nodes()` already exists precisely so a report can point at nodes rather than state a number nobody can act on; it should also point at *why*.
  - **`is_placement()` is renamed to what it computes.** `is_recognisably_a_placement()` is ugly, and the ugliness is the point: it makes a caller ask "recognisably by whom?".
  - **The structural spelling becomes the obvious one.** `scaled(point(), e)` is already shorter than the lambda; it is simply undiscoverable, because the slot advertises a callable. That is an API problem and not a documentation one.

- ~~**[want]** Does collapsing the dust into one node remove `content_hash`'s O(vertices) cost?~~ — **measured 2026-09-17, and the answer is that the cost was never there.** The premise was wrong twice over. The dust is not a `Literal`: `flake()` builds a `Space` node, hashed by its exact form's type plus sixteen chart samples, so the O(vertices) path touches exactly two nodes in the whole scene (the cube and the table). The "158 000 vertices hashed across the dust" figure was inherited from when a speck was `literal(dust_speck(...))` and never re-derived.

  And the measurement settles it regardless of the premise:

  ```
  trace nodes   35 231 -> 34
  cook()         114.6 -> 105.6 ms
  ```

  A thousandfold fewer nodes and **8% less time**. So neither "proportional to the node count" nor "stronger than that" — hashing was not the cost at all. What dominates is per-*object* work, which the collapse does not reduce: 46 196 objects still each evaluate a placement and resolve a material through two opaque fields. The first line to look at, if cooking ever needs to be faster, is that — not the hash.

  The real win is the one the entry was not looking for: the trace went from 35 231 × 456 B ≈ **16 MB** to 34 × 456 B ≈ **15 KB**, which is what makes two million reachable at all.

- **[course]** **An exact form that does not agree with its own chart, created that way by the factory.** `docs/api-reference.md` already states the invariant — *"the exact form must agree with the map, so any operation that can move them apart clears it"* — and `.moving()` honours it. A factory can violate it at the moment of creation, and one does.

  | factory | exact form | chart it is recorded beside | agree |
  |---|---|---|---|
  | `torus(R, r)` | `Torus{R, r}` | `make_torus` | yes |
  | `cylinder(r, h)` | `BoundedQuadric::cylinder_z(r, 0, h)` | `make_cylinder` — lateral surface, no caps | yes |
  | `sphere(r)` | `BoundedQuadric::sphere(r)` | `chart_of(Sphere<2,T>)` | yes |
  | `flake(half, bulge)` | a sphere clipped to `Box{-half, half}` | a spherical cap out to `v_max`, derived from `rim = min(hx, hy)` | **no** |

  `flake`'s chart is the disc *inscribed* in the clip box; the exact form keeps the box's corners. For the donut's `flake({0.010, 0.010, 0.003})` the square's corner reaches 0.0141 against the disc's 0.010 — so the exact form is a rounded square and the tessellation is a circle. The exact form covers **more**, which is the opposite direction from the guess that prompted this entry.

  **Fixed 2026-09-17, and the measurement corrected the diagnosis twice on the way.** The disagreement is not in the corners — sampling both shows x and y agreeing to the bit — it is in **z**: the chart spans `[-0.001431, 0.003]` and the old clip claimed `[-0.003, 0.003]`. Below where the cap ends the sphere keeps widening past the slab's half-width and is cut by its sides, so the exact form carried a skirt the tessellation had never heard of.

  Making the chart match was not available: `ParametricSurface`'s domain is a rectangle in `(u, v)` and "sphere ∩ box" is not one. Making the clip match was, and is exactly right — clip to `rim` and to the height where the cap reaches it, and then at any `z` above that the sphere's radius is at most `rim`, the sides never cut, and the bottom cuts precisely where the cap ends. Sphere ∩ box is then *equal* to the cap. The clip stops describing the size the caller asked for and starts describing the region the surface occupies, which is what a clip is for.

  **The visual consequence, decided rather than absorbed.** Removing the skirt makes the dust finer: on a dust-heavy frame 10.5% of pixels change and 95% of those get brighter, which is less speck covering more background. Kept, for three reasons in order of weight. The skirt was surface the tessellated form never had, so keeping it would keep a disagreement between two renderings of one node — the thing this fix exists to remove. Finer specks read as dust rather than as gravel, which is the direction the scene has been moving anyway. And if the cloud wants more visual mass, the honest lever is the flake's `half`, not a loose clip. A frame hash cannot referee this: the representation changed honestly, so the pixels *should* differ, and only an eye can say whether the new shape is the wanted one.

  Guarded by a test over every factory that records both parts, comparing the sampled bounding box of the chart against the exact form's own. Bounding boxes rather than surfaces, deliberately: cheap, dependent on nothing that can itself be wrong, and it catches the failure that matters — one of the pair covering ground the other does not. A divergence is a bug until someone writes a documented exception, and there are none. Checked against the old clip before being trusted: four assertions fail there.

  The reason this is worth an entry rather than a fix-in-place: the invariant exists, is documented, and is enforced on the one operation that was thought to threaten it. Nothing checks it where the two are first written down together — which is also where it is cheapest to check, since both are right there.

  **Withdrawn along the way:** the claim that this class was demonstrated by `cylinder()` losing its end caps. It is not. `make_cylinder` is the lateral surface with no caps and `BoundedQuadric::cylinder_z` matches it exactly. The donut's sprinkles became open tubes because the demo swapped its own capped `solid_cylinder` helper for `scene.cylinder()`, which is a different *shape*, not a different *representation* of the same one.

- **[course]** **`scatter_frame` builds its frame from one normal, and orienting *along* a surface needs a full one.** Found 2026-09-17, by checking a claim rather than by hitting a bug: the plan said the donut's sprinkles would instance for free once the renderer read a cooked scene, since `cylinder()` records a `BoundedQuadric`. They did not, because the demo does not use `cylinder()` — the sprinkle is a `literal()` of a hand-built mesh whose axes have been permuted so its length points sideways.

  That permutation is the finding. A frame derived from a normal is enough for anything that points *along* the normal, or that has no orientation worth speaking of. A sprinkle lies *across* the surface, and there is no way to say so: `scatter()` binds the item's local z to the normal and offers no alternative, so the only remaining move is to rotate the mesh by hand — which a closed form cannot follow, since the exact `BoundedQuadric` is a cylinder about z and stays one. Hence `Literal`, hence no instancing.

  Two distinct gaps sit under that, worth separating because they block different things:

  - **The in-plane directions are arbitrary.** `basis_from_normal` (Duff et al. 2017) picks *some* pair perpendicular to the normal, and `scatter_frame` then turns it by a hash of the site index. Perfect for sprinkles, which really do lie every which way. Useless for anything that should follow the surface itself — hair along a curvature direction, scales along a flow, tiles along a parametrization.
  - **The axis binding is fixed.** Local z goes to the normal, always. This is the one the sprinkle is blocked by, and it is the cheaper of the two.

  Until a full frame exists — normal, tangent, bitangent, with the tangent meaning something — scattered items that must lie along a surface stay `Literal` meshes and cannot be instanced. The ceiling that puts on the donut demo is concrete and measured: 11 000 sprinkles are 429 848 triangles that a closed form would have made 11 000 instances of one shape.

- ~~**[want]** Scattered items lying fully flush to the target surface, not slightly proud of it~~ — **done 2026-09-17**, and the cause was not the one this entry assumed. They were not proud because of a tangent-plane approximation; they were proud because `offset()` silently rendered at 48×24 instead of the tessellation it was asked for, and a coarse mesh of a convex surface sits *below* the analytic surface that `scatter()` places against. The gap was lifting every item into view by accident. `scatter_lift()` now raises an item by the depth of its own geometry below its local origin, and `scatter()`'s new `seat` parameter says how deep it sits — 1 rests it on the surface, 0 puts its origin there, between is a press into something soft. Note the trap that made a named parameter necessary: since the lift is derived from the item's own lowest point, shifting the item's mesh down by *d* lowers `min_z` by *d* and raises the lift by *d*, so the two cancel exactly and the nudge does nothing. The donut demo had exactly such a nudge, with a comment explaining what it was for, and it had silently stopped working.

- **[course]** **Per-instance parameters: `MotionEnv` gains the instance's own origin.** The blocking step for two million particles, analysed 2026-09-17, and the analysis is the valuable part because the obvious answer is the wrong one.

  The obvious answer is "carry the exact form through `cook()` so a `Scatter` can be instanced". That is true and necessary and it is **not** what stops the dust, because the dust is not a `Scatter` at all: `Scatter` places an item on a *surface*, while dust flies through space with a per-particle burst direction, target letterform and swirl seed. The recorded measurement "one node, N instances → trace = 3 nodes" was about sprinkles.

  Restructure the dust as a scatter over an invisible sphere — so that `site.position` *is* the burst direction — and the real wall appears one line further in:

  ```cpp
  pl = n.transform.placement_at(MotionEnv<T>{Vec<T, 3>{}, t});   // cook(), once per node
  ```

  **A `Scatter` node has one motion field for all its instances.** For particles to fly differently the field must tell them apart, and the only thing it can tell them apart by today is `p`, the vertex position — which sets `reads_point`, makes the motion a deformation, and makes instancing refuse. Correctly: that test is doing its job. Per-instance variation smuggled through the vertex position is exactly what it exists to catch.

  The fix is the one `io/field.hpp` already wrote down the reason for — *"tomorrow a field is added to this struct and nothing else moves"*. `MotionEnv` gains the instance's **origin**, its rest position. A field reading `env.origin` stays a placement, because an origin is one value per object rather than one per vertex, so the affine test still holds and the instance path still applies. The dust then becomes a single node with N instances, each on its own trajectory.

  **A second reason, found 2026-09-17 while trying to avoid needing this at all.** The plan had a cheaper intermediate step: a purely *radial* burst needs no per-instance parameters, because `A(t) + p·s(t)` is affine in the point and each particle's trajectory falls out of its own site position. The instancing half of that is true — `is_placement()` returns true and the objects collapse onto one node. The animation half is not.

  A placement is `p ↦ b + R·s·p` applied to every vertex, and a scattered instance's vertex is `seat + F·local`. So:

  ```
  particle centre = b + R·s·seat     grows with s
  particle radius = s·|local|        grows with s, the same s
  ```

  One scalar drives both. Expanding the cloud tenfold inflates every fleck tenfold: not a burst, an inflating balloon. They cannot be separated, because `seat` and `local` reach the field already summed into one point — and separating them is exactly what an instance origin is.

  So `MotionEnv::origin` is not only what swirl and letterform convergence need. It is what a *plain radial burst at constant particle size* needs, which removes the cheaper intermediate step from the table: building it would mean building something thrown away one step later.

  **A prediction, registered before the measurement rather than after.** Stage 1 gave 36.4× on stored vertices and 11% off the frame, and the gap between those two is itself explained: traversal is logarithmic, and `log(46 000)` against `log(429 000)` is almost nothing. So at two million the expectation is **7–10× on storage and 10–15% on the frame** — storage *falls* because one node's worth of geometry is amortised over vastly more instances while the instance array itself grows, and the frame moves little for the same logarithmic reason. Written down here so that a frame growing by more than that is a signal to look for a cause rather than a number to accept.

  Scope, so this is not mistaken for a parameter change: `MotionEnv`, a new leaf factory in `VecField` alongside `opaque_of_time`, `is_placement`/`affine_in_point`, `cook()`, `materialize_mesh()`, `Shape::exact`, and the renderer moving from `materialize()` to `cook()`. That last one uncovers a further gap — `cook()` hands back *rest* geometry for a refused (deforming) object with no way to recover its motion, so a renderer driven by `Cooked` alone would draw the exploding cube unexploded. ~~Cost at two million, measured by extrapolation rather than guessed: roughly 350 MB for instances plus tree~~ — **wrong by a factor of three, corrected 2026-09-17 by measuring `sizeof` instead of recalling it**, and the correction is the point of checking before rather than after:

  ```
  sizeof(Object<double>)     192 B  ->  366 MB at 2M
  sizeof(Instanced<BQ>)      112 B  ->  214 MB
  the BVH's own copy of them         ->  214 MB
  nodes_, reserving 2N               ->  244 MB
                                         ~1.0 GB
  ```

  plus about 160 MB of build temporaries (`boxes`, `centroids`, `indices`). The 350 MB estimate forgot `Object<T>` entirely, and forgot that `BVH::build` takes its shapes **by value** — so the instance array exists twice.

  **The full breakdown, measured 2026-09-17, since a figure without its layout is not a figure.** `Object<double>` is 192 B: `shape` 8 at offset 0, `source_node` 8, `translation` 24, `scale` 8, `rotation` **72**, `instanceable` 1 at offset 120, `material` **64** at 128 — members summing to 185 with 7 bytes of padding. That padding is *not* reclaimable by reordering: 185 rounds to 192 under 8-byte alignment in any permutation. A `BVH` node is `{Box<3,double> bounds, uint32 first, uint32 count}` = 56 B, of which the bounds are 48.

  What that makes available, in order of return, with the projected `Object` alongside:

  | change | saves at 2M | `Object` after |
  |---|---|---|
  | ~~`rotation` 72 → 32 (quaternion)~~ — **done**, `Object` only; measured 192 → 152 B | 80 MB | 152 B ✓ |
  | `material` 64 → a `uint16` index into a palette — **measured first, and it is not what the estimate assumed**, see below | ~97 MB | 88 B |
  | `shape`, `source_node` → `uint32` | 16 MB | 80 B |
  | BVH bounds in `float`, rounded outward — a bound only has to be conservative | 96 MB | — |
  | `prim_indices_` `size_t` → `uint32` | 8 MB | — |
  | the BVH owning its shapes: hold indices into `Cooked` instead | 214 MB | — |

  **How the rotation change is staged, decided before it starts.** The rotation is a matrix the whole way down — `scatter_frame` assembles columns from a basis, `rotated()` produces one through `SO3::exp`, their product is one — so a quaternion can only come from `from_matrix`, and the round trip is unavoidable. Measured over four thousand real objects: **zero bit-exact**, worst element 1.72e-15, about 8 ulp. A byte-identical frame after the matrix goes is therefore impossible rather than unlikely.

  What that does *not* explain is any visible change, and the first two attempts to explain one both failed the same way — by sounding like a mechanism without ever being given a number.

  The donut's worst object radius is 10.3 world units, so the worst vertex displacement is 1.8e-14 against a pixel of 5.4e-3: eleven orders below one. *The picture cannot shift.* The second story — that exact ties at a silhouette flip instead — is arithmetic too, and it does not survive it either. The band in which a tie could fall the other way is 3.3e-12 of a pixel wide; a 960×720 frame at four samples a pixel is 2.8 million rays, of which perhaps 5% land on a silhouette, so the expected number of flipped ties is about **5e-7**. Not "occasionally" — not once.

  So the honest prediction is the stronger one: **step two should be byte-identical as well.** If it is not, neither of those stories covers it and the cause is something else — so the response is not a third explanation but a count: how many pixels differ, and where. One on a silhouette is a coincidence to accept and close. A dozen scattered across the frame is systematic, and no displacement of 1e-14 accounts for it under any reading.

  **Both steps came out byte-identical, including the second one** — the strengthened prediction held, and the picture did not move by a pixel even though every rotation now arrives through a quaternion.

  Four preconditions were checked rather than assumed, and each could have made the result meaningless:

  - **The renderer is deterministic.** Three runs of one binary, one hash. Without that, "one pixel differs" could not have been attributed to anything.
  - **The comparison can fail.** Substituting an identity quaternion for one object in forty changes the hash. A byte-identical result from a comparison that cannot detect a difference is not a result.
  - **The expansion stayed cold.** `to_world` runs per vertex and now takes the matrix as an argument; the quaternion is unpacked once per object, beside the instance that keeps it as a matrix. A caller writing `o.rotation_q.to_matrix()` inside the vertex loop would have spent the memory and bought nothing.
  - **The object count is asserted separately from any picture.** A lost rotation makes instances coalesce, and a slightly denser cluster is not a visible defect — neither a hash nor a pixel diff distinguishes forty thousand instances from forty thousand with two on top of each other. The test checks counts, pairwise distinct positions, pairwise distinct orientations, and that none of them is the identity.

  So, two steps and two rules:

  - **Step one: store both, and the matrix stays authoritative.** The renderer reads the matrix; the quaternion is carried and compared. A byte-identical hash then means something, because the thing being drawn has not changed — it tests the plumbing and nothing else. Two sources of truth for one rotation is otherwise exactly how they drift apart in silence.
  - **Step two: remove the matrix and introduce the pixel comparison in the same change.** Splitting them leaves a window where the matrix still exists but the renderer already reads the quaternion and nothing checks the difference. The comparison is what makes the removal safe; without it the removal is blind. It reports a count and locations rather than a pass/fail tolerance, because the expected count is zero and anything else is a lead rather than a threshold to clear.

  **The palette was costed before being built, and the costing changed it.** `Material` carries values rather than functions — `base_color`, `roughness`, `emissive`, `opacity` — so a palette is possible in principle; the question is whether per-object materials *collapse*, and that is a property of the scene rather than of the idea. The demo now reports it:

  ```
  t = 1.5    9 762 distinct materials over 46 196 objects   (21%)
  t = 3.9   13 188                                          (29%)
  ```

  Neither the two-to-five the estimate imagined nor the N that would kill it. Three things follow, and the mechanism matters more than the ratio:

  - **The collapse is saturation, not sharing.** `dust_color` clamps at `d/1.4` and the emissive at `d/0.75`, so every particle past those distances lands on exactly the same grey and exactly zero emission. The particles are not sharing a material because they mean the same thing; they are sharing one because two `clamp` calls flattened them.
  - **The saving is smaller than the 124 MB claimed**, which was computed as though the palette were tiny. At 21% on two million: 128 MB of inline materials becomes 4 MB of indices plus 27 MB of palette, so about **97 MB**.
  - **It is fragile, and cheaply broken.** Replace either falloff with something that does not saturate — an exponential, say — and the collapse disappears entirely, the palette becomes size N, and the change *costs* memory rather than saving it. The 97 MB currently rests on two `clamp`s in a demo.

  The ratio also moves with time (9 762 against 13 188 on two frames of one scene), so a palette built in `cook()` is a different size every frame. Not a defect, but not a constant to plan against either.

  If it is built: `Cooked`'s materials are immutable like everything else `cook()` produces, and that has to be said out loud rather than discovered — sharing one entry across N objects means recolouring one object recolours all of them, and someone will eventually want to flash a single instance. And the test that matters is not "the palette works" but "the palette *collapses*": two objects with equal materials must land on one index, asserted directly.

  **`Instanced<S>` is excluded from that list on purpose**, and `docs/conventions.md` carries the principle: compaction belongs to cold storage, and the hot path keeps whatever computes fastest. `Object`'s rotation is read once a frame; `Instanced`'s is applied at every leaf test the traversal reaches, and a quaternion costs more arithmetic to apply than a matrix. Same saving, traversal untouched.

  The trigger for revisiting it, written down instead of guessed: two million `Instanced` at 112 bytes is 214 MB that traversal walks and no cache holds. *If* traversal at that scale measures as memory-bound, the options are splitting the array — positions apart from rotations, so a ray touches only what it needs — or shrinking the element after all. Against a measurement showing the stall, not before it.

  Landing at roughly `2M × 80 B` for objects plus a much thinner tree. The next boundary after that is SoA — positions in one array, rotations in another — which is a different change and should not be started until the cheap ones are measured.

  **One of them was free and is already taken:** `BVH::build` takes its shapes by value, and the demo handed it an lvalue, so the instance array existed twice for the whole frame. `std::move` at the call site removes one copy — 214 MB at two million — and changes nothing else. Worth recording that checking *what still reads the vector afterwards* caught a real break: the report line counted `insts.size()` after the build and would have silently printed zeros.

  **Two million is a memory problem before it is a speed problem**, and that reframes the work: the reductions available are layout ones the entry below already separates from deduplication. `rotation` as a full 3x3 is 72 of `Object`'s 192 bytes and a quaternion would be 32; `shape` and `source_node` are `size_t` where `uint32_t` would do; `Material` is 64 bytes repeated per instance where most instances of a node share everything but one field. And the BVH taking its shapes by value is a duplication nothing needs. None of that is hard, and none of it should be guessed at either — the number above is what makes it worth doing, and a measured number after each change is what says it worked.

- **[course]** **One object deforming another.** Raised 2026-09-17, from wanting sprinkles to press visibly into liquid icing rather than merely sitting in it. Today a field is a function of `(u, v)` or of a `MotionEnv`, and it cannot see anything else in the scene — so "the glaze closes around each sprinkle" is not expressible, and neither is a footprint, a dent, a contact weld, or a drip running off a placed object. The concrete first case is narrow enough to build: `scatter()` knows its sites and does not expose them, so the target's own thickness field could read them. The cost is what makes it a real item rather than a tuning change — the naive form is a sum over every site per surface sample, 600 × 10 240 on today's donut and around 50 million at the sprinkle counts this scene now uses, so it needs a grid or similar index over the `(u, v)` domain before it is usable at all. Worth naming as the general capability rather than as "dimples": the same mechanism is what a scene needs before objects can be said to interact at all, and it is the one the demo keeps reaching for. Deliberately deferred past the first release.

  **What it breaks, which is not what it looks like it breaks.** Raised 2026-09-17 as a worry that this item ends the "`materialize(trace, t)` is a pure function of `t`" property — and therefore time scrubbing, which falls out of that purity for free. Worth being exact, because the worry is right about there being a cost and wrong about where:

  *Purity survives the dimple.* A scatter's sites are a pure function of (target surface, count, seed) with no history in them, so a thickness field reading those sites is still a pure function of `t`. Rewinding keeps working.

  *What the dimple actually breaks is acyclicity.* The icing's shape would depend on where the sprinkles are, and where the sprinkles are depends on the icing's shape — a cycle in a structure that is a DAG **by construction**, since a node can only name indices that already exist. The resolution is to scatter against the *undented* surface and apply the dent afterwards: one truncated iteration of a fixed point, which is what a solver would do anyway. That has to be stated here rather than discovered, because the alternative reading — make it self-consistent — has no answer, and someone will spend a day finding that out.

  *Purity dies one step later, at the thing anyone asks for next:* a mark that stays. A footprint, a dent that does not spring back, glaze that remembers. That is accumulation over history rather than a function of the current configuration, and it is where rewind genuinely stops — not here. Two different failure modes, one step apart, and only the second one is fatal to the property.

- **[want]** Removing the self-intersection an offset surface produces when its thickness exceeds the local radius of curvature. Seen in the donut's icing as pink fins around the hole, where the drip field pushed glaze down the inner wall of the torus and the offset overran its own centre of curvature. The demo's fix is to stop drips at the inner rim, which is correct for glaze — it runs off the outside, it does not climb into the middle — and is a scene decision, not a library one. The library answer is either a validity check (`offset_shell` refusing, or clamping, a thickness above the base's minimum radius of curvature, which `ParametricSurface` has the derivatives to compute) or self-intersection removal on the resulting mesh. The first is cheap and catches the case at the call that makes it; the second is a general mesh operation this tree does not have.
- ~~Per-instance orientation variance in `scatter()`~~ — **done**, and this entry was stale from the day it was written: `build.hpp`'s scatter branch already gives each instance its own spin about the normal, hashed from `(seed, site index)` so the same trace at the same `t` reproduces exactly. Kept here struck through rather than deleted, because it was listed as open for a week while the code was already in the tree, which is the failure worth remembering.
- **[want]** Generalize `Trace::space()`/`resolve_surface()` beyond `ParametricSurface<T>` to any `Surface` (`ImplicitSurface`, `Sphere`, `Hyperbolic`, ...) — today's `offset()`/`scatter()` chain is hand-specialized to `ParametricSurface`'s `evaluate`/`normal_at`/`area_element`, not the general `Surface` concept.
- ~~Exercise `mesh/conform.hpp` and `mesh/scatter.hpp` in a real demo~~ — **settled 2026-09-15, and the old entry was wrong in a way worth keeping visible.** It said they "compile and pass their own tests"; they had no tests at all, no caller anywhere, and are not reached by the umbrella header either — so their bodies are templates that had never been instantiated, only parsed. Both now have a test in `tests/test_mesh_ops.cpp` that instantiates them and checks the result, and both headers say in their own first lines that nothing calls them and why. They are kept rather than deleted because they serve the case the analytic path cannot — a target that is only ever a mesh — and they are explicitly **not** getting a demo. That is the whole of the maintenance they are owed.
- Explicitly not planned for the scene DSL: dependency-tracked reactive recompute (SolidJS/Vue-signals style). If that gets built at all, `gpu/derive_christoffel.py`'s declarative-derivation → auto-codegen pipeline is the intended home for it (GPU rendering section below), not this graph.

## Polynomial solvers: the degenerate leading coefficient

Design settled 2026-09-15. **Shipped 2026-09-17**, and the three
`[!shouldfail]` pins did exactly the job they were written for: the day
the solvers stopped dividing by a vanishing leading coefficient all three
started passing, Catch2 reported *that* as a failure, and the tags had to
come off. A bug recorded as a comment gets skimmed; one recorded as a
test that must keep failing cannot be quietly forgotten, and cannot be
quietly left marked broken after it is fixed either.

The return type is `UpTo<Complex<T>, N>` in a new `core/up_to.hpp` — at
most N values with a count, no allocation — and the name states the thing
that changed. `solve_linear` exists as a named function because a
degenerate quadratic is exactly that, and each solver now collapses to
the one below it. The donut frame is byte-identical afterwards, which is
the expected result rather than a disappointing one: a ray exactly
parallel to a quadric's null direction is measure-zero among camera rays,
and the defect was always about the cases a renderer does not sample.

**Two decisions the design did not contain, recorded because both had a
plausible alternative.**

*The zero test is exact equality, not a tolerance.* A relative test — "a
is negligible beside b and c" — catches more cases and is the wrong trade:
it turns a genuine but ill-conditioned quadratic into a linear one and
drops a real root that is large but honest, replacing one silent wrong
answer with another. Exact zero is the only threshold that says *this is
not a quadratic* rather than *this quadratic is awkward*, and it is also
exactly what the geometry produces — for a ray along a cone's generator
the t² coefficient is one quantity minus itself, zero to the bit. What
this leaves open, stated so it is not mistaken for covered: catastrophic
cancellation when the leading coefficient is tiny but nonzero. The fix
for that is the numerically stable form (`q = -(b + sign(b)·√disc)/2`,
then `x₁ = q/a`, `x₂ = c/q`), it is a separate defect, and bundling it
would have made one change two.

*`ray_quadric_proximity` checks the answer, not the coefficient.* The
`[!shouldfail]` test's own comment predicted the fix would be a
degeneracy check on the t² coefficient *before* `solve_quadratic`, on the
grounds that checking after would catch the symptom. That prediction was
written when the return type could not say "fewer roots", and it is
superseded rather than merely unmet: asking whether the answer has two
roots is what the function's model actually requires, and it leaves one
place deciding what degenerate means. A coefficient test at the call site
would be a second copy of the solver's rule, free to drift from it.

**The defect.** `solve_quadratic`, `solve_cubic` and `solve_quartic` all
divide by the leading coefficient with no guard. `solve_quadratic(0, 2, -1)`
is a perfectly good equation with the root `0.5` and returns `NaN` and
`-inf`. Because `NaN` reports `is_real() == true`, the garbage then walks
through every caller's filter. Concretely: a ray parallel to one of a
cone's own generators meets the surface exactly once, and the hit is
silently lost.

**Where the risk actually lives, stated as a boundary rather than a file
list.** The leading coefficient is dangerous exactly where it is *data*
and safe exactly where it is *structure*.

- **Data.** A ray's direction is user input and can degenerate — down a
  cylinder's axis, along a cone's generator. `geometry/ray_surface.hpp`,
  `geometry/ray_hit.hpp`, `io/scene.hpp`, `physics/mechanics/rigid_contact.hpp`.
  A third live instance sits at `ray_surface.hpp:394` in
  `ray_quadric_proximity`, and it presents worse than the other two. A ray
  down a cylinder's axis returns `Result` **success** with
  `closest_t = NaN`, a `NaN` point, and `miss = 0` — because the NaN root's
  imaginary part is exactly zero, so `abs()` of it is a perfectly plausible
  number in the one field a caller is most likely to branch on. The other
  two hand back NaN, which is at least conspicuous. This one reports a
  clean grazing hit for a ray that runs down the middle of the tube and
  never approaches the wall.
- **Structure.** A characteristic polynomial is `det(A - λI)` by
  definition; its leading coefficient is 1 as a matter of mathematics and
  no input can change that. `spaces/spd.hpp`, `algebra/eigen_decomp.hpp`,
  `core/precision.hpp` are safe not because they were checked but because
  it could not be otherwise. Those three pass `T{1}` literally.

These are different grades of certainty and the second is the stronger
one. An `assert(leading != T{0})` at the structural sites is still worth
its line, but for a different reason than at the data sites: there it
guards a future refactor changing what gets passed, not an input.

**The return type.** A fixed-capacity container with a count, no
allocation — it must be able to say "one root" where the degree collapsed.

- `size()` is the count; `capacity()` is 2/3/4. Range-`for` takes its
  bounds from the count, so every existing `for (auto& r : roots)` call
  site keeps working unchanged.
- `operator[]` is **unchecked in Release** — no branch in a hot path — so
  indexing at or past `size()` is undefined there. This has to be said
  loudly in the docstring, because a reader arriving from `std::array`
  expects a fixed size where every index is always valid, and that is
  precisely what stops being true. In Debug it carries
  `assert(i < size())`.

  That assert was nearly dropped on the theory that it would fire on
  correct code, `spaces/spd.hpp` indexing `roots[2]` being the worry.
  Checked instead of assumed — `git grep "roots\["` finds exactly four
  library sites. The two in `spd.hpp` are monic, so the count is always 3
  and the assert can never fire there. The other two are both inside
  `ray_quadric_proximity`, where the count really can fall short — and
  they are the ones already producing `miss = 0` today. So the assert does
  not break correct code; the only code it breaks is the code that is
  already wrong.
- `at(i)` returns `Result<Complex<T>>` and **never throws**. Checked
  access signals; unchecked access does not; different contracts, so
  different return types. Exceptions appear nowhere else in this library
  and must not appear here.
- In Debug, slots between `size()` and `capacity()` are poisoned with
  `NaN`. Zero instructions under `NDEBUG`, and under the Debug CI job it
  turns reading past the count from undefined into visibly wrong.

**Tests have to cover the case that changed, not the one that already
worked.** `roots[2]` used to be valid always, by virtue of a fixed-size
array; now it is valid only when the count reaches 3. So both need
covering: indexing within the count, and `at()` refusing beyond it.
Testing only in-range indexing tests the old behaviour.

**The one non-library caller, and the convention question it raises.**
`rsc/include/precision_ops.hpp`'s `flatten_cubic_roots` takes
`const std::array<Complex<double>, 3>&` explicitly and copies three roots
unconditionally, so it needs a signature change and a check.

The check is an `assert`, not a `Result`, and the reasoning is worth
keeping because the first answer here was the wrong one. The registry op
receives `std::span<const double> in` with `in[0]` as the leading
coefficient, which *looks* like a boundary with untrusted input. It is
not. A boundary is where data arrives from outside the process and may be
anything — `load_obj`, `json::parse`. Here the data comes from a harness
living in the same binary, generating inputs from a distribution it
defines itself. That is an internal pipeline invariant, and
`conventions.md` already answers those: `assert`, the same as
`physics/mechanics/`, `spaces/` and `mesh/`.

So: `assert(in[0] != T{0} && "harness contract: monic cubic")` in the op,
`assert(roots.size() == 3)` before the copy, and `flatten_cubic_roots`
keeps returning `void`. A monic cubic always has three roots with
multiplicity; if it has two, either the input or the solver is broken,
and both are bugs rather than failure modes. `Result` here would be
signalling a refusal that cannot legitimately happen.

**`Op::Fn` returns `void`, deliberately, until something needs otherwise.**
An op with no failure channel is the right shape for every op the registry
has, because none of them can legitimately fail. A channel is warranted
when an op appears that genuinely can — an external solver, a network
call, a parse, anything with unbounded input — and not before. Recorded
explicitly so the next reading of "but the convention says `Result`"
does not have to re-derive that the convention is about boundaries rather
than about everything.

**NaN poisoning is a debugging aid, not a channel.** Slots between
`size()` and `capacity()` get NaN under Debug so that reading past the
count is visible rather than undefined. It is the same category as
`_GLIBCXX_ASSERTIONS`: a signal while debugging, never part of the
contract, and zero instructions under `NDEBUG`.

The training checkpoint is **not** affected. `rsc/include/precision_task.hpp`
builds coefficients from three sampled roots by Vieta and returns
`{1.0, b, c, d}` — the leading coefficient is literally 1 in every
sample, the degenerate case never occurs in the task's distribution, and
`solve_cubic_f64 4 4 6`'s pinned in/out sizes stay exactly right. The
distribution is not a guarantee, which is why the op still gets the
check; it is the reason the checkpoint does not have to be retrained.

## Copy-constructibility, leaking upward from `std::function`

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

## GIS

- **[want]** Ellipsoid (WGS84) as Space — geodesic distance on Earth

## Discrete / graph geometry

- **[want]** Graphs as metric spaces embedded in a Riemannian manifold — hyperbolic embeddings for hierarchical/tree-like data, using the existing `Hyperbolic<N>` space directly as the embedding target rather than a separate library.

## Geometry

- **[course]** Boolean ops on concave mesh (BSP tree)

## Native math (working defaults, not dependency removal)

Reframed 2026-09-13. The old heading was "dependency reduction" and the item below claimed native SVD/eigendecomposition was still open — both stale. Dense linear algebra landed natively on 2026-09-06 (`0372794` general-N symmetric eigendecomposition via cyclic Jacobi, `4c77d38` general MxN SVD via eigendecomposition of AᵀA, benchmarked in `cb7978a`), works for any `Scalar` including `Real<Digits>`, and is cross-validated against `Eigen::JacobiSVD` in `tests/test_svd.cpp`.

Corrected 2026-09-15: the sentence here used to add "which passes, contrary to the 'one pre-existing unrelated failure' note repeated in several Completed entries above". Both halves were true and the conclusion drawn from them was wrong. The test did pass, in Release; the earlier note was written from a Debug run where it aborted. The cause was ours, not Eigen's — the test asked `JacobiSVD` for thin U/V on a fixed-size matrix, which Eigen allows only when the column count is dynamic and otherwise asserts. `NDEBUG` removes that assert, and every CI job builds Release, so "it passes now" read as "the note went stale" rather than "the check stopped running". Fixed in #29 by asking for full factors, which the test's singular-value comparison does not care about either way. The general lesson is recorded as an open item below: a Release-only CI matrix cannot see this class of defect at all.

The goal was never to remove dependencies. It is to have a working default for every operation so no capability is gated behind an optional package, while any specialized library can still be absorbed behind the same interface — the pattern `SPATIUM_HAS_EIGEN` / `SPATIUM_HAS_BOOST_MULTIPRECISION` / `SPATIUM_HAS_IPC_TOOLKIT` already implements. Spatium does not need to beat specialized libraries at their own work; it needs to compose with them and still stand alone without them.

- **[course]** Sparse linear algebra — the real remaining gap, and the only place with no native default at all. The heat method, DEC, `differential.hpp` and `mesh/geodesic.hpp` all require `SPATIUM_EIGEN=ON` because there is no sparse matrix type or sparse factorization in the tree. Until that exists, geodesics via the heat method are a capability a default build simply does not have.
- **[want]** `Eigen::Ref<…>` public-API adapters — was blocked on dense SVD/QR, which has since landed, so this is now unblocked rather than waiting.

- **[want]** **Spatial acceleration on a manifold — the one part of "any space" that stops at Euclidean.** Noticed 2026-09-16 while measuring instancing, and it is a gap rather than a defect: nothing is currently wrong, because everything that would hit it routes around it. `examples/hyperbolic_tessellation_demo.cpp:8` says so in its own header — *"not a mesh + BVH::ray_cast() (that stack is Euclidean ray-primitive intersection)"* — and sphere-traces with `space.distance()` instead, linear in marker count, which is why that demo has dozens of markers and not thousands. `blackhole_gr_demo` has no BVH at all.

  An AABB cannot simply be carried over, and not only because it would be loose: a box needs coordinates, a chart does not preserve distances, and a sphere has no global chart without singularities — so a node that misses in the chart can be hit on the manifold. The intrinsic bound is the geodesic ball `B(p,r) = { q ∈ M : d(p,q) ≤ r }`: chart-independent, metric-consistent, and exact.

  **Half of this is thirty-five years old and should be borrowed rather than invented.** Ball trees (Omohundro 1989), vantage-point trees (Uhlmann 1991; Yianilos 1993) and M-trees (Ciaccia 1997) index metric spaces by balls for exactly this reason, with established build heuristics. What they answer is nearest-neighbour and range queries.

  **The open piece is the ray.** Casting a *geodesic* against a metric ball through such a tree is where the literature runs out; the closest relatives are sphere tracing in curved space (Hart / Hawksley / Matsumoto / Segerman, hyperbolic VR, ~2017 — which is what our own hyperbolic demo already does) and GR ray tracing, which brute-forces or adaptively steps rather than using a tree. The ball test has closed forms where the geodesics do: spherical and hyperbolic law of cosines. General manifolds need numerical geodesic integration, which `physics/relativity/geodesic.hpp` already has. SAH needs a volume measure, and on a manifold that is curvature-dependent — `2π(1 − cos r)` on the sphere, `2π(cosh r − 1)` in hyperbolic 2-space, an integral of the metric otherwise.

  Not an optimisation: without it there is no acceleration structure off flat space at all. And explicitly **not** an RSC problem — ROADMAP's own description above has RSC choosing *which existing* call to make, with Spatium as the exact substrate. It selects among methods that exist and cannot invent the ray-ball test; it gets work here once there are two ways to do this, not before.

  **Placed in the order, 2026-09-17, so that "not an RSC problem" does not turn into "not a problem yet".** `[want]` was the right grade for how settled the design is and the wrong signal about when to look at it: an item whose first line is *"the one part of any space that stops at Euclidean"* does not belong in the far drawer. The sequence is the HN post, then RSC, then this — and it comes directly after RSC rather than at some later unnamed point, because RSC is precisely the thing that will have two methods to choose between only once this exists. So the dependency runs the opposite way from what the paragraph above might suggest: RSC does not gate the ball tree; the ball tree is what eventually gives RSC something to select here.

## Materials & structures

- **[want]** Porous/microstructure analysis — geodesics on pore surfaces, effective conductivity/strength, and grain-shape/orientation analysis in metals, all as manifold queries on a reconstructed surface rather than a bespoke materials-science pipeline.

## GPU rendering (CUDA)

- **[course]** Land the actual production render (1920×1080, ~750 frames, Kerr flyby) — kernels are built and cross-validated (see Completed above); `gallery/blackhole_gr.mp4` currently ships a partial preview render, not the full sequence.
- **[want]** `gpu/derive_christoffel.py` already derives the closed-form Christoffel symbols symbolically (sympy) and self-checks them (Kerr at a=0 reduces to Schwarzschild term-by-term) — but only *prints* them for a human to hand-transcribe into `christoffel_closed_form.hpp`, instead of emitting the header directly. See `docs/gpu-abi-design.md` for the concrete fix (sympy's `cxxcode()` printer, write the file, no hand transcription step). That's what turns this from a one-off calculation into a standard, repeatable method.

- **[want]** **The C ABI for the DSL.** What gates it is how much of a scene is structural — not RSC, which an earlier draft of `docs/gpu-abi-design.md` named and which pointed at something unrelated, so nothing about it ever got closer. The document was rewritten 2026-09-18 and now carries the design; this entry records the state and the open question.

  Measured on every run of the donut demo, on `main`:

  ```
  fields: 107 total, 77 structural, 30 opaque; leaves 32
          (32 recognized, 0 unknown), 11 distinct types, 3208 B captured
  ```

  A `Trace` is a flat array of tagged ops addressed by index — the exact shape a C ABI wants, and `gpu/geodesic_kernel.cu`'s own boundary comment ("flat arrays only across this boundary — no structs/STL") says the same thing from the other side. So the *structure* crosses trivially. What never crosses is an opaque leaf, because it is a C++ closure.

  **The figure used to be quoted as 39 612 opaque leaves of 39 631, and later as 7 of 73. Neither is usable.** The first predates the dust collapsing from tens of thousands of per-particle nodes into a single `Scatter`. The second came from a report with two blind spots — it never walked a node's emission slot, and it never counted the closures a `Scale` or `Rotate` carries in `scale_fn`/`rot_fn`, which also left `is_structural()` returning true for a field with a lambda in it. Fixed 2026-09-18 with tests that fail without the fix; the numbers above are the honest ones.

  **The open question is not the count but the classification.** `unknown = 0` says every leaf carries a `type_index`, not that every leaf is a *known operation* — a lambda that happens to be a smoothstep and a lambda that is arbitrary are indistinguishable by that number. A first pass suggests all 32 fall into a small vocabulary (noise, hash-from-origin, smoothstep, trig, lerp, clamp, dot/cross, select, normalize, and an indexed gather from a small constant table) with nothing doing I/O or iterating a solver, but that pass was made by reading call sites rather than walking the trace, and the two disagree by one. Until each of the 32 is named with a vocabulary word, "exportable" is a guess.

  **The first step is a test, not a backend:** a POD interpreter run alongside `eval_into` over the whole donut scene, asserting bit-exact agreement. No epsilon — a tolerance passes under a different order of operations, and order of operations is exactly what a SIMT lane does not guarantee, so it would hide the class the test exists to catch. Passing demonstrates exportability for a real scene; failing names the missing operation. Neither outcome needs a GPU toolchain.

  **What the IR is *not* needed for, recorded 2026-09-18 because it kept being assumed.** Not for a browser: WASM is compiled C++, Emscripten handles `std::function` and lambdas, so an opaque leaf there is a compiled lambda that runs. Not for GPU *traversal*: `Cooked<T>` is already flat POD arrays with no pointers, so the CPU can cook and the device can traverse and shade with none of this work — the IR is needed only when fields must be re-evaluated on the device every frame, which is a real-time-animation requirement rather than a rendering one. And not for a language binding, which has a cheaper shape than an ABI: one entry point reading a *description*. `io/scene.hpp` is already half of that — JSON in, a string-keyed `ShapeRegistry`, and an unregistered `shape_kind` still loading and saving correctly rather than breaking. Extending it from placed objects to trace nodes is the job. Keep parsing and dispatch apart while doing it: a precise description determines its result, so it wants a deterministic parser that can fail cleanly, whereas RSC belongs one step further in deciding *how* — analytic or tessellated, what resolution, which root-finder, what precision, three of its existing trained domains. Parser decides what, RSC decides how. See `docs/gpu-abi-design.md`.

## RSC as search, not classification

Recorded 2026-09-15. Untested, but the pieces it needs already exist and were checked before writing this down.

**The reframing.** RSC today classifies: given a problem, pick the operation. The harder and more valuable question is *synthesis* — find the shortest correct **sequence** of operations. That is not open-ended code generation, it is combinatorial search over a closed, finite action space (a few dozen registered operations, every one already implemented and tested), which is the classical setting for program synthesis over a fixed DSL. No LLM belongs in that loop.

**Why the pieces fit.** `rsc/include/dispatcher.hpp` already emits a softmax distribution over the registry — that is a policy network in the AlphaZero sense, already trained. What is missing is a value head (a second small output estimating how close a partial chain is to the goal, so a tree search prunes instead of enumerating) and the search itself, wrapping the existing `Registry`/`Chain`. Reward shape is the one `rsc/README.md` already cites AlphaDev for: correctness is a hard filter via the existing `grade()`, length is minimised only among candidates that pass it.

**The cheap first experiment, using only what exists.** `rsc/include/chain.hpp` is teacher-forced on add/multiply chains, with a dispatcher measured at 0.312 → 0.972 on 3-step chains (`rsc/README.md:241`), and `tests/test_rsc_chain.cpp` already in place. Replace teacher forcing with free-running search that uses the trained dispatcher as its heuristic, and see whether it finds chains as short as, or shorter than, the teacher-forced path. `chain.hpp:21` records that multi-step credit assignment was never attempted precisely because teacher forcing sidesteps it — which is the gap this closes. Three existing pieces connected, nothing built from scratch; if it fails, the failure points at a specific next thing to diagnose, which is how the rest of this project is already run.

**The part nobody else can do, and the reason it is worth trying.** A chain found this way is expressed over *concepts*, not concrete types. An algorithm discovered by searching on `Sphere` transfers to `Hyperbolic` and to an arbitrary `ParametricSurface` by construction, because it was never written against a particular space. AlphaTensor found faster matrix multiplication, but that discovery is bound to a specific ring and does not generalise itself to another algebraic structure. Space-agnostic discoveries are not something anyone has searched for, because nobody else has a single interface through which such a generalisation could even be *checked*. That is the claim worth testing, and it costs one experiment to start.

**Measured 2026-09-18, and it changes where the search has to happen.** A first attempt built the search over a toy numeric action space -- five integer operations, reach a target -- and the result was a clean negative: a trained policy used to order a best-first search was worse than uninformed breadth-first on every axis at once, solving less often (90.8% against 100% at depth 6), expanding 1.6-1.7x *more* nodes, and returning a shortest chain half as often, while costing 600 ns per expansion on top. See `rsc/tools/search_heuristic.cpp`.

The reason is structural rather than a training failure, and it generalises: breadth-first is the provably optimal order for shortest paths under uniform cost, so a heuristic there competes with an optimal ordering rather than improving a poor one. It has to carry real information merely to break even.

But the deciding fact is about the **state type**, not the method. Numbers collide: `10+3+3` and `2*8` reach the same state, so the search space collapses. Measured on that toy, 9 765 625 paths at depth 10 reduce to 13 836 distinct states -- a 706x collapse, and an effective branching factor of 2.13 rather than 5. A space that small is searched exhaustively for nothing, and no value head changes that.

`Registry::chainable_ops()` selects on arity -- `in_size==2 && out_size==1`, so the output feeds back as the next accumulator -- which for the current registry means `add`/`multiply` over a *scalar* accumulator. That was always labelled a mechanism check, and the collapse is a property of scalars, not of chains.

**Geometries do not collide.** Two different sequences of mesh operations essentially never produce the identical mesh, so generalising the accumulator from a scalar to an object of the substrate leaves the branching factor intact and puts the search in the regime this entry is about. The selection rule does not change at all; only the type of the state does. (That non-collision is reasoning from the structure, not yet measured -- hashing meshes after distinct sequences and counting collisions is the cheap check, the same one that produced the 706x figure for integers.)

It also inverts the economics, which the toy got backwards. Expanding a node over numbers is five integer operations and a map probe -- tens of nanoseconds -- so a 600 ns policy evaluation costs more than ten nodes and can never pay. Expanding a node over geometry means *running* a mesh operation: milliseconds. Against that a learned heuristic is free, and the distilled 8 ns tree (see `rsc/tools/distill_tree.cpp`) is free several times over.

**What the search is actually for, stated so the ambition is calibrated.** It finds *compositions* of registered operations; it cannot invent an operation. "A new algorithm" here means a sequence nobody wrote down -- the same sense in which AlphaDev found shorter sorting sequences out of existing instructions. What makes that more than a curiosity is the concept-generic part: a chain of concept-constrained ops is a program over the concept, so running it on `Sphere` and on `Hyperbolic` is two instantiations of one algorithm rather than two algorithms.

**And the existing design already predicts which half transfers.** This file's own dispatch/calibration split -- dispatch is the trained model, calibration is classical optimization and not the model -- maps exactly onto it: the *dispatch* transfers, since which operations in which order is space-agnostic by construction, while the *calibration* does not, because a step size tuned on a unit sphere will not hold on a hyperbolic space. Structure travels, constants are re-fitted per space.

**The experiment that would settle it:** search for a chain meeting a specification on one space, instantiate the same chain on two others, and measure the transfer rate -- and separately whether the transferred chain is any worse than one found by searching the target space directly. That is the claim, run on its own machinery.

Related and already in this file: manifold-native RL (Object model section), Fisher-Rao as a `RiemannianManifold` (Manifold applications) — the closest to buildable, since it needs no new abstraction — stochastic processes on manifolds, and hyperbolic embeddings for hierarchical data (Discrete / graph geometry).

## Two searchable domains, and why owning the library is what makes them searchable

Recorded 2026-09-18, pulling together threads that were already here separately.

**The reason everything is hand-rolled turns out to be a training reason, not only an independence one.** A search needs to grade its candidates, and grading at the resolution a search needs means running instrumented code you own. On someone else's library an operation is a black box: it can be called, not scored. Here every registered operation is simultaneously a legal action and a measurable outcome, which is the whole substrate. `Vec`, the solvers, the eigendecomposition, `json.hpp`, `wav.hpp` — the argument for those was self-containment, and this is the second one.

There are exactly two places where an oracle is free, exact, and ours.

**Rewriting the IR.** "The trace is the IR" (see the Declarative scene DSL section). A rewritten op graph must compute the identical result, and the oracle for that is already queued rather than new work: the POD interpreter asserting bit-exact agreement with `eval_into`. The test that would demonstrate exportability *is* the grader for a rewrite — one artifact, two uses. This is also, precisely, "arrive at one algorithm's results from a different set of primitives": a rewritten graph computes the same thing by a different composition, which is the shape AlphaTensor's result had.

**Building trees.** The oracle is measured query cost — build the hierarchy, then run queries and time them.

The difference between the two is worth keeping straight, because it decides what each search can be asked for. **IR search can only equal**: bit-exact equality admits no improvement on correctness, so the objective is cost subject to exact identity. **Tree search can beat**: there is no output to match, only a cost to lower, so improvement over the reference is on the table.

**Domain-specific models are not a hypothesis here.** Measured 2026-09-18: one distilled tree per domain loses 0.005 accuracy against the network at 8.43 ns, where a single tree covering all four domains loses 0.008 at 11.18 ns — better and roughly half the cost, because the domain is known at the call site and should not be costing the tree two levels to re-derive. The base-plus-custom split already in the RSC design is the same idea one storey up.

**And the manifold tree is where the transfer claim gets tested cheaply.** A BVH and a ball tree differ only in the bound; construction and traversal are the same. Given a bound concept — bound a set, merge two bounds, test a query against a bound — one hierarchy serves axis-aligned boxes, geodesic balls, OBBs and k-DOPs alike. A split policy learned over that concept on Euclidean boxes can then be run against geodesic balls, and whether it still helps is a measurement, not an argument. That is the space-agnostic discovery claim in its cheapest concrete form, and the compiler does the checking.

**Dependencies, so the order is not confused.** IR search needs the op vocabulary and the POD interpreter, both already queued. Tree search needs the bound concept, which does not exist yet. Neither needs the other, so either can be taken first.

## Contact physics / RSC

- **[course]** Full multi-config sweep for the implicit-contact Newton solver, matching the XPBD investigation's own 20-point discipline, plus a performance pass (current single traced config takes multiple CPU-minutes under ipc-toolkit's TBB-based collision detection).
- **[course]** RSC's calibration-search against the implicit-contact pipeline — not yet built.
- **[dare]** RSC "real-time control of complex dynamics" domain (illustrative: an underwater drone) — depends on the base+custom deployment split actually working.

## Interop / ecosystem

- **[dare]** Heat-method log map, CGAL-grade exact polyhedral geodesics (geometry-central and CGAL each cover one half of this; Spatium currently ships neither on top of Dijkstra/heat-distance).
- **[want]** Own Vec/Matrix creates impedance mismatch with the Eigen ecosystem — interop adapters beyond the current `to_eigen`/`from_eigen`/`eigen_view`. Previously gated on native SVD/eigendecomposition, which landed 2026-09-06 (see Backlog → Native math), so nothing blocks this now.

---

# Competitive Position

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
