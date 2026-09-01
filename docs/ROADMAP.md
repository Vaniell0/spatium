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
- **Table renderer fix** (`src/io/table.cpp`) — separator width off-by-one (`w+1` → `w+2`); console tables now render with aligned borders
- **Eigen integration decision** — remain optional via `SPATIUM_EIGEN=ON`; heat method + differential + interop stay isolated. Vec/Matrix not migrated to Eigen; `Eigen::Ref<…>` planned for public APIs only once dense QR/SVD lands (see Backlog → Native math)
- ~~CPU raytracer: Quadric → viewer texture (analytical render without mesh)~~ done: `ray_quadric`/`ray_quadric_proximity` in `geometry/ray_surface.hpp`, exercised by `examples/parametric_analytical_demo.cpp`'s `glow_sphere.png`; ships PNG output rather than a live viewer texture, which covers the same "analytical render without mesh" goal for the gallery use case
- ~~Ray-ParametricSurface: Newton UV iteration for arbitrary f(u,v) surfaces~~ done: `geometry/ray_parametric.hpp` — no tessellation, no BVH, Newton-solves `S(u,v) = o + t·d`; exercised by `examples/parametric_analytical_demo.cpp` producing `klein_analytical.png`/`mobius_analytical.png`/`bumpy_analytical.png`/`parametric_gallery.png`

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

Docs truing-up pass: `architecture.md` concept-hierarchy fix + a new "header-only spine, and three principled exceptions" section, `CHANGELOG.md` deleted (redundant with this file), the "Block A/B/C/D" physics-milestone naming stripped from source comments and — completed in a follow-up pass — from `ROADMAP.md`/`api-reference.md`/`concept-driven-physics.md` as well, `cloth_sphere_probe.cpp` removed (superseded investigation tool, predates several since-changed APIs). `Result<T>`/camera-module dedup deferred.

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

## C++23 modules: caught up to a working state (2026-08-28)

The 2026-04-24 migration (see above) had gone stale — this closes the gap rather than deleting the subsystem, since the missing partitions were genuinely mechanical, not a design problem.

- **New algebra partitions**: `dual`, `ode`, `linear_solve`, `calculus` — added to `modules/algebra/`, registered in `algebra.cppm`'s export list and `CMakeLists.txt`'s `PARTITIONS`. `eigen_interop.hpp` stays deliberately header-only (SSE intrinsics from `Eigen/Core` conflict with `:vec_simd`'s BMI — this was already documented, not new).
- **`physics/relativity/*` folded into `spatium.physics`**: `schwarzschild.hpp`, `kerr.hpp`, `geodesic.hpp`, `accretion_disk.hpp` added to `physics.cppm`'s single-TU `#include` block (no cross-includes among the four, any order works; `geodesic.hpp`'s `Dual<T>`/`ode`/`linear_solve` needs are covered by the algebra partitions above). `physics/atomic/`'s own stale include paths (left over from its 2026-08-28 split out of `physics/`'s flat top level) were also fixed in passing — a real, separate bug: `modules/physics.cppm` still referenced the pre-split `physics/orbital.hpp`-style paths, breaking the modules build specifically since the header-tree build never exercises that file. Drive-by cleanup: the "Block A/B/C/D" milestone-naming comments the 2026-08-26 architecture audit was supposed to have stripped everywhere had survived inside this one file; removed.
- **`spatium.render`: new module**, `modules/render.cppm` — single partition covering `camera.hpp`, `color.hpp`, `spectral.hpp`, `parallel_for_rows.hpp`, `supersample.hpp`, `sky.hpp`. `write_image.hpp` deliberately stays header-only: its whole contract is "exactly one translation unit defines `STB_IMAGE_WRITE_IMPLEMENTATION` before including it," which is a link-time property, not a module one — folding it into the module would bake `stb_image_write`'s externally-linked function bodies into `spatium_render_module` itself, and any TU that also defines the macro (every example does, so does `tests/test_render.cpp`) and links against the module in the same binary would hit duplicate-symbol errors. Same category of documented exclusion as `eigen_interop.hpp` and `mesh/dec.hpp`/`mesh/primitives.hpp` — a real technical call, not neglect.
- **Real GCC modules-ts bug found and fixed, not routed around**: linking a consumer of `spatium.render` against `Sky`/`make_starfield()` failed with `undefined reference` to `std::vector<int, allocator<int>>`'s special members (`_Vector_impl`'s destructor, the default constructor, the move constructor) — reproduced from a clean isolation (removing `Sky` usage from the test made the link succeed; adding real `#include <vector>` in both the producing module and the consuming TU did not fix it). Root cause narrowed to `Sky::star_buckets`'s element type specifically: switching `std::vector<std::vector<int>>` to `std::vector<std::vector<std::uint32_t>>` (matching the index-type convention already used elsewhere in the codebase — `BVH::Hit::index`, mesh vertex/face ids — not a workaround invented for this bug) made it link cleanly in both the modules and header-tree builds. A genuine, narrow GCC/libstdc++-modules-ts limitation with plain `int` crossing a module boundary through a nested container inside an inlined function, not an application-level mistake.
- **CI**: `.github/workflows/ci.yml` gained a `modules-build` job (configure/build/test via the existing `modules` CMake preset) — the gap ROADMAP had flagged as the actual reason this was allowed to rot the first time.
- **Verified**: 825/825 tests green on the `modules` preset (`build-modules/`, fresh configure) and 825/825 on the existing `build/` tree (which already had `SPATIUM_USE_MODULES=ON` cached from earlier work). `scripts/gen_dependency_graph.py --check` and `scripts/check_claude_md_layout.py` both still pass (neither script inspects `modules/`, so unaffected either way, checked to be sure).

---

# Backlog

Grouped by topic, not by version — an item sits here until it's ready to become real work, then moves to Completed above with the date it landed.

## Mesh processing

- Mesh simplification (edge collapse + QEM)
- Isotropic remeshing on Surface
- PLY import/export

## Analytical rendering

- Fragment shader raymarcher: GPU-native analytical render (GLSL quadric math)
- Ray-quartic: torus intersection via solve_quartic
- Dual Quaternion

## Manifold applications

- Fiber bundles (tangent/cotangent)
- Geodesic FEM (Laplace-Beltrami, heat equation)

## GIS

- Ellipsoid (WGS84) as Space — geodesic distance on Earth

## Geometry

- Boolean ops on concave mesh (BSP tree)

## Native math (dependency reduction)

- SVD / eigendecomposition — a real native-implementation candidate (unlike the heat method's sparse-Cholesky step or ipc-toolkit, whose cost/benefit doesn't favor a from-scratch rewrite). Blocks `Eigen::Ref<…>` public-API adapters noted above and dense-solve paths that currently require `SPATIUM_EIGEN=ON`.

## GPU rendering (CUDA)

- Land the actual production render (1920×1080, ~750 frames, Kerr flyby) — kernels are built and cross-validated (see Completed above), but the render's completion status was last confirmed "in flight" on 2026-08-28; treat as open until confirmed landed.
- `gpu/derive_christoffel.py` already derives the closed-form Christoffel symbols symbolically (sympy) and self-checks them (Kerr at a=0 reduces to Schwarzschild term-by-term) — but only *prints* them for a human to hand-transcribe into `christoffel_closed_form.hpp`, instead of emitting the header directly. See `docs/gpu-abi-design.md` for the concrete fix (sympy's `cxxcode()` printer, write the file, no hand transcription step). That's what turns this from a one-off calculation into a standard, repeatable method.
- `gpu/` itself is not tracked in git yet (part of the uncommitted pile) and won't be until the above makes it a standard mechanism rather than a one-off port.

## Contact physics / RSC

- Full multi-config sweep for the implicit-contact Newton solver, matching the XPBD investigation's own 20-point discipline, plus a performance pass (current single traced config takes multiple CPU-minutes under ipc-toolkit's TBB-based collision detection).
- RSC's calibration-search against the implicit-contact pipeline — not yet built.
- RSC "real-time control of complex dynamics" domain (illustrative: an underwater drone) — depends on the base+custom deployment split actually working.

## Interop / ecosystem

- Heat-method log map, CGAL-grade exact polyhedral geodesics (geometry-central and CGAL each cover one half of this; Spatium currently ships neither on top of Dijkstra/heat-distance).
- Own Vec/Matrix creates impedance mismatch with the Eigen ecosystem — interop adapters beyond the current `to_eigen`/`from_eigen`/`eigen_view` are planned once SVD/eigendecomposition (above) lands.

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
