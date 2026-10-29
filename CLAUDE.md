# Spatium

C++23 header-only math library for arbitrary mathematical spaces, geometric primitives, and mesh operations.

## Build

```bash
nix develop
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug -DSPATIUM_BOOST=ON
ninja -C build
ctest --test-dir build
./build/examples/showcase
```

`-DSPATIUM_BOOST=ON` is what keeps the full test set in play: without it
`test_precision.cpp` and the symmetry registry are dropped. Leave it off
only to reproduce what someone with no system packages gets --
`cmake --preset noboost`, the counterpart to `noeigen`, does exactly that
and is the configuration a bare `cmake -B build` with no flags now
produces.

A second config, `build-release/` (`-DCMAKE_BUILD_TYPE=Release -DSPATIUM_EIGEN=ON`,
modules off), exists alongside `build/` for heavy runs. Configure/build it the
same way, pointed at `build-release` instead of `build`. Measured no real
wall-clock difference vs. Debug on the actual project build (a standalone
`-O0` vs `-O2` comparison suggested a ~5.7x win that didn't reproduce here,
cause not tracked down, possibly modules-related). Keep the config anyway;
don't assume it's a speed fix without re-measuring.

## Architecture

### Spaces (concept hierarchy)

```
Set → TopologicalSpace → MetricSpace → NormedSpace → InnerProductSpace
                    ↓                                        ↓
              Manifold → RiemannianManifold          Complete + Inner = Hilbert
                    ↓ ↘ LorentzianManifold (interval, causal; no distance)
              Surface (+ project/normal)               Complete + Normed = Banach
```

Concrete: Euclidean<N>, Sphere<N>, Hyperbolic<N>, Minkowski<N> (Lorentzian)

### Geometry (primitives + operations)

```
Shape concepts: Shape, ClosedShape, Measurable, Bounded, DistanceQueryable, BoundedRegion
Primitives:     Line, Ray, Segment, Hyperplane, Triangle, Polygon, Circle, Disk, Box, Simplex, Quadric
Operations:     intersect() free functions, distance() free functions (all shape pairs), operator| pipe syntax
                clip() — constrain results to shape bounds (point/line/segment → shape)
                intersect_via_subspace() — generic intersect for BoundedRegion pairs
                intersection_region() — coplanar boolean (Sutherland-Hodgman)
                difference_region(), difference_area(), symmetric_difference_area()
                operator& / operator- / operator+ on Polygon (intersection/difference/union)
                union_of() — named polygon union alternative
                ray_quadric() — analytical ray-quadric hit, ray_quadric_proximity() — miss via complex roots
                lazy() — deferred transform composition, Result<T> pipe-unwrap through Transform/Morphism
Subspaces:      Triangle→Plane, Segment→Line, Disk→Plane, Circle→Plane, Polygon→Plane, Ray→Line
Factories:      tri(), seg(), box(), line(), ray(), plane(), circle(), disk(), poly()
                Quadric::sphere(), ::cylinder_z(), ::cone_z(), ::ellipsoid()
```

### Mesh

```
Mesh<Surface> — indexed triangle mesh for any Surface
subdivide_once / subdivide — midpoint subdivision with surface projection
LodChain — multi-level LOD
icosahedron / tetrahedron — sphere starting meshes
GeodesicMethod::Dijkstra | Heat — method selection enum
HeatSolver<S> — pre-factored heat method (Crane 2013), O(h^2), requires Eigen
differential.hpp — cotangent Laplacian, mass matrix, face gradients, divergence
```

### Extras

- `Point<Space>` — type-safe point wrapper
- `Morphism<From, To>` — runtime maps with pipe composition: `point | f | g`
- `Complex<T>` — complex arithmetic, sqrt, cbrt, constexpr, std::format
- `solve_quadratic / solve_cubic / solve_quartic` — polynomial solvers → Complex roots
- `algebra::power / commutator / adjoint / poly_eval` — generic functions (Group/LieGroup/Ring concepts)
- `Dual<T>` — forward-mode autodiff, satisfies Scalar (drops into existing ADL-style code unchanged)
- `Function<F,Domain,Codomain>` concept, `gradient() / integrate() / minimize()` — calculus over plain callables, no wrapper type
- `raise_gradient() / project_tangent() / riemannian_minimize()` — Riemannian gradient descent on any RiemannianManifold+Surface (Euclidean, Sphere, Hyperbolic): index-raise the ambient covector via the space's own metric_at(), project to the tangent space, retract via exp_map
- `verify_metric / verify_inner_product / verify_exp_log` — axiom verification
- `Real50 / Real100` — arbitrary precision (Boost.Multiprecision, requires `SPATIUM_BOOST=ON`)
- `Table / Svg` — structured output and 2D visualization
- `std::format` support for all types
- UDLs: `_deg`, `_pi`, `_x`, `_y`, `_z`

### CMake targets

- `Spatium::sdk` (INTERFACE) — header-only, link this
- `Spatium::core` (INTERFACE, will become SHARED with .cpp sources)
- `spatium` — backward-compat alias for sdk

## Directory layout

A map, not a manual: one line per directory. What each header is for and what it declares is generated from the headers themselves into [`docs/capabilities/`](docs/capabilities/README.md) (one page per domain, `scripts/gen_capabilities.py`, checked in CI) -- search there first.

- `include/spatium/core/` — concepts (`Set` ... `RiemannianManifold`, `LorentzianManifold` as a different concept with no distance), `access.hpp` (the customization points in `spatium::spaces`: `exp_map`/`log_map`/`metric`/`distance`/`project`/`normal` found as a member, then by ADL, then derived; `geodesic`, `midpoint`; **new generic algorithms call these, not members** -- see `docs/architecture.md`), `Dimension` (finite, dynamic or infinite), `Causal`, `Result<T>`, the axiom verifiers (`verify_metric`, `verify_exp_log`, `verify_lorentzian`), `Real50`/`Real100`, `UpTo<T, N>`
- `include/spatium/algebra/` — `Vec`, `Matrix`, `Quaternion`, `Octonion`, `Complex`, `Dual` (autodiff), `Series` (limits and infinity as arithmetic: `limit`, `AtInfinity`), calculus (`gradient`, `minimize`, Riemannian descent and Frechet mean), quadrature (**`integrate(f, domain, options)` is the one door for the integral**: the domain is a type -- `Finite`, `HalfLine`, `WholeLine<T>` -- `checked(domain)` proves the ends first, `IntegralOptions{.method, .witness}`; the answer is an `IntegralResult`), `monte_carlo`, ODE solvers (`integrate_extrapolated`), linear solve, polynomial solvers. `inline namespace algebra` -- see `docs/conventions.md`
- `include/spatium/algebra/groups/` — `SO3`, `SE3`, non-inline `spatium::algebra::` (see `docs/conventions.md` for why these stay qualified-only); `HyperoctahedralGroup<N>` (signed permutations, acting on integer vectors and on the cells of a cube) and `action.hpp` (`GroupAction`, `FiniteGroup`, `orbit`, `stabilizer`, `orbits`, `verify_action`)
- `include/spatium/algebra/models/` — models as functions with their gradients, templated on the scalar: `LinearSoftmax<K, T>` (a softmax over a linear score; the distribution into the caller's buffer, sample, choose, cross-entropy and gradient at a temperature -- distillation from a teacher's soft distribution -- with no allocation; a `Chooser`), `TanhNet<F, H, T>` (one hidden tanh layer to a number in (-1, 1); squared loss and its backward pass, held against `Dual`), `ScoreTree<K, T>` (a binary tree of column thresholds scoring a row, with its description length in bits), `softmax` of scores
- `include/spatium/algebra/analysis/` — step rules that move parameters given a gradient: `AdaGrad<N, T>`, `Adam<N, T>`
- `include/spatium/statistics/` — fitting parameters to a sample: `fit_pass` (one pass in a drawn order, the loss and gradient from the caller's objective, the step from a rule); how much a distribution says and how far two are apart: `entropy`, `cross_entropy`, `kl_divergence`, `in_bits` (`information.hpp`); growing a regression tree by CART (`tree_fit.hpp`: `fit_tree`); `spatium::statistics::`
- `include/spatium/spaces/` — `Euclidean`, `Sphere`, `Hyperbolic`, `ConstantCurvature<N>{kappa}` (all three as one family), `Minkowski<N>` (Lorentzian, c = 1), `L2Interval` (infinite-dimensional), `LieGroupManifold<G>`, `ProductSpace`, typed `ParametricSurface`/`ImplicitSurface` (exact partials through `Dual`; `.erased()` is the boundary), `MetricChart` (a space from its metric alone: exp/log derived), `pullback_metric`, `integrate_volume`/`volume_of`, `cos_sinc.hpp`. `chart_of(space)` is the ADL extension point into the scene DSL. Includes nothing of `mesh/`
- `include/spatium/geometry/` — primitives (`Line`, `Ray`, `Segment`, `Triangle`, `Polygon`, `Box`, `Simplex`, `Quadric`, ...), `intersect`/`distance`/`clip` free functions, polygon booleans, lazy transform chains, `ray_surface`
- `include/spatium/mesh/` — `Mesh<Surface>`, subdivision (at the space's own midpoint), LOD, topology, geodesics (Dijkstra, heat method), transport (Schild's ladder), Voronoi, conform, scatter, the materialisers `tessellate` and `marching_cubes` (a surface becomes triangles only in `mesh/`, so `spaces/` includes nothing of it)
- `include/spatium/spatial/` — BVH over any `Bound`; `GeodesicBall`/`GeodesicBallTree` for manifolds with no coordinates a box could be drawn in
- `include/spatium/render/` — the CPU raytracer engine (camera, supersampling, PNG, sky, texture) and the GPU scene layouts and shaders (`cooked_scene`, `gpu_scene`, `lbvh`, `blackhole_glsl`, `gpu_instances`)
- `include/spatium/discrete/` — `FiniteSet`, `GeometricSet`, overflow-checked combinatorics, and the rules of chess as an oracle that reads no text: the 8x8 board in one word (`bitboard.hpp`: `Square`, `Bitboard`, the board's D4 as signed 2x2 matrices), the squares a piece reaches (`attacks.hpp`), the position and its rules (`position.hpp`: `PositionSetup` checked by `Position::from_setup`, `legal_moves`, `perft`, `outcome` for mate, stalemate, fifty moves, dead positions and repetition), and a board of any number of axes with pieces as orbits of the cube's symmetries (`lattice.hpp`), and a move as a row of numbers (`move_features.hpp`: `steps_taken`, `geometric_features`, `join`)
- `include/spatium/physics/` — `Element` (the one compiled, non-header-only translation unit)
- `include/spatium/physics/atomic/` — visualization-support atom models
- `include/spatium/physics/mechanics/` — compile-time SI units, bodies, forces, integrators (Euler/RK4/Lie-group/LGVI), symplectic and variational integrators, continuum mechanics, IPC barrier + XPBD, **certified continuous collision** (`narrow_phase.hpp` `distance_bound`, `rigid_contact.hpp`, `surface_ccd.hpp`), wave PDEs. Does not yet use `Result<T>` for fallible ops (see `docs/architecture.md`)
- `include/spatium/physics/relativity/` — metric-agnostic geodesic integration: `Schwarzschild`, `Kerr` (templated callables, `Dual`-substitutable), `MetricField` (a metric as ten IR fields), GLSL emitters, `SpacetimeScene`, the accretion disk
- `include/spatium/io/` — the scene field IR (`field.hpp`, lowered to plain data by `field_pod.hpp` and to GLSL by `field_glsl.hpp`), `build.hpp` (`spatium::io::build`, the declarative scene DSL: a `Trace` of Space/Offset/Scatter/Compose operations, analytic end to end; see `docs/getting-started-dsl.md`), JSON, scene, OBJ/STL/SVG/WAV writers, `Table`, chess notation (`chess.hpp`: FEN, SAN, `e2e4`; a FEN is read into a `PositionSetup` and checked by the rules)
- `include/spatium/viewer/` — the Vulkan app, headless compute (`Context`, `Kernel`), `DeviceLbvh`
- `include/spatium/` — `Point`, `Morphism`, umbrella headers
- `include/spatium/vendor/` — stb_image_write, stb_image (public domain; a caller defines the `*_IMPLEMENTATION` macro in one translation unit)
- `tests/` — Catch2; `ctest -N` reports the count, which depends on which optional dependencies are enabled. `tests/connectivity/` is the connectivity matrix (scalars x spaces x generic algorithms, graded L0 compiles / L1 finite / L2 axiom or symmetry / L3 agrees across scalars; `expected.hpp` generated, held by `test_connectivity.cpp`)
- `examples/` — see `examples/CMakeLists.txt` and each file's header; `donut_demo` is the DSL getting-started (`docs/getting-started-dsl.md`), `blackhole_live` the black-hole renderer (`--check`, `--bench`, `--live`, `--frame`, `--video`), `napkins` and `cloth_live` the cloth scenes, `chess_agent` an agent that learns to choose among legal moves from a player it watches or from won games (`--observe`, `--play`, `--reinforce`, `--distill` for the agent distilled into trees of growing depth beside its own description length; `--observe 64` is the ctest `chess_agent_observes`), `chess_zero` the tree search of `rsc/core` learning chess by playing itself, measured against a random, a greedy and a two-ply material player (`--iterations`, `--games`, `--sims`, `--threads`), `scripts/external/ipc_compare/` the out-of-tree ipc-toolkit cloth-on-sphere comparison (offline/Vulkan raytracers and diagrams: `spatium_add_example()` in `cmake/SpatiumTarget.cmake`)
- `benchmarks/` — Google Benchmark (vec, intersection, mesh, bvh, orbital, raycast)
- `rsc/` — the synthesis engine, being rebuilt: `rsc/core` is generic (`graph.hpp`: a typed graph and the query for the routes to a quantity, with its own JSON reader; `tree_search.hpp`: a tree search guided by a learned prior and value over any `SearchGame`, through an `Evaluator`, and `self_play`; it includes nothing of `<spatium/...>`, which a CI grep holds), `rsc/spatium` binds it to this library (`synthesis.hpp`: particular analytic solutions synthesised as IR fields on an interval -- expressions over the interval's variable built bottom up, those with equal values merged, a sum of a few fitted by least squares and certified on a dense grid -- and piecewise where one does not reach the tolerance; `query.hpp`: a question with a free symbol -- `solve` a root on a bracket, `integrate` over an interval -- written as a field of U and V, its meaning computed by the library and its answer, a field of the parameter, synthesised and certified, refused where the question has no value; `gen_graph.py` writes `graph.json`, the library as a graph over its concept hierarchy -- concepts parsed, type-satisfies-concept compiled, operations from `access.hpp` with the connectivity matrix's levels -- `--check` in CI; `rsc_route <type> <operation>` prints the routes, built by name; the graph has edges between types -- a space over `Real50`/`Dual`, a product from its factors, a Lie group manifold from its group). `rsc_find.py` finds and certifies: each route becomes C++ (`emit.py`, one function per edge), is compiled, run against a reference on the matrix's samplers, and its outcome recorded in `find-report.md`; routes that agree and are independent of the reference become pairs in `tests/symmetry/pairs/rsc_*.cpp` (generated, `--all --check` in CI). The 2026-09 prototype is at tag `rsc-prototype-final` (`docs/roadmap/2026-10-20-01-the-rsc-prototype-leaves-the-tree.md`)
- `cmake/` — `SpatiumConfig.cmake.in`, `SpatiumModule.cmake` (C++23 modules helper), `SpatiumTarget.cmake` (`spatium_add_example()` helper)
- `scripts/` — `gen_connectivity.py` (builds every connectivity cell alone, writes `docs/connectivity.md` and `tests/connectivity/expected.hpp`, `--check` in CI); `gen_capabilities.py` (writes `docs/capabilities/`, one page per domain, so changes in different domains never edit the same lines; `--check` in CI); `gen_consumers.py` (the consumer rule, counted: which declared names nothing else mentions, by identifier; `docs/consumers-baseline.json` holds today's unused names, `--check` in CI fails on one added); `gen_roadmap.py` (a roadmap history entry is one file `docs/roadmap/<date>-<nn>-<slug>.md` with a front matter -- date, title, optional stage, `closes` (findings rows deleted), `files` (must exist) -- not rendered into `docs/ROADMAP.md`, so two changes in flight never edit the same lines; an entry's `open:` list says what it did NOT do and a later entry's `resolves:` closes items, so `--open` prints what is still not done, computed; `--check` validates, no argument prints them as one page); the doc-freshness checks, which run in CI so docs cannot silently drift: `gen_dependency_graph.py --check` (real `#include` graph vs the committed `.dot`), `check_claude_md_layout.py` (every `include/spatium/` subdirectory mentioned above), `check_doc_file_refs.py` (every source file the docs name exists; names with no counterpart -- recorded deletions, proposals, generated files, other projects' -- are listed with a reason in the script's `NOT_IN_TREE`)

## Principles (read `docs/architecture.md` before new work)

- **Interfaces, relations, capabilities -- not implementations.** Algorithms are written against concepts and customization points (`core/access.hpp`: member, then ADL, then derivation); descriptions stay lazy (expression templates, the Field IR, callables kept as types -- `std::function` only as an explicit erased form at a boundary); a scene is a program of programs that ends in a picture.
- **Physics runs on spaces and charts, not on meshes.** Contact, distance and continuous collision are computed on the analytic or parametric description through certified bounds (`narrow_phase.hpp`, `surface_ccd.hpp`). A mesh is a materialiser at the edge -- for drawing, and as the baseline a method is compared against -- not the representation physics is built on. This is the limitation IPC has and we do not: its barrier and CCD are defined on a distance but computed only between triangles. Deformables should be smooth patches (control points as degrees of freedom), not triangle soups; `napkin_scene.hpp` is still triangles and is the comparison, not the design.
- **No module without a consumer.** A module with no real caller that exercises it counts as no code and no tests, however good the idea: it is written consumer first, and the consumer and the test that goes through it land together. A result nothing calls is a note (findings, ROADMAP), not something in the tree. It applies to every directory, `rsc/` included.
- **Consistency comes from automatic checks, and the cheapest check is a symmetry**: two paths to one answer that must agree (member / ADL / derived, `double` / `Real50`, `Dual` / finite differences, host / device, erased / typed, C++ / GLSL). Prefer such a pair over a single-path assert; every verifier gets a case it must fail. More reuse through the interfaces means more paths that meet, so more free checks. See `docs/findings-2026-09.md` sections 0 and 11.

## Working on the repo

- **Generated files and their order.** `docs/connectivity.md` and `tests/connectivity/expected.hpp` come from `scripts/gen_connectivity.py`; `rsc/spatium/graph.json` from `rsc/spatium/gen_graph.py`; `rsc/spatium/find-report.md` and the witness pairs from `rsc/spatium/rsc_find.py`. Each reads the one before it, so regenerate in this order: matrix, graph, find. CI regenerates all three and compares; running them out of order is the usual cause of a red "up to date" job. `gen_capabilities.py`, `gen_dependency_graph.py`, `gen_consumers.py` and `check_claude_md_layout.py` have `--check` modes that CI runs too.
- **Contract of `distance`.** `spaces::distance(s, p, q)` is the geodesic distance `|log_p q|_g`, or it does not exist: erased surfaces (`ParametricSurface<double>`, `ImplicitSurface<double>`) have only `chord_distance`, and Lorentzian spaces (`Minkowski`) have none. Code that needs a length on such a space says `chord_distance` or takes the metric explicitly; never derive a distance from a retraction.
- **Branches and PRs.** Code goes through a PR from a worktree (`git worktree add ../spatium-wt/<branch> origin/main`); the CI jobs in `.github/workflows/ci.yml` must all be green, including the Nix-free gcc15 job that installs its compiler from a PPA (a PPA outage turns five jobs red without any code cause; rerun later). Documentation without code changes is committed straight to `main`.
- **Words in public text.** Code comments, commits, PR titles and bodies, `docs/roadmap/` entries and docs never name internal plan items, phases or task labels. They say what the code does and why, with the measured numbers.
- **Builds.** `nix develop`, at most `-j6`; look at the exit code of a build before reading its output (a pipe into `tail` hides a failure). The Real50 cells of the matrix are heavy: a full `[connectivity]` run belongs to Release, and each cell has a 600 s timeout.

## Conventions

- Namespace: `spatium::`, `spatium::algebra::` (inline — see `docs/conventions.md`), `spatium::geometry::`, `spatium::mesh::`, `spatium::io::`, `spatium::spatial::`, `spatium::render::`, `spatium::physics::`, `spatium::viewer::`, `spatium::statistics::`; `rsc::` under `rsc/core`
- PascalCase classes, snake_case functions, trailing underscore for private members
- `Result<T> = std::expected<T, Error>` for fallible operations
- constexpr where possible, ADL-friendly math (using std::sqrt etc)
- Header-only, Boost optional (`SPATIUM_BOOST=ON`, multiprecision only), Eigen optional (`SPATIUM_EIGEN=ON` for heat method). Both follow the same pattern: a CMake option sets `SPATIUM_HAS_<DEP>=0/1`, and the headers that need the dependency compile to nothing when it is off — so a bare `-Iinclude` with no packages installed still builds `<spatium/core.hpp>` and `<spatium/spatium.hpp>`
- Catch2 v3 for tests
- Template params: `<std::size_t N, Scalar T = double>`
- Clean constructors: `Triangle3(a, b, c)` not `{{{a, b, c}}}`
- Operator convention: `|` = intersect/pipe, `&` = boolean intersect, `+` = union, `-` = difference
- Measure naming: the `Measurable` concept uses `measure()` as the dimension-generic name (length / area / volume / Hausdorff k-measure). Concrete 2D shapes (Triangle, Polygon, Disk, Circle, Box-2D) ship `area()` as a convenience alias; 1D shapes (Segment, Line) ship `length()`; 3D shapes will ship `volume()`. Aliases must always forward to `measure()`, never re-implement the formula.
