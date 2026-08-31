# Spatium

C++23 header-only math library for arbitrary mathematical spaces, geometric primitives, and mesh operations.

## Build

```bash
nix develop
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
ninja -C build
ctest --test-dir build
./build/examples/showcase
```

A second config, `build-release/` (`-DCMAKE_BUILD_TYPE=Release -DSPATIUM_EIGEN=ON`,
modules off), exists alongside `build/` for RSC training-heavy work.
Configure/build it the same way, pointed at `build-release` instead of `build`.
Note: measured no real wall-clock difference vs. Debug on the actual project
build for RSC's REINFORCE training loop (~22s either way) — a standalone
`-O0` vs `-O2` comparison suggested a ~5.7x win that didn't reproduce here
(cause not tracked down, possibly modules-related). Keep the config anyway;
don't assume it's a speed fix without re-measuring.

## Architecture

### Spaces (concept hierarchy)

```
Set → TopologicalSpace → MetricSpace → NormedSpace → InnerProductSpace
                    ↓                                        ↓
              Manifold → RiemannianManifold          Complete + Inner = Hilbert
                    ↓                                Complete + Normed = Banach
              Surface (+ project/normal)
```

Concrete: Euclidean<N>, Sphere<N>, Hyperbolic<N>

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
- `Real50 / Real100` — arbitrary precision (Boost.Multiprecision)
- `Table / Svg` — structured output and 2D visualization
- `std::format` support for all types
- UDLs: `_deg`, `_pi`, `_x`, `_y`, `_z`

### CMake targets

- `Spatium::sdk` (INTERFACE) — header-only, link this
- `Spatium::core` (INTERFACE, will become SHARED with .cpp sources)
- `spatium` — backward-compat alias for sdk

## Directory layout

- `include/spatium/core/` — concepts, error, verify, precision
- `include/spatium/algebra/` — Vec, Matrix, Quaternion, Complex, Dual (autodiff), calculus (gradient/integrate), general IVP solvers (ode.hpp), N×N linear solve (linear_solve.hpp: solve_direct/solve_jacobi), Eigen interop (eigen_interop.hpp, opt-in), polynomial solvers, generic functions (power/commutator/adjoint/poly_eval), literals, format. `inline namespace algebra` — see `docs/conventions.md`
- `include/spatium/algebra/groups/` — SO3 (rotation group, Rodrigues formula), SE3 (rigid-body motions), non-inline `spatium::algebra::` (see `docs/conventions.md` for why these stay qualified-only, unlike the rest of `algebra/`)
- `include/spatium/spaces/` — Euclidean, Sphere, Hyperbolic, ParametricSurface, ImplicitSurface
- `include/spatium/geometry/` — primitives, intersection, distance, clip, boolean (polygon ops), make, transform (lazy chains), ray_surface (Quadric, ray_quadric), format
- `include/spatium/mesh/` — Mesh, subdivision, LOD, primitives, topology, geodesic, voronoi, transport, operations
- `include/spatium/spatial/` — BVH (SAH build, ray_cast, nearest, query_box)
- `include/spatium/render/` — CPU-raytracer engine shared across `examples/`: `Camera`/`make_camera_basis`/`camera_ray_dir` (pinhole camera), `parallel_for_rows()` (work-stealing row parallelism), `supersample_pixel()` (NxN jittered-grid antialiasing), `write_png_rgb`/`write_png_rgba` (image output), `blackbody_to_rgb255()` (physical temperature→RGB), `hsv_to_rgb255()` (non-physical color-picker model), `Sky`/`make_starfield()`/`sample_sky_color()` (procedural starfield + nebulae; `wide_sky` toggles the whole-sky brightness gradient vs. a flat background for close single-object framing)
- `include/spatium/discrete/` — FiniteSet, GeometricSet
- `include/spatium/physics/` — Element (the one compiled, non-header-only translation unit)
- `include/spatium/physics/atomic/` — AtomModel, AtomPalette, BohrModel, orbital, atom_svg — visualization-support models, split out from physics/'s top level per `docs/conventions.md`'s subdivision rule
- `include/spatium/physics/mechanics/` — units (compile-time SI), PointMass/RigidBody, forces, integrators (Euler/RK4/Lie-group/LGVI), symplectic manifolds, variational integrators, geometric continuum mechanics, IPC contact barrier + XPBD. Doesn't yet use `Result<T>` for fallible ops (see `docs/architecture.md`'s "Header-only spine, and three principled exceptions")
- `include/spatium/physics/relativity/` — metric-agnostic geodesic integration: `schwarzschild.hpp` and `kerr.hpp` (both metrics as templated callables, Dual<T>-substitutable; Kerr's non-diagonal g_tphi needed zero changes downstream since the metric inverse uses the general `solve_direct()`, not a diagonal-only shortcut), `geodesic.hpp` (exact Christoffel symbols via Dual<T> partials, `Vec<T,8>` state through `algebra/ode.hpp`'s `rk4_step`, Killing-vector conserved quantities), `accretion_disk.hpp` (Schwarzschild thin-disk redshift; Kerr's own BPT-1972 equatorial-orbit/ISCO/photon-orbit/redshift formulas live in `kerr.hpp` itself)
- `include/spatium/io/` — Table, SVG, OBJ, STL
- `include/spatium/viewer/` — Vulkan App (multi-mesh, point clouds, per-mesh color, ImGui), Camera
- `include/spatium/` — Point, Morphism, umbrella header
- `include/spatium/vendor/` — stb_image_write (screenshot export)
- `tests/` — Catch2 (790+ tests)
- `examples/` — geometry_demo, showcase, sets_demo (text); atom_demo + primitives_demo (Vulkan; primitives_demo dispatches `--scene primitives|torus|klein`); blackhole_demo + blackhole_gr_demo, wormhole_demo, tumbling_body_demo, parametric_analytical_demo, geodesic_curvature_grid_demo, geodesic_procgen_demo, hyperbolic_tessellation_demo, wave_ca_demo, collatz_demo, burning_ship_demo (offline/Vulkan raytracers and diagrams — see `spatium_add_example()` in `cmake/SpatiumTarget.cmake`)
- `benchmarks/` — Google Benchmark (vec, intersection, mesh, bvh, orbital, raycast)
- `cmake/` — SpatiumConfig.cmake.in, SpatiumModule.cmake (C++23 modules helper), SpatiumTarget.cmake (spatium_add_example() helper)
- `scripts/` — doc-freshness checks, both run in CI so this file and `docs/dependency-graph.dot` can't silently drift from the real code again: `gen_dependency_graph.py --check` (real `#include` graph vs. the committed `.dot`), `check_claude_md_layout.py` (every `include/spatium/` subdirectory mentioned somewhere above)

## Conventions

- Namespace: `spatium::`, `spatium::algebra::` (inline — see `docs/conventions.md`), `spatium::geometry::`, `spatium::mesh::`, `spatium::io::`, `spatium::spatial::`, `spatium::render::`, `spatium::physics::`, `spatium::viewer::`
- PascalCase classes, snake_case functions, trailing underscore for private members
- `Result<T> = std::expected<T, Error>` for fallible operations
- constexpr where possible, ADL-friendly math (using std::sqrt etc)
- Header-only, Boost optional (multiprecision only), Eigen optional (SPATIUM_EIGEN=ON for heat method)
- Catch2 v3 for tests
- Template params: `<std::size_t N, Scalar T = double>`
- Clean constructors: `Triangle3(a, b, c)` not `{{{a, b, c}}}`
- Operator convention: `|` = intersect/pipe, `&` = boolean intersect, `+` = union, `-` = difference
- Measure naming: the `Measurable` concept uses `measure()` as the dimension-generic name (length / area / volume / Hausdorff k-measure). Concrete 2D shapes (Triangle, Polygon, Disk, Circle, Box-2D) ship `area()` as a convenience alias; 1D shapes (Segment, Line) ship `length()`; 3D shapes will ship `volume()`. Aliases must always forward to `measure()`, never re-implement the formula.
