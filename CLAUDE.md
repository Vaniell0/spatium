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
`test_precision.cpp` and every `test_rsc_*` are dropped (RSC's registry
reaches `Real50` through `checkpoint.hpp`). Leave it off only to reproduce
what someone with no system packages gets — `cmake --preset noboost`, the
counterpart to `noeigen`, does exactly that and is the configuration a bare
`cmake -B build` with no flags now produces.

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
- `Real50 / Real100` — arbitrary precision (Boost.Multiprecision, requires `SPATIUM_BOOST=ON`)
- `Table / Svg` — structured output and 2D visualization
- `std::format` support for all types
- UDLs: `_deg`, `_pi`, `_x`, `_y`, `_z`

### CMake targets

- `Spatium::sdk` (INTERFACE) — header-only, link this
- `Spatium::core` (INTERFACE, will become SHARED with .cpp sources)
- `spatium` — backward-compat alias for sdk

## Directory layout

- `include/spatium/core/` — concepts, error, verify, precision, up_to.hpp (`UpTo<T, N>` — at most N values with a count, no allocation; the return type for an answer whose size is *bounded* rather than fixed, so a degenerate quadratic reports one root instead of dividing by zero. `operator[]` is unchecked under `NDEBUG`, `at()` returns `Result<T>` and never throws, Debug poisons the slots past the count through an ADL `debug_poison`)
- `include/spatium/algebra/` — Vec, Matrix, Quaternion, Complex, Dual (autodiff), calculus (gradient/integrate), general IVP solvers (ode.hpp), N×N linear solve (linear_solve.hpp: solve_direct/solve_jacobi), Eigen interop (eigen_interop.hpp, opt-in), polynomial solvers, generic functions (power/commutator/adjoint/poly_eval), noise.hpp (`PerlinNoise` — seeded gradient noise, deterministic per seed, no shared mutable state), literals, format. `inline namespace algebra` — see `docs/conventions.md`
- `include/spatium/algebra/groups/` — SO3 (rotation group, Rodrigues formula), SE3 (rigid-body motions), non-inline `spatium::algebra::` (see `docs/conventions.md` for why these stay qualified-only, unlike the rest of `algebra/`)
- `include/spatium/spaces/` — Euclidean, Sphere, Hyperbolic, ParametricSurface, ImplicitSurface, offset.hpp (`offset_surface()` — compose a new ParametricSurface pushed out along a base surface's own normal, analytic, no mesh), sample.hpp (`sample_surface_uniform()` — area-weighted rejection sampling via the first fundamental form's `area_element()`, analytic placement, no mesh/Voronoi graph), chart.hpp (`Chart` concept — `evaluate`/`normal_at`/`area_element`/`domain`/`periodic_u`/`periodic_v`, i.e. what the scene DSL's operations actually require, as opposed to `Surface`, which none of them calls; `chart_of(space)` is the ADL extension point that lets a space outside this library enter the DSL, and gives `Sphere<2,T>` its first parametrization)
- `include/spatium/geometry/` — primitives, intersection, distance, clip, boolean (polygon ops), make, transform (lazy chains), ray_surface (Quadric, ray_quadric), format
- `include/spatium/mesh/` — Mesh, subdivision, LOD, primitives, topology, geodesic, voronoi, transport, operations, conform.hpp (`conform_to_surface()` — drape a guide mesh onto a target Surface via project()+normal()), scatter.hpp (`scatter_on_surface()` — mesh-graph geodesic-Voronoi placement; spaces/sample.hpp's analytic version is preferred when the target is a ParametricSurface)
- `include/spatium/spatial/` — BVH (SAH build, ray_cast/ray_test with `t_max`, nearest, query_box; `BVH<Shape, Bound>` over any bound, a `Box` by default), bound.hpp (`Bound` concept — merge, lower_distance, sah_measure, ray_interval, center, overlaps, bound_of — with the box adapter and `Ball<N,T>`, the ambient ball a chart search or point cloud wants and the Euclidean case of a geodesic ball), geodesic_ball.hpp (`GeodesicBall<Space>` — a ball in a manifold's own metric; `ray_interval` for a unit-speed geodesic in closed form on Euclidean, Sphere (periodic, first entry) and Hyperbolic (the stable root in e^t), `lower_distance`, `merge` along the geodesic inside the injectivity radius, `sah_measure` by the boundary's measure, `geometry_bounds(space)` — injectivity radius and curvature range), ball_tree.hpp (`GeodesicBallTree<Space>` — geodesic balls split by two pivots with `distance` alone, nodes merged along geodesics; `ray_cast` of a unit-speed geodesic with `t_max` and `nearest`, exact on the sphere and in hyperbolic space)
- `include/spatium/render/` — CPU-raytracer engine shared across `examples/`: `Camera`/`make_camera_basis`/`camera_ray_dir` (pinhole camera), `parallel_for_rows()` (work-stealing row parallelism), `supersample_pixel()` (NxN jittered-grid antialiasing), `write_png_rgb`/`write_png_rgba` (image output), `blackbody_to_rgb255()` (physical temperature→RGB), `hsv_to_rgb255()` (non-physical color-picker model), texture.hpp (`Texture` — an RGB image loaded through the vendored stb_image or generated from a function, bilinear `sample(u, v)` with repeat or clamp: the sampling half of the ROADMAP's "a texture is a field plus a sampling call"), `Sky`/`make_starfield()`/`sample_sky_color()` (procedural starfield + nebulae; `wide_sky` toggles the whole-sky brightness gradient vs. a flat background for close single-object framing); cooked_scene.hpp (`lay_out()` — a `Cooked` scene as a world-space triangle tree plus a tree of instanced closed forms, the one layout both the CPU and GPU paths trace), gpu_types.hpp (the fp32 std430 layouts — `Node`, `Triangle`, `Quadric`, `Instance`, `LNode` — on their own so packing and tree building need not include each other), gpu_scene.hpp (`pack()` — the triangle tree as fp32 std430 arrays in leaf order and the instances under a linear BVH, plus `trace()`/`shade()`, the compute shader's specification written in plain fp32 on the host), gpu_trace_glsl.hpp (that shader, transcribed function for function), lbvh.hpp (`build_lbvh()` — a linear BVH (Karras 2012) over instances that move every frame: Morton codes, a parallel radix sort and a hierarchy decided per node, culling instances scaled to zero; host-parallel, the CPU path's tree and the spec a device build follows; `trace_lbvh()` walks it), blackhole_glsl.hpp (`blackhole_shader()` — the compute shader that renders a `SpacetimeScene`: rays traced backwards from an observer's tetrad through the metric's outgoing Kerr-Schild form, a step sized by the distance to the nearest horizon, capture by the Boyer-Lindquist radius, a disk whose emission falls as r^-3 and is shifted by g = (p.u_obs)/(p.u_emit), stars as pixel-wide Gaussians and a galactic band; `observer_tetrad()`, `blackbody_table()`), gpu_instances.hpp (`make_instance_kernel()` — a compute shader generated from a Scatter node's own fields that writes each instance's `gpu::Instance` slot on the device, placement and material, starting from `io::build::scatter_spots()`, the sites `cook()` uses)
- `include/spatium/discrete/` — FiniteSet, GeometricSet, combinatorics.hpp (`factorial`, `binomial_coefficient`, `permutations_count` — overflow-checked, `Result<T>`; `k_combinations` — `vector<FiniteSet<T>>`)
- `include/spatium/physics/` — Element (the one compiled, non-header-only translation unit)
- `include/spatium/physics/atomic/` — AtomModel, AtomPalette, BohrModel, orbital, atom_svg — visualization-support models, split out from physics/'s top level per `docs/conventions.md`'s subdivision rule
- `include/spatium/physics/mechanics/` — units (compile-time SI), PointMass/RigidBody, forces, integrators (Euler/RK4/Lie-group/LGVI), symplectic manifolds, variational integrators, geometric continuum mechanics, IPC contact barrier + XPBD, certified continuous collision (narrow_phase.hpp `distance_bound` — a lower bound on distance, exact for sphere/torus, branch and bound for a `LipschitzChart`, `|f|/L` for a `LipschitzImplicit`; rigid_contact.hpp `sweep_sphere_surface` advances by it, a miss only when proved), wave PDE time-stepping (wave_string.hpp — 1D clamped string, wave_membrane.hpp — 2D clamped rectangular membrane; both explicit central-difference/"leapfrog" steppers with closed-form fundamental-frequency helpers, purpose-built rather than routed through `algebra/ode.hpp` — see either header's own comment for why). Doesn't yet use `Result<T>` for fallible ops (see `docs/architecture.md`'s "Header-only spine, and three principled exceptions")
- `include/spatium/physics/relativity/` — metric-agnostic geodesic integration: `schwarzschild.hpp` and `kerr.hpp` (both metrics as templated callables, Dual<T>-substitutable; Kerr's non-diagonal g_tphi needed zero changes downstream since the metric inverse uses the general `solve_direct()`, not a diagonal-only shortcut), `geodesic.hpp` (exact Christoffel symbols via Dual<T> partials, `Vec<T,8>` state through `algebra/ode.hpp`'s `rk4_step`, Killing-vector conserved quantities), `metric_field.hpp` (`MetricField` — a metric as ten IR fields of `Field::coord(0..3)`, lowered and merged into one deduplicated pool, evaluated on any scalar so `Dual<T>` gives geodesic.hpp its derivatives unchanged; `kerr_boyer_lindquist()` matches `KerrMetric` entry for entry and `kerr_schild()` is the axis-regular Cartesian form, held to the vacuum equations by `vacuum_residual()` -- a finite-difference Ricci check -- against a non-vacuum control), `metric_glsl.hpp` (`emit_metric_glsl()` — a MetricField's pool as GLSL, each op a value and its four partials, forward-mode; `geodesic_glsl()` — Christoffel symbols through GLSL's inverse(mat4), the geodesic equation and an RK4 step, transcribed from geodesic.hpp, whose run on Dual<float> is the host specification), `spacetime_scene.hpp` (`SpacetimeScene` — a black-hole scene described by what is in it: holes with paths as motion fields, `binary()` circling or shrinking by the quadrupole formula to a stop separation, the metric as Minkowski plus one Kerr-Schild term per hole at its centre at the event's time, `residual_at()` reporting how far that superposition is from vacuum; disk, dust, sky and camera settings for the renderer), `accretion_disk.hpp` (Schwarzschild thin-disk redshift; Kerr's own BPT-1972 equatorial-orbit/ISCO/photon-orbit/redshift formulas live in `kerr.hpp` itself)
- `include/spatium/io/` — field.hpp (`Field`/`VecField` — a scene field is an expression over a flat pool of tagged ops in topological order, children before parents, with an opaque callable as a *leaf* rather than as an alternative to the expression; `is_structural()` is derived by looking for opaque leaves, never declared; a `Scale`'s or `Rotate`'s factor may be a callable or a `ScalarField<T>` expression, reading time as the field's first parameter (`ScalarField<T>::t()`) — the expression form is what lets a growing object stay structural, and the scalar vocabulary carries `Min`/`Max` for it, since without a clamp there is no `smoothstep` and therefore no growth curve that is not a lambda; a motion reads a named `MotionEnv` rather than a bare `t`, so an object's own clock or a solver-written pose later adds a struct field instead of changing every signature; `field_report()` counts recognized and unknown separately and prices captured state against L3; a field evaluated per instance reads `FieldInputs` -- the instance `id` (set per Scatter site in `MotionEnv`), its `origin` and the point -- through `Id`/`Origin`/`Point`, with `Hash` (an integer hash of the id, the same bits on host and device), `Sqrt` and `Gather` (a component of a point from a shared table), and `VecField::make()` builds a vector from three such fields, which is how the donut's dust became closure-free; `Coord` reads a spacetime point for a metric, and `in_spacetime()` puts a motion's time on the metric's clock -- coordinate 0 -- so a path written as a motion enters a metric unchanged), field_pod.hpp (`lower()` turns a `Field`/`VecField` into plain arrays of op codes, child indices and numbers plus its noise tables, failing with `Result` on an opaque leaf and naming it; `interpret()` evaluates that data alone, and `tests/test_field_pod.cpp` holds it to bit-exact agreement with `eval_into` over every field of the donut scene), field_glsl.hpp (`emit_scalar`/`emit_vector`/`emit_placement` — a lowered field written out as straight-line GLSL for a driver to compile, with Perlin noise, the instance hash, `SO3::exp` and the gather transcribed from their C++ definitions; checked on a device against cook() within a tolerance, since a GPU compiler may reorder). Table, SVG, OBJ, STL, JSON (json.hpp — hand-rolled parser/serializer, no external dependency), Scene (scene.hpp — JSON scene format: placed objects + optional camera, shape kinds resolved through an open registry rather than a closed enum, also home of `Material`), WAV (wav.hpp — hand-rolled 16-bit PCM writer, no external dependency, writer-only), build.hpp — `spatium::io::build`'s declarative scene DSL: a `Trace` is a flat, inspectable record of Space/Offset/Scatter/Compose operations over real spaces, not a tree of opaque closures. `materialize()` returns `Placed<T>` — a view onto a node, not a mesh: `surface()` gives the real `ParametricSurface` with the node's `.moving()` composed into the map, `mesh()` builds triangles only when something asks. Analytic end to end, so every manifold operation in the library applies to a scene object directly; see `docs/getting-started-dsl.md` and `examples/donut_demo.cpp`
- `include/spatium/viewer/` — Vulkan App (multi-mesh, point clouds, per-mesh color, ImGui), Camera, compute.hpp (headless Vulkan compute: `Context` picks a real GPU and refuses llvmpipe, host-visible storage `Buffer`s, a `Kernel` compiled from GLSL at run time, GPU time from timestamps; `Presenter` copies a kernel's pixel buffer into a window's swapchain and draws ImGui on top), gpu_lbvh.hpp (`DeviceLbvh` — the linear BVH built on the device over instances already there: boxes and Morton keys in one submission, the sort on the host in shared memory, then leaves, nodes and bounds; kernels in `render/gpu_lbvh_glsl.hpp`; instances below a projected size in pixels are left out of the tree and drawn by `render/gpu_splat_glsl.hpp` instead — projected, depth-tested against the trace, coverage summed as optical depth and composited with transmittance exp(-sum))
- `include/spatium/` — Point, Morphism, umbrella header
- `include/spatium/vendor/` — stb_image_write (screenshot export), stb_image (texture loading; both public domain, a caller defines the `*_IMPLEMENTATION` macro in one translation unit)
- `tests/` — Catch2; `ctest -N` reports the count, which depends on which optional dependencies are enabled
- `examples/` — geometry_demo, showcase, sets_demo, spd_relationship_demo, scene_demo (text/SVG); donut_demo (spatium::io::build's declarative DSL getting-started demo — console by default, `--photo` for PNG; see `docs/getting-started-dsl.md`; the scene itself is built in `donut_scene.hpp`, so a test can build the same one); atom_demo + primitives_demo (Vulkan; primitives_demo dispatches `--scene primitives|torus|klein`); blackhole_live (black holes described with `SpacetimeScene` and traced on the device; `--check` sends rays through the metric on the device, on the host in float and in double and reports the end points' spread and the device time; `--shadow-check` measures a Schwarzschild shadow's angular radius against 3 sqrt(3) M / D sqrt(1 - 2M/D); `--frame PATH` renders one frame), donut_live (the donut traced by a Vulkan compute shader: headless by default, checked against the host's fp32 trace and the fp64 trees; `--live` opens a window with a free-flying camera and ImGui; `--dust-check [--dust N]` moves the scene's Scatters on the device and compares every slot with `cook()`, then builds the instance tree on the host and on the device and checks the device tree against brute force; in `--live` every structural Scatter is moved by its kernel each frame and the tree rebuilt on the device, with play/rate in ImGui); blackhole_demo + blackhole_gr_demo, wormhole_demo, tumbling_body_demo, parametric_analytical_demo, geodesic_curvature_grid_demo, geodesic_procgen_demo, hyperbolic_tessellation_demo, wave_ca_demo, collatz_demo, burning_ship_demo, native_collision_demo, ball_pit_demo, sound_synthesis_demo, terminal_donut_demo (ASCII, `render::terminal_canvas.hpp`, not Vulkan), cloth_sphere_probe (offline/Vulkan raytracers and diagrams — see `spatium_add_example()` in `cmake/SpatiumTarget.cmake`)
- `benchmarks/` — Google Benchmark (vec, intersection, mesh, bvh, orbital, raycast)
- `cmake/` — SpatiumConfig.cmake.in, SpatiumModule.cmake (C++23 modules helper), SpatiumTarget.cmake (spatium_add_example() helper)
- `scripts/` — doc-freshness checks, all three run in CI so the docs can't silently drift from the real code again: `gen_dependency_graph.py --check` (real `#include` graph vs. the committed `.dot`), `check_claude_md_layout.py` (every `include/spatium/` subdirectory mentioned somewhere above), `check_doc_file_refs.py` (every source file the docs name actually exists; names with no counterpart — recorded deletions, proposals, generated files, other projects' — are listed with a reason in the script's `NOT_IN_TREE`)

## Conventions

- Namespace: `spatium::`, `spatium::algebra::` (inline — see `docs/conventions.md`), `spatium::geometry::`, `spatium::mesh::`, `spatium::io::`, `spatium::spatial::`, `spatium::render::`, `spatium::physics::`, `spatium::viewer::`
- PascalCase classes, snake_case functions, trailing underscore for private members
- `Result<T> = std::expected<T, Error>` for fallible operations
- constexpr where possible, ADL-friendly math (using std::sqrt etc)
- Header-only, Boost optional (`SPATIUM_BOOST=ON`, multiprecision only), Eigen optional (`SPATIUM_EIGEN=ON` for heat method). Both follow the same pattern: a CMake option sets `SPATIUM_HAS_<DEP>=0/1`, and the headers that need the dependency compile to nothing when it is off — so a bare `-Iinclude` with no packages installed still builds `<spatium/core.hpp>` and `<spatium/spatium.hpp>`
- Catch2 v3 for tests
- Template params: `<std::size_t N, Scalar T = double>`
- Clean constructors: `Triangle3(a, b, c)` not `{{{a, b, c}}}`
- Operator convention: `|` = intersect/pipe, `&` = boolean intersect, `+` = union, `-` = difference
- Measure naming: the `Measurable` concept uses `measure()` as the dimension-generic name (length / area / volume / Hausdorff k-measure). Concrete 2D shapes (Triangle, Polygon, Disk, Circle, Box-2D) ship `area()` as a convenience alias; 1D shapes (Segment, Line) ship `length()`; 3D shapes will ship `volume()`. Aliases must always forward to `measure()`, never re-implement the formula.
