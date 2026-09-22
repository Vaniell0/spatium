<p align="center"><img src="gallery/hero.png" alt="Spatium" width="800"></p>

# Spatium

**Any mathematical space, one library.**

Geometry libraries hardcode their space: CGAL's kernels, Eigen's linear algebra, GLM's vectors all assume flat Euclidean R^N wired into every type. Need geodesics on a sphere, distances in hyperbolic space, mesh operations on some other manifold? That's a different, specialized library each time — or a parallel hand-written stack duplicating the one you already have.

Spatium doesn't hardcode a space. It has a concept hierarchy — Set → TopologicalSpace → MetricSpace → NormedSpace → InnerProductSpace → Manifold → RiemannianManifold → Surface — and any type satisfying a concept's requirements gets the whole library for free:

```cpp
struct FlatTorus { /* distance(), exp_map(), log_map(), project()... */ };
static_assert(spatium::RiemannianManifold<FlatTorus>);
// Mesh<FlatTorus>, subdivision, geodesics, morphisms — all work automatically.
```

One exception to "all work automatically", stated here rather than left to
be discovered: **spatial acceleration is still flat.** `spatial/`'s BVH
bounds with axis-aligned boxes, so ray casting and nearest-neighbour
queries are accelerated in Euclidean space and unaccelerated off it. The
operations still give correct answers on any space; they just walk
everything. A ball tree over geodesic balls is the fix and is an open item
in the [roadmap](docs/ROADMAP.md), not an oversight.

C++23, in large part header-only — three deliberate exceptions exist where real complexity made that the wrong tradeoff, not an oversight: the Vulkan viewer needs genuine C linkage, the periodic-table data backs a single compiled translation unit, and the `physics/mechanics` research track plus optional CUDA/ipc-toolkit integrations sit outside the header-only spine on purpose. See [Architecture](docs/architecture.md#header-only-spine-and-three-principled-exceptions) for the honest breakdown, not a marketing gloss.

## Gallery

- [Kerr black hole](gallery/blackhole_gr.png) — full 4-coordinate geodesic integration, GPU-rendered (CUDA) at 1920x1080 ([video](gallery/blackhole_gr.mp4))
- [A donut, declaratively](gallery/donut_dsl.png) — built entirely from `torus()`/`offset()`/`scatter()`, see the [getting-started guide](docs/getting-started-dsl.md) ([build-up video](gallery/donut_dsl.mp4))

The donut is also where the scene DSL's argument is easiest to check. A
scene is described as *spaces* rather than as meshes, and the description
stays a small inspectable graph: **34 nodes describing 2,021,984 objects**.
Geometry that would be 64,654,768 vertices if every object carried its own
copy is stored as 39,272 — about **1646x** — and the frame renders in under
a gigabyte. That is not a trick in the renderer; it is what having
described the scene as spaces buys.

More in [`gallery/`](gallery/).

## Features

**Spaces.** The concept hierarchy — Set, TopologicalSpace, MetricSpace,
NormedSpace, InnerProductSpace, Manifold, RiemannianManifold, Surface — and
the spaces that satisfy it: Euclidean\<N\>, Sphere\<N\>, Hyperbolic\<N\>,
ParametricSurface, ImplicitSurface, product spaces. `chart_of()` is the ADL
extension point that lets a space declared outside this library into the
scene DSL.

**Lie groups and matrix manifolds.** SO(3) via Rodrigues, SE(3) for
rigid-body motions, and SPD(n) — the manifold of symmetric positive-definite
matrices — under two metrics, log-Euclidean and affine-invariant. All
templated on the scalar type.

**Geometry.** Line, Ray, Segment, Hyperplane, Triangle, Polygon, Circle,
Disk, Box, Simplex, Quadric. Intersection (Möller–Trumbore, slab method,
analytical ray–quadric and ray–torus), distance between every shape pair,
polygon booleans, clipping, and a `|` pipe syntax for composing them.

**Mesh and geodesics.** `Mesh<Surface>` for any surface, subdivision with
surface projection, LOD chains, geodesic distance by Dijkstra or the heat
method (Crane 2013), geodesic Voronoi, parallel transport, and discrete
exterior calculus — cotangent Laplacian, mass matrix, face gradients,
divergence.

**Spatial acceleration.** BVH with SAH construction, ray casting, nearest
queries, box queries, over triangles or over analytic primitives directly.

**Numerics.** `Dual<T>` forward-mode autodiff that satisfies `Scalar` and so
substitutes into ordinary code; polynomial solvers through the quartic;
native SVD and symmetric eigendecomposition; N×N linear solve; ODE
integrators; `Complex<T>`; calculus over plain callables — gradient,
integrate, minimize. Arbitrary precision through Boost.Multiprecision
(`Real50`, `Real100`, any digit count), optional.

**Riemannian optimization.** Index-raise an ambient covector through the
space's own metric, project to the tangent space, retract by `exp_map` — so
gradient descent works on any RiemannianManifold that is also a Surface.

**Physics.** Compile-time SI units, point masses and rigid bodies,
symplectic and Lie-group and variational integrators, geometric continuum
mechanics, the IPC contact barrier with continuous collision detection
against any Surface, XPBD, wave PDEs for strings and membranes, and
metric-agnostic geodesic integration with exact Christoffel symbols —
Schwarzschild and Kerr, with the accretion-disk redshift for both.

**Declarative scene DSL** (`io::build`). `torus()`/`offset()`/`scatter()`/
`compose()` build a flat, inspectable `Trace` rather than a tree of opaque
closures, analytic until the last mile: offset surfaces and area-weighted
placement are real function composition, and no mesh exists until something
asks for triangles. Tutorial in [`docs/getting-started-dsl.md`](docs/getting-started-dsl.md),
runnable in [`examples/donut_demo.cpp`](examples/donut_demo.cpp).

**Rendering and I/O.** A CPU ray tracer — pinhole camera, row-parallel
work stealing, jittered supersampling, PNG output, physical blackbody
colour — plus a Vulkan viewer with ImGui. Hand-rolled readers and writers
with no external dependency: JSON, WAV, OBJ, STL, SVG.

**RSC.** A trained dispatcher that picks which implementation to use for a
problem — Newton or bisection, double or fifty digits, analytic or
tessellated — with ground truth generated by running the candidates rather
than hand-labelled. Seven domains; see the [roadmap](docs/ROADMAP.md) for
what is and is not connected.

**Verification.** `verify_metric`, `verify_inner_product`, `verify_exp_log`
and a symplecticity check, so a space you define yourself can be tested
against the axioms it claims.

**N-dimensional and zero-cost.** Templated on dimension and scalar type,
concepts checked at compile time, no virtual dispatch anywhere.

## Quick Start

```cpp
#include <spatium/spatium.hpp>
#include <print>

using namespace spatium;
using namespace spatium::geometry;

int main() {
    // Geometry — clean factory syntax
    auto t = tri(Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0});
    std::println("area = {:.4f}, normal = {}", t.area(), t.normal());

    // Intersection via pipe
    auto r = *ray(Vec3{0.25, 0.25, 5}, Vec3{0, 0, -1});
    if (auto hit = r | t)
        std::println("hit at {}", *hit);

    // Morphism pipeline
    auto scale = morph<E3, E3>([](const Vec3& p) { return p * 2.0; });
    auto proj  = morph<E3, E2>([](const Vec3& p) -> Vec2 { return {p[0], p[1]}; });
    auto result = pt<E3>(Vec3{1, 2, 3}) | scale | proj;
    std::println("{}", result);  // P(2, 4)

    // Sphere geodesics
    S2 sphere;
    auto north = pt<S2>(Vec3{0, 0, 1});
    auto east  = pt<S2>(Vec3{1, 0, 0});
    auto tangent = north.log(east, sphere);
    auto midpoint = north.exp(tangent, 0.5, sphere);
    std::println("geodesic midpoint: {}", midpoint);

    // Mesh subdivision
    auto mesh = mesh::icosahedron(sphere);
    auto refined = mesh::subdivide(mesh, sphere, 3);
    std::println("{}", refined);  // Mesh{V=642 F=1280 E≈1920}
}
```

## Build

Requires C++23 (GCC 15+ or Clang 19+), CMake 3.28+, Catch2 v3 for tests.

```bash
# With Nix (recommended)
nix develop
cmake --preset default
cmake --build --preset default
ctest --preset default

# Without Nix
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
ninja -C build
```

`CMakePresets.json` has presets beyond `default` for common configurations — `release` (Eigen, for RSC training-heavy work), `modules` (the C++23 modules build path), `noeigen`, `vulkan-dev`, `cuda`, and two benchmark-harness presets. `cmake --list-presets` shows all of them.

### CMake options

| Option                       | Default | Notes |
|------------------------------|:-------:|-------|
| `SPATIUM_BUILD_TESTS`        | `ON`*   | Catch2 v3 unit tests (`ctest --preset default`) |
| `SPATIUM_BUILD_EXAMPLES`     | `ON`*   | All `examples/*` binaries |
| `SPATIUM_BUILD_BENCHMARKS`   | `OFF`   | Google Benchmark suite (`benchmarks/`) |
| `SPATIUM_BUILD_VIEWER`       | `ON`*   | Vulkan + GLFW + shaderc viewer. Emits a `WARNING` and skips the target if any of those packages are missing. |
| `SPATIUM_EIGEN`              | `OFF`   | Required by the heat-method geodesic solver and the cotangent-Laplacian DEC operators. |
| `SPATIUM_NATIVE_ARCH`        | `OFF`   | Adds `-march=native`. Resulting binaries are not portable across CPUs — use only for local performance work. |
| `SPATIUM_USE_MODULES`        | `OFF`   | C++23 named-modules build. Currently behind the header tree (see the option's own comment in `CMakeLists.txt`) — not a compiler-bug wait, real catch-up work. |
| `SPATIUM_IPC_TOOLKIT`        | `OFF`   | Implicit contact physics via [ipc-toolkit](https://github.com/ipc-sim/ipc-toolkit) (Newton + log-barrier + CCD). `FetchContent`-based, pulls its own dependency tree. |
| `SPATIUM_CUDA`                | `OFF`   | CUDA GPU kernels (`gpu/`) for GR ray tracing. Requires nvcc; not part of a default build. |
| `SPATIUM_BUILD_RSC_TOOLS`     | `ON`*   | RSC training tools (`rsc/tools/train_base`, ...). |
| `IMGUI_DIR` (env or `-D`)    | unset   | Source path of Dear ImGui; enables the in-viewer panel when set. |

\* Defaults to `ON` only when Spatium is the top-level CMake project (built
standalone, as above). Pulled in via `add_subdirectory()` or `FetchContent`
from another project, these four default to `OFF` instead, so a downstream
consumer gets just `Spatium::sdk` without forcing a Vulkan/Catch2/example
build it never asked for -- see "Using in Your Project" below.

## Using in Your Project

### CMake FetchContent

```cmake
include(FetchContent)
FetchContent_Declare(spatium
    GIT_REPOSITORY https://github.com/Vaniell0/spatium.git
    GIT_TAG v1.0.0
)
FetchContent_MakeAvailable(spatium)

target_link_libraries(your_target PRIVATE Spatium::sdk)
```

This pulls in only the header-only `Spatium::sdk` interface target -- the
Vulkan viewer, examples, tests, and RSC tools all default OFF when Spatium
isn't the top-level CMake project, so nothing beyond `Spatium::sdk` and its
one required dependency (Boost headers, for `Real50`/`Real100`) gets built.
See [`examples/external-consumer/`](examples/external-consumer/) for a
complete, independently-buildable project using exactly this snippet.

### After Install

```cmake
find_package(Spatium REQUIRED)
target_link_libraries(your_target PRIVATE Spatium::sdk)
```

## Defining Your Own Space

Any struct with the right methods satisfies the concepts automatically:

```cpp
struct FlatTorus {
    using ScalarType = double;
    using PointType = Vec<double, 2>;
    using TangentVector = Vec<double, 2>;
    static constexpr std::size_t dimension = 2;
    static constexpr bool is_complete = true;

    bool contains(const PointType& p) const { /* ... */ }
    ScalarType distance(const PointType& a, const PointType& b) const { /* ... */ }
    PointType exp_map(const PointType& p, const TangentVector& v, ScalarType t) const { /* ... */ }
    TangentVector log_map(const PointType& p, const PointType& q) const { /* ... */ }
    ScalarType metric_at(const PointType& p, const TangentVector& u, const TangentVector& v) const { /* ... */ }
    PointType project(const PointType& p) const { /* ... */ }
    TangentVector normal(const PointType& p) const { /* ... */ }
};

static_assert(spatium::RiemannianManifold<FlatTorus>);
static_assert(spatium::Surface<FlatTorus>);
// Mesh<FlatTorus>, subdivision, morphisms — all work automatically.
```

## Arbitrary Precision

```cpp
#include <spatium/core/precision.hpp>

using namespace spatium;

// 50-digit precision
Euclidean<3, Real50> space;
Vec<Real50, 3> a{Real50{0}, Real50{0}, Real50{0}};
Vec<Real50, 3> b{Real50{3}, Real50{4}, Real50{0}};
auto d = space.distance(a, b);  // 5.000...000 (50 digits)
```

## Documentation

- [Architecture](docs/architecture.md) — concept hierarchy, design decisions, the real dependency graph
- [Conventions](docs/conventions.md) — namespace/subdivision/error-handling rules, and the known violations being fixed
- [API Reference](docs/api-reference.md) — all types, methods, concepts
- [Quick Start Guide](docs/quickstart.md) — getting started
- [Getting Started: The Declarative Scene DSL](docs/getting-started-dsl.md) — zero-barrier-to-entry, build a donut in three declarative steps
- [Extending Spatium](docs/extending.md) — defining custom spaces and primitives
- [Roadmap](docs/ROADMAP.md) — what's done, what's planned, project history
- [Concept-Driven Physics](docs/concept-driven-physics.md) — how `physics/mechanics/` fits the concept hierarchy

## Project Structure

```
include/spatium/
    core/                concepts, error, verify, precision
    algebra/             Vec, Matrix, Quaternion, Complex, Dual (autodiff), calculus,
                         ODE solvers, linear solve, polynomial solvers, Eigen interop
    algebra/groups/      SO3, SE3
    spaces/              Euclidean, Sphere, Hyperbolic, ParametricSurface, ImplicitSurface
    geometry/            primitives, intersection, distance, boolean ops, ray_surface (Quadric)
    mesh/                Mesh, subdivision, LOD, topology, geodesic, voronoi, DEC
    spatial/             BVH (SAH build, ray_cast, nearest, query_box)
    discrete/            FiniteSet, GeometricSet
    render/              supersample_pixel(), camera, write_image, parallel_for_rows
    io/                  Table, SVG, OBJ, STL
    physics/             periodic-table element data (the one compiled TU)
    physics/atomic/      atom/orbital models, Bohr model, SVG rendering
    physics/mechanics/   integrators, symplectic/Lie-group/variational structure, contact
    physics/relativity/  Schwarzschild/Kerr geodesic integration, accretion disks
    viewer/              Vulkan app (multi-mesh, point clouds, ImGui)
    point.hpp, morphism.hpp, spatium.hpp
rsc/                     RSC — trained dispatcher on top of Spatium (7 domains, see docs/ROADMAP.md)
tests/, examples/, benchmarks/
```

### Demos

```bash
nix run .#primitives                    # unified primitives + BVH raycast, interactive Vulkan
nix run .#tumbling                      # Dzhanibekov-effect rigid-body tumble (LGVI), frame sequence
```

More demos exist in `examples/` — analytical ray tracing, a Schwarzschild/Kerr GR raytracer, an Ellis wormhole flythrough, and others; `cmake --list-presets` and `nix flake show` list every buildable target.

## License

Apache License 2.0 — see [`LICENSE`](LICENSE). Contributing guide: [`CONTRIBUTING.md`](CONTRIBUTING.md).
Project history: [`docs/ROADMAP.md`](docs/ROADMAP.md).
