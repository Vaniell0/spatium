# Architecture

## Overview

Spatium models mathematics through two orthogonal axes:

1. **Algebraic structure** — what operations exist (addition, inner product, group composition)
2. **Topological/geometric structure** — what "shape" the space has (distance, curvature, geodesics)

These are expressed as C++23 concepts that compose freely. A space becomes a Hilbert space not by inheriting from a `HilbertSpace` class, but by satisfying both `InnerProductSpace` and `Complete`.

## Header-only spine, and three principled exceptions

Spatium is not one uniform thing — it's a header-only core (the concept
hierarchy above, plus `algebra/`, `spaces/`, `geometry/`, `mesh/`,
`spatial/`, `discrete/`, `render/`) unified by the concepts on this page,
plus three parts that sit outside that spine for concrete, documented
reasons rather than by accident:

- **`physics/mechanics/`, `spaces/`, and `mesh/`'s internal operations**
  — all three use `assert()` on internal invariants rather than
  `Result<T>`, unlike `geometry/`, `io/`, and `discrete/`, which use
  `Result<T>` at their boundary functions (construction/parsing from
  arbitrary input). This is a real split, not an oversight: see
  `docs/conventions.md`'s error-handling section for the boundary-vs-
  internal-invariant rule that explains it. `physics/mechanics/` was
  the first domain where this was named (see `docs/
  concept-driven-physics.md`) — `spaces/` and `mesh/` follow the same
  logic without previously being listed here.
- **`viewer/`** — Vulkan/GLFW/ImGui require genuine C linkage, which is
  incompatible with the `SPATIUM_EXPORT`/C++ modules pattern every other
  header uses. It's deliberately excluded from the modules build (see
  the top-level `CMakeLists.txt`), not an oversight.
- **`physics/elements.hpp`** — backed by one compiled translation unit,
  `src/physics/elements.cpp` (linked via the `spatium_core` shared
  library), the one non-header-only file in the entire tree. A large
  static periodic-table dataset is exactly the kind of thing header-only
  won't always fit cleanly, and this is a deliberate, narrow exception
  rather than a reason to force everything into one mold.

Naming the real shape here is the point: an agent or contributor running
into one of these three seams should find it documented, not surprising.

## Concept Hierarchy

### Space Concepts (`core/concepts.hpp`)

```
Set                         — PointType, ScalarType, dimension, equality
 └── TopologicalSpace       — contains(point)
      ├── MetricSpace        — distance(p, q)
      │    └── (+ VectorSpace + norm) → NormedSpace
      │         └── (+ inner) → InnerProductSpace
      │              └── (+ Complete) → HilbertSpace
      │                   └── (+ finite dim) → EuclideanSpace
      └── Manifold           — TangentVector, exp_map, log_map
           ├── RiemannianManifold — metric_at(point, u, v)
           └── Surface            — project(point), normal(point)
```

Orthogonal composition (all in `core/concepts.hpp`):
- `VectorSpace` — `PointType` supports linear ops via `VectorType`
- `Complete` — opt-in tag: `S::is_complete == true`
- `BanachSpace = NormedSpace + Complete`
- `HilbertSpace = InnerProductSpace + Complete`
- `EuclideanSpace = HilbertSpace + finite dimension`

### Geometry Concepts (`geometry/concepts.hpp`)

```
Shape                — ScalarType, PointType, ambient_dimension, centroid()
 ├── ClosedShape     — + contains(point)
 ├── Measurable      — + measure()
 ├── Bounded         — + bounding_box()
 ├── DistanceQueryable — + distance(point), project(point)
 └── BoundedRegion   — + subspace() → Result<SubspaceType>
```

Used as constraints:
- `Shape<A> && Shape<B>` in `intersection_region()` (`boolean.hpp`)
- `Measurable<A> && Measurable<B>` in `difference_area()` (`boolean.hpp`)
- `BoundedRegion<A> && BoundedRegion<B>` in `intersect_via_subspace()` (`intersection.hpp`)

### Algebra Concepts (`algebra/concepts.hpp`)

```
Magma        — compose(a, b)
 └── Semigroup (associative)
      └── Monoid     — + identity()
           └── Group — + inverse(a)
                └── LieGroup — + AlgebraType, exp(v), log(a)

Ring         — add, negate, zero, multiply, one
 └── Field  — + reciprocal(a)

LieAlgebra   — bracket(x, y), add(x, y), scale(s, x)
```

Used as constraints in generic functions (`algebra/functions.hpp`):
- `power(g, a, n)` requires `Group<G>`
- `commutator(g, a, b)` requires `Group<G>`
- `adjoint(g, elem, v)` requires `LieGroup<G>`
- `poly_eval(ring, coeffs, x)` requires `Ring<R>`

Concrete implementations: `SO3` (`algebra/groups/so3.hpp`), `SE3` (`algebra/groups/se3.hpp`).

### Physics Concepts (`physics/mechanics/*.hpp`)

```
SymplecticManifold    — Configuration, Momentum, State, dimension
                        (physics/mechanics/symplectic.hpp)
DiscreteLagrangian    — callable (q_k, q_{k+1}, h) → scalar + derivatives
                        (physics/mechanics/variational.hpp)
StrainEnergy          — W(F) energy + PK1 stress
                        (physics/mechanics/continuum.hpp)
ContactSurface<S, T>  — free point_to(p, s) → ContactQuery<T>
                        (physics/mechanics/narrow_phase.hpp)
```

The `SymplecticManifold<S>` concept is satisfied by `CotangentBundle<M>` for any `Manifold M`; the canonical 2-form ω = dq ∧ dp is synthesised there. `verify_symplecticity_drift` is the runtime axiom check (symplectic forms cannot be expressed in `requires`).

`DiscreteLagrangian` captures the Marsden-West discrete variational principle. `SeparableMidpointLagrangian` is the canonical implementation; `variational_step_separable` is Störmer-Verlet derived in closed form from the discrete Euler-Lagrange equation of the trapezoidal discrete action.

`ContactSurface` is the dispatch vehicle: `point_to(p, s)` has overloads for `Sphere<2>`, `Torus`, `ParametricSurface`, plus a generic `Surface`-fallback using `project` + `normal`. The compiler selects statically by argument type — no virtuals, no runtime tag. User-defined shapes join the concept by supplying their own `point_to` overload (ADL extension). See `docs/concept-driven-physics.md` for the full specification and rationale.

Dispatch shape:

```cpp
template<Scalar T>
ContactQuery<T> point_to(const Vec<T,3>& p, const Sphere<2,T>& s);      // ~3.7 ns/op
template<Scalar T>
ContactQuery<T> point_to(const Vec<T,3>& p, const Torus<T>& t);         // ~17.7 ns/op
template<Scalar T>
ContactQuery<T> point_to(const Vec<T,3>& p, const ParametricSurface<T>& s); // ~5 µs/op
template<Surface S>
ContactQuery<typename S::ScalarType> point_to(
    const Vec<typename S::ScalarType,3>& p, const S& surf);             // generic Newton
```

## Design Decisions

### Why Concepts (not Inheritance)

C++ concepts use **structural typing**: if a type has the required methods, it satisfies the concept. No `class MySpace : public MetricSpace<...>` boilerplate.

- Mathematical spaces combine properties **orthogonally** — inheritance forces a single hierarchy.
- Zero runtime cost — concepts are compile-time only.
- Error messages say "concept X not satisfied" instead of template error walls.
- Users define a struct, implement methods, everything works.

### Space = Type + Value

```cpp
Sphere<2> unit_sphere;                    // radius = 1.0 (default)
Sphere<2> big_sphere{.radius = 5.0};      // radius = 5.0
```

Both are the same type (`Sphere<2, double>`) but different instances. The type satisfies concepts at compile time; the value carries runtime parameters. `Mesh<Sphere<2>>` works for any radius.

### Points in Ambient Coordinates

`Sphere<N>` stores points as `Vec<T, N+1>` — coordinates in the embedding Euclidean space. Intrinsic coordinates (spherical coords) introduce chart singularities (poles) that complicate everything. Ambient coordinates keep `project()`, `distance()`, and vector operations uniform.

### Geometry Primitives are Standalone

`Triangle<3, double>` is not tied to a space. It works with `Vec<T, N>` directly. Shapes parameterized by dimension, not by space. A triangle "on" a sphere is a mesh face, not a Triangle type. The mesh system handles curved surfaces; primitives stay simple.

### Intersection as Free Functions

```cpp
auto hit = intersect(ray, triangle);  // Result<Vec3>      (intersection.hpp)
auto line = intersect(plane1, plane2); // Result<Line3>     (intersection.hpp)
```

Not methods — because the return type depends on *both* operands. Free functions with overloads handle this naturally. `operator|` in `geometry/make.hpp` dispatches to `intersect()`.

### Result<T> = std::expected<T, Error>

All fallible operations return `Result<T>` (`core/error.hpp`). Error codes in `ErrorCode` enum: `NoIntersection`, `DegenerateInput`, `SingularMatrix`, etc.

Pipe-unwrap: `Result<A> | B` auto-unwraps A before applying B:
- `Result<A> | B` for intersection (`geometry/make.hpp`)
- `Result<Vec> | AffineTransform` (`geometry/transform.hpp`)
- `Result<Point> | Morphism` (`morphism.hpp`)

### ADL-Friendly Math

```cpp
using std::sqrt;
auto n = sqrt(x);  // finds std::sqrt for double, boost::multiprecision::sqrt for Real50
```

All types work with arbitrary precision without changing a single line.

## Module Map

```
core/                                   — concepts, error, verify, precision, epsilon
  concepts.hpp          161 lines        14 space concepts (Set through Surface)
  error.hpp              56 lines        ErrorCode, Error, Result<T>
  verify.hpp            189 lines        verify_metric(), verify_inner_product(), verify_exp_log()
  precision.hpp          19 lines        Real50, Real100
  epsilon.hpp            50 lines        epsilon<T>(), approx_equal()

algebra/                                — linear algebra, complex, polynomials, groups
                                          (inline namespace algebra — see docs/conventions.md)
  vector.hpp            282 lines        Vec<T,N> — constexpr, SIMD, expression templates
  matrix.hpp            248 lines        Matrix<T,R,C> — column-major, inverse()→Result
  quaternion.hpp        168 lines        Quaternion<T> — axis-angle, slerp, from_matrix
  complex.hpp            96 lines        Complex<T> — arithmetic, sqrt, cbrt, from_polar
  polynomial.hpp        202 lines        solve_quadratic/cubic/quartic → Complex roots
  concepts.hpp           95 lines        Magma through Field, LieGroup, LieAlgebra
  functions.hpp          80 lines        dot/cross/normalize/lerp/distance, power(), commutator(), adjoint(), poly_eval()
  vec_expr.hpp          219 lines        VecBinExpr, VecScalarExpr — expression templates
  vec_simd.hpp          126 lines        SSE2/AVX2 dispatch
  dual.hpp              129 lines        Dual<T> — forward-mode autodiff, satisfies Scalar
  calculus.hpp          127 lines        gradient()/integrate()/minimize() over plain callables
  ode.hpp                70 lines        General IVP (Cauchy problem): euler_step, rk4_step
  linear_solve.hpp      150 lines        solve_direct()/solve_jacobi() — general N×N Ax=b
  eigen_interop.hpp      98 lines        Eigen interop for Vec/Matrix, opt-in via SPATIUM_EIGEN
  literals.hpp           61 lines        _deg, _pi, _x, _y, _z
  format.hpp             47 lines        std::format for Vec
  verify.hpp             53 lines        verify_matrix_group()
  groups/so3.hpp        119 lines        SO3 — rotation group, Rodrigues formula
  groups/se3.hpp        144 lines        SE3 — rigid body motions

spaces/                                 — concrete mathematical spaces
  euclidean.hpp          94 lines        Euclidean<N,T> — flat space, satisfies all concepts
  sphere.hpp            109 lines        Sphere<N,T> — great-circle geodesics
  hyperbolic.hpp         97 lines        Hyperbolic<N,T> — hyperboloid model
  product.hpp           124 lines        ProductSpace<S1,S2> — Cartesian product
  parametric.hpp        298 lines        ParametricSurface<T> — f(u,v)→R³, tessellate()
  implicit.hpp          566 lines        ImplicitSurface<T> — F(x,y,z)=0, marching_cubes()

geometry/                               — primitives, operations, ray casting
  geometry.hpp            20 lines        Domain umbrella — includes the rest of geometry/
  concepts.hpp           109 lines        Shape, ClosedShape, Measurable, Bounded, BoundedRegion
  line.hpp              171 lines        Line, Ray, Segment
  hyperplane.hpp        104 lines        Hyperplane (Plane = 3D alias)
  triangle.hpp          187 lines        Triangle with subdivision
  polygon.hpp           244 lines        Polygon — ear-clipping, winding number
  box.hpp               152 lines        Box — AABB, contains(), intersects()
  circle.hpp            204 lines        Circle, Disk — containment, distance
  simplex.hpp           236 lines        Simplex<N,K,T> — K-simplex in N-D
  intersection.hpp      369 lines        11 intersect() overloads, Moller-Trumbore, slab method
  distance.hpp          206 lines        distance() for all shape pairs
  clip.hpp              218 lines        clip() point/line/segment against shapes
  boolean.hpp           329 lines        intersection_region(), difference_region(), operator&/-/|
  convex_hull.hpp        59 lines        convex_hull() — Andrew's algorithm (2D)
  ray_surface.hpp       326 lines        Quadric<T>, ray_quadric(), ray_quadric_proximity()
  ray_hit.hpp            74 lines        ray_hit() — concept-driven unified dispatch (RayHittable<S>)
  ray_parametric.hpp    303 lines        Ray x ParametricSurface via Newton UV -- depends on spaces/
  transform.hpp         264 lines        AffineTransform, lazy(), TransformLeaf/TransformExpr
  surface_adapter.hpp   111 lines        ShapeSurface — geometry shape → Surface concept
  make.hpp              113 lines        Factory functions + operator| for intersection
  format.hpp            213 lines        std::format for all geometry types

mesh/                                   — indexed triangle mesh on any Surface
  mesh.hpp                              Mesh<S> — vertices + indices
  subdivision.hpp                       subdivide_once(), subdivide() with projection
  lod.hpp                               LodChain — multi-level LOD
  primitives.hpp                        icosahedron(), tetrahedron(), grid_mesh()
  topology.hpp                          MeshTopology — edges, adjacency, boundary
  geodesic.hpp                          geodesic_distances() (Dijkstra single + multi-source), shortest_path(), GeodesicMethod enum
  geodesic_types.hpp                    DistanceField<S>, GeodesicPath<S>, no_vertex sentinel
  heat_geodesic.hpp                     Heat method (Crane 2013) — HeatSolver<S>, requires SPATIUM_EIGEN
  differential.hpp                      Cotangent Laplacian, mass matrix, face gradients, divergence
  dec.hpp                               Discrete Exterior Calculus — Form0/1/2, exterior_derivative, hodge_star, laplace_beltrami_dec, DecHeatSolver
  transport.hpp                         Parallel transport (Schild's ladder)
  voronoi.hpp                           Geodesic Voronoi — multi-source Dijkstra + face_labels()
  operations.hpp                        merge, transform, flip_normals, compute_normals
  quality.hpp                           aspect_ratio(), min_angle()
  face_metrics.hpp                      Per-face area/normal helpers

spatial/                                — acceleration structures
  bvh.hpp               366 lines        BVH<Shape> — SAH build, ray_cast(), nearest()

render/                                 — CPU-raytracer utilities
  supersample.hpp                       supersample_pixel() — NxN jittered-grid antialiasing for any per-pixel raytracer

discrete/                               — finite set theory
  finite_set.hpp        172 lines        FiniteSet<T> — ∪ ∩ ∖ △ ⊆, power_set, cartesian

physics/                                — atomic physics + classical mechanics
  elements.hpp                          Element DB (118/118) — backed by src/physics/elements.cpp, the one compiled (non-header-only) translation unit in the library, see "Header-only spine" below

physics/atomic/                         — visualization-support models, split out from physics/'s top level per docs/conventions.md's subdivision rule
  atom_model.hpp                        AtomModel — rejection sampling on orbitals
  atom_palette.hpp                      CPK-style colour palette per element
  bohr_model.hpp                        BohrModel — electron shells, transition diagrams
  orbital.hpp                           Radial + spherical-harmonic wavefunctions
  atom_svg.hpp                          SVG rendering of atoms

physics/mechanics/                      — classical/geometric mechanics; see docs/concept-driven-physics.md
  units.hpp                             Compile-time SI Quantity<M,L,T,I,K,N,J> + literals + constants
  state.hpp                             State containers shared by integrators
  body.hpp                              PointMass<N,T>, RigidBody<N,T>
  force.hpp                             UniformGravity, PointGravity, Spring, Damper composables
  integrator.hpp                        euler_step, semi_implicit_euler_step, verlet_step, rk4_step
  symplectic.hpp                        SymplecticManifold concept, CotangentBundle<M>, verify_symplecticity_drift, Yoshida4
  lie_integrator.hpp                    lie_euler_step, lie_midpoint_step, lie_rkmk4_cf_step — commutator-free RKMK on Lie groups
  lgvi.hpp                              Lie-group variational integrator on SO(3) via Cayley 1-cut
  manifold_body.hpp                     PointOnManifold<M> + analytical geodesic_step<Sphere<N>>
  variational.hpp                       DiscreteLagrangian concept, SeparableMidpointLagrangian, variational_step_separable
  continuum.hpp                         StrainEnergy concept, deformation_gradient, right_cauchy_green, green_strain, SaintVenantKirchhoff
  contact.hpp                           ipc_barrier, ipc_contact_energy, ipc_contact_force_on
  narrow_phase.hpp                      ContactSurface concept + point_to(p, surf) overload set (Sphere/Torus/ParametricSurface/generic)
  xpbd.hpp                              XpbdParticle, XpbdDistanceConstraint, xpbd_step, build_distance_constraints

io/                                     — import/export
  table.hpp              30 lines        ASCII table
  svg.hpp               215 lines        SVG 2D writer — wireframe, filled, colored
  obj.hpp                89 lines        OBJ load/save
  stl.hpp               124 lines        STL load/save (binary + ASCII)

viewer/                                 — Vulkan 1.3 real-time viewer
  app.hpp               221 lines        Multi-mesh, point clouds, ImGui, screenshots
  camera.hpp             43 lines        Spherical orbital camera

root/
  point.hpp              62 lines        Point<Space> — type-safe wrapper
  morphism.hpp          113 lines        Morphism<From,To> — pipe composition
  spatium.hpp            51 lines        umbrella header
```

**Total: ~15.3K LOC headers (grows steadily — treat as a snapshot, not a target).**

The listing order above reads top-to-bottom as if it were a dependency
order (foundational domains first, consumer domains last) — it isn't
one, and shouldn't be read as one. The real, current dependency graph is
generated from actual `#include` statements by `scripts/
gen_dependency_graph.py` into `docs/dependency-graph.dot`, checked by CI
so this can't drift silently again. Two things that listing order gets
wrong: `spaces/` and `mesh/` have a genuine mutual dependency
(`spaces/parametric.hpp` and `spaces/implicit.hpp` need `mesh/mesh.hpp`
for tessellation; `mesh/primitives.hpp` needs `spaces/euclidean.hpp` and
friends for the surfaces its generators build on) — neither is really
"below" the other. `discrete/` (billed above as independent finite-set
theory) reaches into both `algebra/` and `geometry/`, so it isn't the
zero-dependency base its position implies either. Regenerate the graph
(`python3 scripts/gen_dependency_graph.py`) whenever this section
changes, rather than hand-editing a claim about dependency order here.

## Data Flow

### Geometry Pipeline

```
Factory (make.hpp)  →  Shape  →  intersect/distance/clip  →  Result<T>
    tri(), ray()         │              │                         │
                         │     operator| pipe syntax              │
                         │              │                    pipe-unwrap
                         v              v                         v
              BoundedRegion    intersect_via_subspace()    next operation
                   │
                   v
             Boolean ops (boolean.hpp)
             operator& / operator- / operator|
```

### Mesh Pipeline

```
Surface concept  →  Mesh<S>  →  subdivide(mesh, space, n)  →  refined Mesh<S>
     │                 │              │
     │          icosahedron()    project each midpoint onto Surface
     │                 │
     v                 v
 MeshTopology    LodChain::build()
     │
     v
 Geodesic algorithms (Dijkstra, Voronoi, transport)
```

### Ray Casting Pipeline

```
Ray  →  BVH::ray_cast(ray)        →  Hit{index, t, point}     (mesh-based)
Ray  →  ray_quadric(ray, quadric)  →  vector<RayHit>           (analytical)
Ray  →  ray_quadric_proximity()    →  RayProximity{miss_dist}  (near-miss)
```

### Viewer Pipeline

```
Mesh<S>  →  App::add_mesh(mesh, color)  →  Vulkan vertex/index buffers
                                                  │
Point cloud  →  App::add_points(pts)    →  VK_PRIMITIVE_TOPOLOGY_POINT_LIST
                                                  │
                                           Rasterization pipeline
                                           (solid + wireframe + points)
                                                  │
                                           ImGui overlay (optional)
                                                  │
                                           Swapchain → Window / Screenshot
```

## Expression Templates

Vec arithmetic uses lazy evaluation (`algebra/vec_expr.hpp`):

```
3.0_x + 2.0_y + 1.0_z
  │       │       │
  v       v       v
Vec3    Vec3    Vec3
  \      |      /
   VecBinExpr<OpAdd, Vec3, Vec3>
         \         /
          VecBinExpr<OpAdd, VecBinExpr, Vec3>
                    │
            (materialized only on assignment to Vec3)
```

Transform chains use similar pattern (`geometry/transform.hpp`):

```
lazy(translate) * lazy(rotate) * lazy(scale)
        │              │             │
  TransformLeaf  TransformLeaf  TransformLeaf
         \            |           /
          TransformExpr<L, R>
                  \       /
              TransformExpr<Outer, Inner>
                      │
               .apply(point)   — sequential, no matrix multiply
               .collapse()     — eager matrix combine for batch
```
