# API Reference

## Concepts

> Concept tables list each concept by name only; jump to the
> definitions in the cited header. Line numbers were dropped on
> purpose — they shift with every refactor and were a recurring
> source of stale documentation.

### Space Concepts (`core/concepts.hpp`)

| Concept | Requires | Adds |
|---------|----------|------|
| `Scalar<T>` | regular, ordered, `+`, `-`, `*`, `/`, `-unary`, `T{0}`, `T{1}` | Scalar type |
| `Set<S>` | `PointType`, `ScalarType`, `dimension`, equality | Base space |
| `TopologicalSpace<S>` | Set + `contains(point) -> bool` | Membership |
| `MetricSpace<S>` | Topological + `distance(p, q) -> Scalar` | Distance |
| `VectorSpace<S>` | Set + `VectorType`, `p + v`, `v + w`, `a * v`, `-v` | Linear ops |
| `NormedSpace<S>` | VectorSpace + MetricSpace + `norm(v)` | Magnitude |
| `InnerProductSpace<S>` | NormedSpace + `inner(u, v)` | Dot product |
| `Complete<S>` | `S::is_complete == true` | Cauchy tag |
| `BanachSpace<S>` | NormedSpace + Complete | |
| `HilbertSpace<S>` | InnerProductSpace + Complete | |
| `EuclideanSpace<S>` | HilbertSpace + finite dimension | R^N |
| `Manifold<S>` | Topological + `TangentVector`, `exp_map`, `log_map` | Geodesics |
| `RiemannianManifold<S>` | Manifold + MetricSpace + `metric_at(p, u, v)` | Curvature |
| `Surface<S>` | Manifold + `project(p)`, `normal(p)` | Embedding |

### Geometry Concepts (`geometry/concepts.hpp`)

| Concept | Requires |
|---------|----------|
| `Shape<S>` | `ScalarType`, `PointType`, `ambient_dimension`, `centroid()` |
| `ClosedShape<S>` | Shape + `contains(point) -> bool` |
| `Measurable<S>` | Shape + `measure() -> Scalar` |
| `Bounded<S>` | Shape + `bounding_box()` |
| `DistanceQueryable<S>` | Shape + `distance(point)`, `project(point)` |
| `BoundedRegion<S>` | Shape + `subspace() -> Result<Subspace>` |

### Algebra Concepts (`algebra/concepts.hpp`)

| Concept | Requires |
|---------|----------|
| `Magma<G>` | `ElementType`, `compose(a, b)` |
| `Semigroup<G>` | Magma (associativity by convention) |
| `Monoid<G>` | Semigroup + `identity()` |
| `Group<G>` | Monoid + `inverse(a)` |
| `AbelianGroup<G>` | Group (commutativity by convention) |
| `Ring<R>` | `add`, `negate`, `zero`, `multiply`, `one` |
| `Field<F>` | Ring + `reciprocal(a)` |
| `LieGroup<G>` | Group + `AlgebraType`, `exp(v)`, `log(a)` |
| `LieAlgebra<A>` | `bracket(x, y)`, `add(x, y)`, `scale(s, x)` |

---

## Types

### `Vec<T, N>` (`algebra/vector.hpp`)

N-dimensional vector. Constexpr arithmetic with expression templates and optional SIMD.

```cpp
Vec3 a{1, 2, 3};
auto v = 3.0_x + 2.0_y + 1.0_z;  // via UDLs (literals.hpp)
```

| Method | Line | Returns | Description |
|--------|------|---------|-------------|
| `operator[](i)` | 63 | `T&` | Element access |
| `operator+/-` | 67-73 | `VecBinExpr` | Lazy arithmetic (expression template) |
| `operator*(scalar)` | 79 | `VecScalarExpr` | Lazy scalar multiply |
| `operator/(scalar)` | 83 | `VecScalarExpr` | Lazy scalar divide |
| `operator-()` | 75 | `VecNegExpr` | Lazy negate |
| `operator+=/−=/∗=/÷=` | 93-156 | `Vec&` | In-place (SIMD-dispatched for N=4) |
| `dot(rhs)` | 161 | `T` | Dot product (SIMD for N=4) |
| `norm_squared()` | 183 | `T` | a² + b² + ... (constexpr) |
| `norm()` | 185 | `T` | √(norm_squared) |
| `normalized()` | 192 | `Vec` | Unit vector |
| `cross(rhs)` | 252 | `Vec` | Cross product (3D only) |
| `lerp(other, t)` | 199 | `Vec` | Linear interpolation |
| `distance_to(other)` | 207 | `T` | Euclidean distance |
| `reflect(normal)` | 212 | `Vec` | v - 2(v·n)n |

Aliases: `Vec2`, `Vec3`, `Vec4`, `Vec2f`, `Vec3f`, `Vec4f` (line 265-271).

### `Matrix<T, R, C>` (`algebra/matrix.hpp`)

R×C matrix. Column-major storage.

| Method | Returns | Description |
|--------|---------|-------------|
| `operator()(r, c)` | `T&` | Element access |
| `operator*` | `Matrix` / `Vec` | Matrix multiply / matrix-vector |
| `transpose()` | `Matrix<T,C,R>` | Transpose |
| `inverse()` | `Result<Matrix>` | Inverse (returns error if singular) |
| `identity()` | `Matrix` | Static: identity matrix |
| `determinant()` | `T` | Determinant (up to 4×4) |

Aliases: `Mat3 = Matrix<double,3,3>`, `Mat4 = Matrix<double,4,4>`.

### `Quaternion<T>` (`algebra/quaternion.hpp`)

Unit quaternion for 3D rotations.

```cpp
auto q = Quat::from_axis_angle(Vec3{0,0,1}, pi/2);
auto rotated = q.rotate(Vec3{1,0,0});  // → (0, 1, 0)
```

| Method | Returns | Description |
|--------|---------|-------------|
| `from_axis_angle(axis, angle)` | `Quaternion` | Static factory |
| `from_matrix(mat)` | `Quaternion` | From 3×3 rotation matrix |
| `to_matrix()` | `Mat3` | Convert to rotation matrix |
| `rotate(vec)` | `Vec3` | Rotate vector |
| `slerp(other, t)` | `Quaternion` | Spherical linear interpolation |
| `conjugate()` | `Quaternion` | q* |
| `inverse()` | `Result<Quaternion>` | q⁻¹ (error if zero norm) |
| `norm()` | `T` | |w² + x² + y² + z²| |

Aliases: `Quat = Quaternion<double>`, `Quatf = Quaternion<float>`.

### `Complex<T>` (`algebra/complex.hpp`)

Complex number with constexpr arithmetic.

```cpp
Complex64 z{3.0, 4.0};        // 3 + 4i
auto w = Complex64::from_polar(5.0, pi/4);
auto sq = sqrt(Complex64{-1});  // → (0, 1)  i.e. i
```

| Method | Line | Returns | Description |
|--------|------|---------|-------------|
| `operator+/-/*/÷` | 25-36 | `Complex` | Arithmetic (constexpr) |
| `conjugate()` | 50 | `Complex` | a - bi |
| `magnitude()` | 52 | `T` | √(a² + b²) |
| `magnitude_sq()` | 51 | `T` | a² + b² (constexpr) |
| `phase()` | 53 | `T` | atan2(im, re) |
| `is_real(eps)` | 54 | `bool` | |im| ≤ eps |
| `from_polar(r, theta)` | 22 | `Complex` | Static: r·e^(iθ) |

Free functions: `sqrt(Complex)` (line 63), `cbrt(Complex)` (line 71).
Aliases: `Complex64`, `Complex32` (line 79-80).

### `Point<S>` (`point.hpp`)

Type-safe point wrapper. Prevents mixing points from different spaces.

```cpp
auto p = pt<E3>(Vec3{1, 2, 3});
auto q = p | morph<E3, E2>([](const Vec3& v) -> Vec2 { return {v[0], v[1]}; });
```

| Method | Line | Constraint | Returns |
|--------|------|-----------|---------|
| `raw()` | — | — | `PointType&` |
| `distance_to(other, space)` | — | `MetricSpace<S>` | `Scalar` |
| `exp(tangent, t, space)` | — | `Manifold<S>` | `Point` |
| `log(other, space)` | — | `Manifold<S>` | `TangentVector` |

### `Morphism<From, To>` (`morphism.hpp`)

Runtime map between spaces.

```cpp
auto f = morph<E3, E2>([](const Vec3& p) -> Vec2 { return {p[0], p[1]}; });
auto q = pt<E3>(Vec3{1, 2, 3}) | f;  // q = P(1, 2)
```

| Operation | Line | Syntax | Description |
|-----------|------|--------|-------------|
| Apply | 24 | `f(point)` | Map raw point |
| Apply typed | 29 | `f(Point<From>)` | Map typed Point |
| Pipe apply | 67 | `point \| f` | Apply morphism to point |
| Pipe compose | 73 | `f \| g` | Pipe order: f then g |
| Math compose | 48 | `g * f` | Math order: f first, then g |
| Result pipe | 80 | `Result<Point> \| f` | Auto-unwrap Result |
| Identity | 100 | `identity<S>()` | No-op morphism |
| With inverse | 90 | `morph<A,B>(fwd, inv)` | Bijection |

---

## Spaces

### `Euclidean<N, T>` (`spaces/euclidean.hpp`)

Flat N-dimensional space. Satisfies **all** concepts through EuclideanSpace + Surface.
Static asserts at line 67-87 verify concept satisfaction for N=1,2,3,100 and T=float.

Aliases: `E1`, `E2`, `E3`, `E4`.

### `Sphere<N, T>` (`spaces/sphere.hpp`)

N-sphere embedded in R^{N+1}. Satisfies RiemannianManifold + Surface + Complete.

```cpp
S2 sphere;
Sphere<2> big{.radius = 5.0};
```

| Method | Line | Description |
|--------|------|-------------|
| `distance(a, b)` | 37 | Great-circle distance |
| `exp_map(p, v, t)` | 55 | Walk along geodesic |
| `log_map(p, q)` | 68 | Tangent vector from p toward q |
| `project(p)` | 80 | Normalize to sphere surface |
| `normal(p)` | 84 | Outward unit normal (= normalized p) |
| `project_tangent(p, v)` | 88 | Remove radial component |

Aliases: `S1`, `S2`, `S3`.

### `Hyperbolic<N, T>` (`spaces/hyperbolic.hpp`)

Hyperboloid model: -x₀² + x₁² + ... + xₙ² = -1, x₀ > 0.

| Method | Line | Description |
|--------|------|-------------|
| `minkowski(a, b)` | 23 | Minkowski inner product (static) |
| `distance(a, b)` | 32 | Hyperbolic distance (arccosh) |
| `exp_map(p, v, t)` | 41 | Geodesic motion |
| `log_map(p, q)` | 56 | Inverse exponential |
| `origin()` | 93 | North pole (1, 0, ..., 0) (static) |

Aliases: `H1`, `H2`, `H3`.

### `ProductSpace<S1, S2>` (`spaces/product.hpp`)

Cartesian product. Product metric: d = √(d₁² + d₂²). Conditional Manifold/Riemannian/Surface methods when both components satisfy them (lines 71, 84, 96, 107).

### `ParametricSurface<T>` (`spaces/parametric.hpp`)

User-defined surface via f(u,v) → R³. Auto-computes metric, normal, exp/log via finite differences and Newton iteration.

```cpp
auto torus = make_torus(2.0, 0.5);            // see spaces/parametric.hpp
auto mesh = tessellate(torus, 32, 16);
```

| Method | Line | Description |
|--------|------|-------------|
| `operator()(u, v)` | 39 | Evaluate f(u,v) |
| `find_params(point)` | 112 | Newton UV iteration (8×8 grid + 5 iterations) |
| `du(u, v)` / `dv(u, v)` | 95-103 | Partial derivatives (finite diff) |
| `normal_at(u, v)` | 105 | Cross product of partials |

Factories (in `spaces/parametric.hpp`): `make_torus`, `make_cylinder`, `make_cone`, `make_mobius`, `parametric(fn, domain)`.

### `ImplicitSurface<T>` (`spaces/implicit.hpp`)

Level set F(x,y,z) = 0. Auto gradient, Newton projection, marching cubes.

```cpp
auto sphere = make_implicit_sphere(3.0);      // see spaces/implicit.hpp
auto mesh = marching_cubes(sphere, 64);
```

| Method | Line | Description |
|--------|------|-------------|
| `operator()(x,y,z)` | 37 | Evaluate F |
| `project(p)` | 52 | Newton iteration to surface (20 iters max) |
| `normal(p)` | 66 | Gradient via central differences |
| `gradient(p)` | 94 | ∇F (private, h = ε·1000) |

Factories (in `spaces/implicit.hpp`): `make_implicit_sphere`, `make_implicit_torus`, `make_gyroid`.

---

## Geometry Primitives

All in `namespace spatium::geometry`. Template: `<std::size_t N, Scalar T = double>`.

### `Box<N, T>` (`geometry/box.hpp`)

| Method | Line | Returns | Description |
|--------|------|---------|-------------|
| `extents()` | 42 | `Vec` | Dimensions |
| `centroid()` | 46 | `Vec` | Center |
| `measure()` | 54 | `T` | Area/volume |
| `contains(point)` | 66 | `bool` | Point inside |
| `contains(box)` | 73 | `bool` | Box inside |
| `intersects(box)` | 80 | `bool` | Overlap test |
| `distance(point)` | 111 | `T` | Clamp-based |
| `project(point)` | 107 | `Vec` | Clamp to box |

### `Line<N,T>`, `Ray<N,T>`, `Segment<N,T>` (`geometry/line.hpp`)

| Type | Line | Domain | Key methods |
|------|------|--------|-------------|
| `Line` | 18 | t ∈ (-∞,+∞) | `at(t)`, `parameter(p)`, `project(p)`, `distance(p)` |
| `Ray` | 56 | t ≥ 0 | Same + `subspace() → Result<Line>` |
| `Segment` | 93 | t ∈ [0,1] | Same + `length()`, `midpoint()`, `bounding_box()` |

Factory functions in `make.hpp`: `line(o, d) → Result<Line>` (39), `ray(o, d) → Result<Ray>` (46), `seg(a, b) → Segment` (25).

### `Hyperplane<N, T>` (`geometry/hyperplane.hpp`)

| Method | Line | Returns |
|--------|------|---------|
| `from_normal_and_point(n, p)` | 28 | `Result<Hyperplane>` (static) |
| `from_points(a, b, c)` | 40 | `Result<Hyperplane>` (3D, static) |
| `signed_distance(p)` | 60 | `T` — positive on normal side |
| `distance(p)` | 64 | `T` — unsigned |
| `project(p)` | 68 | `Vec` |
| `contains(p, eps)` | 72 | `bool` |
| `side(p, eps)` | 76 | `int` — -1, 0, +1 |

Alias: `Plane<T> = Hyperplane<3, T>`. Factory: `plane(n, p)` (53), `plane(a, b, c)` (59).

### `Triangle<N, T>` (`geometry/triangle.hpp`)

| Method | Line | Returns |
|--------|------|---------|
| `area()` | 65 | `T` — cross product (3D) or shoelace (2D) |
| `normal()` | 80 | `Vec` (3D only) |
| `centroid()` | 62 | `Vec` |
| `barycentric(p)` | 97 | `Vec<T,3>` |
| `contains(p, eps)` | 111 | `bool` — barycentric check |
| `distance(p)` | 139 | `T` — via project |
| `project(p)` | 126 | `Vec` — closest point |
| `subdivide()` | 145 | `array<Triangle, 4>` — midpoint split |

Factory: `tri(a, b, c)` (18).

### `Polygon<N, T>` (`geometry/polygon.hpp`)

| Method | Line | Returns |
|--------|------|---------|
| `area()` | 30/42 | `T` — shoelace (2D) / triangle fan (3D+) |
| `perimeter()` | 56 | `T` |
| `centroid()` | 63 | `Vec` |
| `normal()` | 71 | `Vec` (3D) |
| `contains(p)` | 79/100 | `bool` — winding (2D) / triangulation (3D+) |
| `distance(p)` | 116 | `T` — min edge distance |
| `triangulate()` | 141/189 | `vector<Triangle>` — ear-clipping (2D) / fan (3D+) |

Factory: `poly({v0, v1, ...})` (80).
Operators: `poly & poly` → intersection, `poly - poly` → difference, `poly + poly` → union (2D, boolean.hpp). `|` is reserved for intersect/pipe and is not overloaded for `Polygon`.

### `Circle<N, T>`, `Disk<N, T>` (`geometry/circle.hpp`)

Circle: curve only (line 18). Disk: filled region (line 105).

| Method | Type | Line | Returns |
|--------|------|------|---------|
| `circumference()` | Circle | 27 | `T` |
| `distance(p)` | Circle | 35 | `T` — distance to ring |
| `area()` | Disk | 113 | `T` — πr² |
| `contains(p)` | Disk | 120 | `bool` — inside check |

Factories: `circle(c, r, n)` (66), `disk(c, r, n)` (73).

### `Simplex<N, K, T>` (`geometry/simplex.hpp`)

K-simplex in N-D. K=0: point, K=1: segment, K=2: triangle, K=3: tetrahedron.

| Method | Line | Returns |
|--------|------|---------|
| `face(i)` | — | `Simplex<N, K-1, T>` |
| `measure()` | — | `T` — via Gram determinant |
| `barycentric(p)` | — | `Vec<T, K+1>` |
| `contains(p, eps)` | — | `bool` |

### `Quadric<T>` (`geometry/ray_surface.hpp`)

General quadric surface via 4×4 symmetric matrix Q. p^T Q p = 0.

```cpp
auto sphere = Quadric<>::sphere(3.0);          // line 46
auto cyl = Quadric<>::cylinder_z(2.0);          // line 53
auto cone = Quadric<>::cone_z();                // line 60
auto ell = Quadric<>::ellipsoid(5.0, 3.0, 2.0); // line 67
```

| Method | Line | Returns |
|--------|------|---------|
| `operator()(p)` | 23 | `T` — evaluate (positive outside) |
| `normal(p)` | 31 | `Vec3` — gradient, normalized |

---

## Operations

### Intersection (`geometry/intersection.hpp`)

Free functions returning `Result<T>`. Also available via `operator|` (`make.hpp:92-96`).

| Pair | Line | Returns |
|------|------|---------|
| `Line + Hyperplane` | 14 | `Result<Vec>` |
| `Ray + Hyperplane` | 29 | `Result<Vec>` (t ≥ 0 check) |
| `Segment + Hyperplane` | 43 | `Result<Vec>` (t ∈ [0,1]) |
| `Line + Line (2D)` | 56 | `Result<Vec>` |
| `Segment + Segment (2D)` | 72 | `Result<Vec>` |
| `Ray + Triangle (3D)` | 88 | `Result<Vec>` — Moller-Trumbore |
| `Line + Triangle (3D)` | 118 | `Result<Vec>` |
| `Ray + Box` | 145 | `Result<Vec>` — slab method |
| `intersect_parameters(Ray, Box)` | 145 | `Result<pair<T,T>>` — t_entry, t_exit |
| `Hyperplane + Hyperplane (3D)` | 179 | `Result<Line>` |
| `Triangle + Triangle (3D)` | 225 | `Result<Segment>` |
| BoundedRegion + BoundedRegion | 277 | `Result<Vec>` — generic via subspace |

### Distance (`geometry/distance.hpp`)

All return `T`. Point-to-shape delegates to member `.distance()`:

| Pair | Line |
|------|------|
| `Vec + Line/Ray/Segment/Hyperplane/Triangle/Box` | 18-46 |
| `Vec + Circle/Disk/Polygon` | 110-120 |
| `Line + Line (3D)` | 50 |
| `Segment + Segment` | 75 |
| `Box + Box` | 125 |
| `Circle + Circle (2D)` | 133 |
| `Triangle + Triangle` | 138 |
| `Polygon + Polygon` | 154 |
| `Segment + Triangle` | 170 |

### Clip (`geometry/clip.hpp`)

`clip(entity, shape) → Result<entity>` — pass through if inside, error if outside.

Point clips: Triangle (22), Segment (28), Disk (43), Box (49), Polygon (55).
Line clips: Triangle (65), Segment (123), Disk (138).
Segment clips: Triangle (167), Disk (194).

### Boolean Operations (`geometry/boolean.hpp`)

| Function/Operator | Line | Returns | Constraint |
|-------------------|------|---------|-----------|
| `intersection_region(A, B)` | 142 | `Result<Polygon>` | `requires Shape<A> && Shape<B>` |
| `difference_region(A, B)` | 219 | `Result<Polygon>` | Polygon only |
| `difference_area(A, B)` | 193 | `Result<T>` | `requires Measurable<A> && Measurable<B>` |
| `symmetric_difference_area(A, B)` | 206 | `Result<T>` | `requires Measurable<A> && Measurable<B>` |
| `poly & poly` | 287 | `Result<Polygon>` | Intersection |
| `poly - poly` | 293 | `Result<Polygon>` | Difference (2D) |
| `poly + poly` (2D) | 299 | `Result<Polygon>` | Union via convex hull |
| `union_of(poly, poly)` | 312 | `Result<Polygon>` | Named alternative for `+` |

Algorithm: Sutherland-Hodgman clipping (line 67-131). Convex polygons only.

### Ray-Quadric (`geometry/ray_surface.hpp`)

Substitutes p(t) = o + t·d into quadric equation → at² + bt + c = 0 → `solve_quadratic`.

| Function | Line | Returns |
|----------|------|---------|
| `ray_quadric(ray, quadric)` | 100 | `vector<RayHit<T>>` — hits sorted by t |
| `ray_quadric_proximity(ray, quadric)` | 126 | `Result<RayProximity<T>>` — miss distance |
| `ray_quadric_full(ray, quadric)` | 155 | `variant<hits, proximity>` |

**RayHit** (line 82): `{t, point, normal}`.
**RayProximity** (line 89): `{closest_t, miss_distance, closest_point}`.
`miss_distance` = |imaginary part| of complex root — geometrically proportional to how far the ray misses the surface.

### Polynomial Solvers (`algebra/polynomial.hpp`)

| Function | Line | Returns | Algorithm |
|----------|------|---------|-----------|
| `solve_quadratic(a, b, c)` | 18 | `array<Complex<T>, 2>` | Discriminant |
| `solve_cubic(a, b, c, d)` | 48 | `array<Complex<T>, 3>` | Cardano |
| `solve_quartic(a, b, c, d, e)` | 92 | `array<Complex<T>, 4>` | Ferrari |
| `real_roots_quadratic(a, b, c)` | 33 | `vector<T>` | Filter |im| < eps |
| `real_roots_cubic(...)` | 148 | `vector<T>` | Filter |
| `real_roots_quartic(...)` | 155 | `vector<T>` | Filter |

### Generic Algebra (`algebra/functions.hpp`)

| Function | Line | Constraint | Description |
|----------|------|-----------|-------------|
| `algebra::power(g, a, n)` | 33 | `Group<G>` | a^n via square-and-multiply |
| `algebra::commutator(g, a, b)` | 43 | `Group<G>` | a·b·a⁻¹·b⁻¹ |
| `algebra::adjoint(g, elem, v)` | 51 | `LieGroup<G>` | Ad_g(v) (numerical) |
| `algebra::poly_eval(ring, coeffs, x)` | 60 | `Ring<R>` | Horner's method |

### Transform (`geometry/transform.hpp`)

**AffineTransform<N, T>** (line 17): (N+1)×(N+1) homogeneous matrix.

| Method/Factory | Line | Description |
|----------------|------|-------------|
| `operator()(point)` | 24 | Apply to Vec |
| `operator*(other)` | 40 | Eager composition |
| `inverse()` | 45 | `Result<AffineTransform>` |
| `translation(offset)` | 98 | Static factory |
| `scaling(factors)` | 104 | Static factory |
| `rotation(axis, angle)` | 117 | Rodrigues (3D) |
| `rotation(angle)` | 138 | 2D rotation |

Shorthand factories (line 174-184): `translate()`, `scale()`, `rotate_x/y/z()`.

**Lazy chains** (line 191+):

```cpp
auto chain = lazy(translate(1,0,0)) * lazy(rotate_z(0.5)) * lazy(scale(2.0));
chain.apply(point);     // sequential, no matrix multiply
chain.collapse();       // eager matrix combine
pt | chain;             // pipe syntax
```

**Pipe-unwrap** (line 159-171): `Result<Vec> | Transform` and `Result<Point> | Transform` auto-unwrap.

---

## Mesh (`mesh/`)

### `Mesh<S>` (`mesh/mesh.hpp`)

Indexed triangle mesh for any `Surface S`.

```cpp
S2 sphere;
auto m = mesh::icosahedron(sphere);
auto refined = mesh::subdivide(m, sphere, 3);
```

| Method | Line | Returns |
|--------|------|---------|
| `vertex_count()` | — | `size_t` |
| `face_count()` | — | `size_t` |
| `triangle(i)` | — | `array<PointType, 3>` |
| `area(space)` | — | `Scalar` — total surface area |

### Subdivision (`mesh/subdivision.hpp`)

`subdivide_once(mesh, space)` — midpoint split + project onto Surface (line 18).
`subdivide(mesh, space, n)` — n iterations (line 55).

### MeshTopology (`mesh/topology.hpp`)

Build from Mesh via `MeshTopology(mesh)` or `MeshTopology(shared_ptr<Mesh>)`.
Edges, vertex neighbors, face-edge map, boundary detection.

### Geodesics (`mesh/geodesic.hpp`)

Dijkstra with `space.distance()` as edge weights — works on ANY Surface.

| Function | Description |
|----------|-------------|
| `geodesic_distances(topo, space, source)` | Single-source Dijkstra |
| `geodesic_distances(topo, space, sources)` | Multi-source Dijkstra (`std::span<const uint32_t>`) |
| `geodesic_distances(topo, space, source, GeodesicMethod)` | Method selector — `Dijkstra` or `Heat` (Heat requires `SPATIUM_EIGEN=ON`) |
| `shortest_path(topo, space, src, dst)` | Dijkstra + predecessor trace → `GeodesicPath` |

### Voronoi (`mesh/voronoi.hpp`)

`geodesic_voronoi(mesh, space, seeds)` — multi-source Dijkstra + label propagation.
Returns `VoronoiResult` with `face_labels` for per-face viewer coloring.

### Transport (`mesh/transport.hpp`)

`parallel_transport(mesh, space, path, v)` — Schild's ladder (exp/log only, any Manifold).

### Operations (`mesh/operations.hpp`)

`merge()`, `transform()`, `flip_normals()`, `compute_face_normals()`, `compute_vertex_normals()`, `centered()`.

---

## Spatial (`spatial/bvh.hpp`)

### `BVH<Shape>`

SAH-built bounding volume hierarchy. Stack-based traversal.

```cpp
auto bvh = BVH<Triangle3>::build(triangles);
auto hit = bvh.ray_cast(ray);   // optional<Hit>
```

| Method | Line | Returns | Description |
|--------|------|---------|-------------|
| `ray_cast(ray)` | 63 | `optional<Hit>` | Closest intersection |
| `ray_test(ray)` | 120 | `bool` | Any intersection (early exit) |
| `nearest(point)` | — | `optional<NearestResult>` | Closest primitive |
| `query_box(box)` | — | `vector<size_t>` | Primitives in box |

**Hit**: `{index, t, point}`. SAH build at construction (12 bins).

---

## IO (`io/`)

| Module | File | Functions |
|--------|------|-----------|
| Table | `io/table.hpp` | `Table("Col1", "Col2").row(a, b).print()` |
| SVG | `io/svg.hpp` | `mesh_wireframe()`, `mesh_filled()`, `mesh_colored()`, `mesh_to_svg()` |
| OBJ | `io/obj.hpp` | `load_obj(path) → Result<Mesh<E3>>`, `save_obj(mesh, path)` |
| STL | `io/stl.hpp` | `load_stl(path) → Result<Mesh<E3>>`, `save_stl(mesh, path)` |

---

## Discrete (`discrete/finite_set.hpp`)

### `FiniteSet<T>`

| Operator | Line | Description |
|----------|------|-------------|
| `\|` | 35 | Union (∪) |
| `&` | 44 | Intersection (∩) |
| `-` | 53 | Difference (∖) |
| `^` | 62 | Symmetric difference (△) |
| `<=` | 71 | Subset (⊆) |
| `<` | 77 | Proper subset (⊂) |

Methods: `size()`, `empty()`, `contains()`, `insert()`, `erase()`, `power_set()`, `cartesian()`, `map()`, `filter()`.

---

## Verification (`core/verify.hpp`)

| Function | Line | Constraint | Checks |
|----------|------|-----------|--------|
| `verify_metric` | 38 | `MetricSpace<S>` | non-negativity, symmetry, triangle inequality |
| `verify_inner_product` | 82 | `InnerProductSpace<S>` | symmetry, positive-definiteness, linearity |
| `verify_exp_log` | 115 | `Manifold<S>` | exp(p, log(p, q), 1) ≈ q |
| `verify_norm_consistency` | 140 | `InnerProductSpace<S>` | norm(v) = √(inner(v, v)) |

---

## Literals (`algebra/literals.hpp`)

```cpp
using namespace spatium::literals;
```

| Literal | Line | Result |
|---------|------|--------|
| `45_deg` | 7 | `double` — radians (π/4) |
| `2_pi` | 15 | `double` — 2π |
| `3.0_x` | 24 | `Vec3{3, 0, 0}` |
| `2.0_y` | 26 | `Vec3{0, 2, 0}` |
| `1.0_z` | 28 | `Vec3{0, 0, 1}` |

Combination: `3.0_x + 2.0_y + 1.0_z` → `Vec3{3, 2, 1}` (via expression templates).

---

## Precision (`core/precision.hpp`)

```cpp
Euclidean<3, Real50> space;
Vec<Real50, 3> a{Real50{0}, Real50{0}, Real50{0}};
```

All types work: `Vec<Real50, 3>`, `Triangle<3, Real50>`, `Morphism<Euclidean<3, Real50>, ...>`.
ADL-friendly: `using std::sqrt; sqrt(x)` finds correct overload.

---

## Physics — Mechanics (`physics/mechanics/*.hpp`)

### Units (`units.hpp`)

Compile-time SI dimensional analysis via `std::ratio` exponents.

```cpp
template<int M, int L, int T, int I, int K, int N, int J, Scalar S = double>
struct Quantity;
```

Derived aliases: `Mass`, `Length`, `Time`, `Velocity`, `Acceleration`, `Force`, `Energy`, `Power`, `Momentum`. Literals (`_kg`, `_m`, `_s`, `_N`, `_J`, `_W`). Constants: `g_earth`, `G_newton`, `c_light`, `h_planck`, `k_boltz`.

### Bodies, forces, integrators (`body.hpp`, `force.hpp`, `integrator.hpp`)

```cpp
template<std::size_t N, Scalar T = double> struct PointMass { mass, x, v; };
struct UniformGravity; struct PointGravity; struct Spring; struct Damper;
void euler_step        (body, force, dt, t=0);
void semi_implicit_euler_step(body, force, dt, t=0);
void verlet_step       (body, force, dt, t=0);           // O(dt²) symplectic
void yoshida4_step     (body, force, dt, t=0);           // O(dt⁴) composition
void rk4_step          (body, force, dt, t=0);           // O(dt⁴), not symplectic
```

### Symplectic manifolds (`symplectic.hpp`)

```cpp
template<typename S> concept SymplecticManifold = ...;
template<Manifold M> struct CotangentBundle;       // T*M with ω = dq ∧ dp
template<SymplecticManifold S>
double verify_symplecticity_drift(state, step_fn, eps, dt);
```

### Lie-group integrators (`lie_integrator.hpp`)

```cpp
auto lie_euler_step   (state, rhs, h);                 // O(h)
auto lie_midpoint_step(state, rhs, h);                 // O(h²)
auto lie_rkmk4_cf_step(state, rhs, h);                 // commutator-free 4-th order
```

### Manifold-embedded mechanics (`manifold_body.hpp`)

```cpp
template<Manifold M> struct PointOnManifold { M space; point; velocity; };
PointOnManifold<Sphere<N,T>> geodesic_step(body, dt);    // exact great-circle
```

### Discrete Exterior Calculus (`mesh/dec.hpp`)

```cpp
template<int K, Surface S> struct DiscreteForm { coeffs; };
using Form0 = DiscreteForm<0>;     // vertex-valued
using Form1 = DiscreteForm<1>;     // edge-valued (one-form)
using Form2 = DiscreteForm<2>;     // face-valued (two-form)

Form1 exterior_derivative_0(Form0, topo);
Form2 exterior_derivative_1(Form1, topo);
Form2 hodge_star_0(Form0, topo);   // ⋆₀ = mass matrix
Form1 hodge_star_1(Form1, topo);   // ⋆₁ = cotangent weights
Form0 hodge_star_2(Form2, topo);

SparseMatrix laplace_beltrami_dec(topo);  // Δ = δd + dδ for 0-forms
DecHeatSolver make_dec_heat_solver(topo, dt);  // backward-Euler (M + dt·L)
```

### Variational integrators (`variational.hpp`)

```cpp
template<typename LD, std::size_t N, typename T>
concept DiscreteLagrangian = requires (LD ld,
                                       const Vec<T,N>& qk,
                                       const Vec<T,N>& qk1,
                                       T dt) {
    { ld.action(qk, qk1, dt)   } -> std::convertible_to<T>;
    { ld.dL_d_dq1(qk, qk1, dt) } -> std::convertible_to<Vec<T,N>>;
    { ld.dL_d_dq2(qk, qk1, dt) } -> std::convertible_to<Vec<T,N>>;
};

template<std::size_t N, typename T,
         typename Kinetic, typename Potential, typename GradPotential>
struct SeparableMidpointLagrangian;                  // L(q, q̇) = T(q̇) − V(q)

template<std::size_t N, typename T,
         typename Kinetic, typename Potential, typename GradPotential>
auto make_separable_midpoint_lagrangian(K, V, gradV, mass = T{1});

template<std::size_t N, Scalar T, typename GradPotential>
void variational_step_separable(PointMass<N,T>& body,
                                GradPotential&& grad_V, T dt);   // → verlet_step
```

### LGVI — Lie-Group Variational Integrator (`lgvi.hpp`)

```cpp
struct LGVIRigidBodyState { SO3::Element R; Vec3 Pi; };

LGVIRigidBodyState lgvi_rigid_body_step(state, J_diag, h);  // implicit Cayley
Vec3   lgvi_spatial_angular_momentum(state);
double lgvi_kinetic_energy(Pi, J_diag);
```

Exact orientation (`R^T R = I`), exact discrete-Noether (spatial momentum), exact Casimir (`|Π|²`). Energy drift ~3 % over 5 000 steps at `h = 1e-3` — Cayley 1-cut limitation; full Lee-Leok-McClamroch Newton on joint `(F, Π)` is a documented follow-up.

### Geometric continuum (`continuum.hpp`)

```cpp
template<typename MFrom, typename MTo, typename Phi, Scalar T = double>
struct DeformationMap { Phi phi; ...; };

template<typename MFrom, typename MTo, typename F, Scalar T = double>
auto make_deformation_map(F&& phi);                      // factory + CTAD

Matrix<T,N,N> deformation_gradient(DeformationMap, q);   // F = Dφ via FD
Matrix<T,N,N> right_cauchy_green(F);                     // C = F^T F
Matrix<T,N,N> green_strain(F);                           // E = ½(C − I)

template<typename W, Scalar T, std::size_t N>
concept StrainEnergy = ...;

template<Scalar T, std::size_t N>
struct SaintVenantKirchhoff;                             // W = ½λ(tr E)² + μ tr(E²)
```

### IPC barrier kernel (`contact.hpp`)

```cpp
T  ipc_barrier         (T d, T d_hat);     // = -(d-d̂)²·log(d/d̂) for 0 < d < d̂
T  ipc_barrier_grad    (T d, T d_hat);     // ≤ 0 in active band
T  ipc_barrier_hessian (T d, T d_hat);
T  ipc_default_stiffness(T d_hat);
T  ipc_potential       (T d, T d_hat);     // κ·barrier
```

C² smooth across the threshold `d̂`. Diverges as `d → 0⁺` — provably no-penetration when the barrier is integrated implicitly (Li-Kaufman 2020, SIGGRAPH).

### XPBD (`xpbd.hpp`)

```cpp
template<std::size_t N, Scalar T = double>
struct XpbdParticle            { Vec x, x_prev; T w; };  // w = 1/mass

template<std::size_t N, Scalar T = double>
struct XpbdDistanceConstraint  { i, j; T rest, compliance, lambda; };

template<std::size_t N, Scalar T = double>
using  XpbdBendingDistanceConstraint = XpbdDistanceConstraint<N, T>;

void xpbd_solve_distance(c, particles, dt);
void xpbd_step(particles, c_begin, c_end, apply_accel, dt, n_iter = 4);
std::vector<XpbdDistanceConstraint> build_distance_constraints(
    particles, faces, compliance = 0);
std::vector<XpbdBendingDistanceConstraint> build_bending_distance_constraints(
    particles, faces, compliance = 0);
```

### Narrow-phase + ContactSurface (`narrow_phase.hpp`)

```cpp
template<Scalar T = double>
struct ContactQuery {
    T distance;                 // ≥ 0, unsigned distance
    Vec<T, 3> closest_point;
    Vec<T, 3> normal;           // outward unit
    bool inside;
    T signed_distance() const;
};

// Closed-form helpers (3.7 ns / 17.7 ns):
ContactQuery<T> point_to_sphere(p, centre, radius);
ContactQuery<T> point_to_torus (p, const geometry::Torus<T>&);

// Concept-driven unified API. `point_to(p, surf)` resolves
// statically to the closed-form path for Sphere<2> / Torus, and to
// the generic Surface fallback (Newton-on-UV via project + normal,
// ~µs) for ParametricSurface, ImplicitSurface, and any user type
// satisfying the Surface concept.
template<typename S, typename T = double>
concept ContactSurface = requires(const S& s, const Vec<T, 3>& p) {
    { point_to(p, s) } -> std::same_as<ContactQuery<T>>;
};

ContactQuery<T> point_to(p, Sphere<2,T>);                // forwards
ContactQuery<T> point_to(p, geometry::Torus<T>);
ContactQuery<T> point_to(p, ParametricSurface<T>);
template<Surface S> ContactQuery<...> point_to(p, S surf);  // fallback

// IPC on top of a ContactQuery:
T        ipc_contact_energy  (ContactQuery, d_hat, kappa = 1);
Vec<T,3> ipc_contact_force   (ContactQuery, d_hat, kappa = 1);

template<typename S, Scalar T>
    requires ContactSurface<S, T>
Vec<T,3> ipc_contact_force_on(p, surface, d_hat, kappa = 1);
```

Benchmarks (i5-1235U, `-O3`): analytical paths 3.7–17.7 ns/op, parametric 5.2 µs/op, full IPC pipeline 12.6–25.5 ns/op. Order-of-magnitude faster than typical GJK (~100-500 ns/op) on the analytical shapes.

See `docs/concept-driven-physics.md` for the full concept-hierarchy specification and dispatch rationale.

## Viewer (`viewer/app.hpp`)

Vulkan 1.3 real-time viewer with multi-mesh, point clouds, ImGui.

| Method | Description |
|--------|-------------|
| `add_mesh(MeshData, color)` | Add indexed triangle mesh |
| `add_mesh(Mesh<S>, color)` | Template: auto-project to MeshData |
| `add_points(data, color, size)` | Add point cloud |
| `update_mesh(idx, data)` | Replace mesh data |
| `update_mesh_vertices(idx, data)` | Update vertices only (animation) |
| `fit_camera(radius)` | Auto-scale camera |
| `save_screenshot(path)` | Export PNG |
| `set_key_callback(fn)` | Key handler |
| `set_frame_callback(fn)` | Per-frame callback |
| `run()` | Main loop |

Rendering pipelines: solid (per-face HSV coloring), wireframe (depth-biased), points.
ImGui integration: `SPATIUM_HAS_IMGUI` compile flag.
