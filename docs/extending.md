# Extending Spatium

## Defining a Custom Space

To create a new space, define a struct with the required type aliases and methods. C++ concepts check the interface at compile time — no inheritance or registration needed.

### Minimal: MetricSpace

A metric space needs distance between points. See the `MetricSpace` concept in `core/concepts.hpp` for requirements.

```cpp
#include <spatium/core/concepts.hpp>
#include <spatium/algebra/vector.hpp>

using namespace spatium;

struct ManhattanPlane {
    using ScalarType = double;
    using PointType = Vec<double, 2>;
    static constexpr std::size_t dimension = 2;

    bool contains(const PointType&) const { return true; }

    ScalarType distance(const PointType& a, const PointType& b) const {
        using std::abs;
        return abs(a[0] - b[0]) + abs(a[1] - b[1]);
    }
};

static_assert(MetricSpace<ManhattanPlane>);
```

This gives you: `Point<ManhattanPlane>`, `pt<ManhattanPlane>(...)`, `verify_metric()`.

### Full: RiemannianManifold + Surface

For mesh subdivision, geodesics, and the viewer to work, satisfy `Surface` (in `core/concepts.hpp`):

```cpp
struct FlatTorus {
    using ScalarType = double;
    using PointType = Vec<double, 3>;       // embedded in R³
    using TangentVector = Vec<double, 3>;
    static constexpr std::size_t dimension = 2;
    static constexpr bool is_complete = true;

    double major_radius = 2.0;
    double minor_radius = 0.5;

    // TopologicalSpace
    bool contains(const PointType& p) const {
        using std::sqrt; using std::abs;
        auto r = sqrt(p[0]*p[0] + p[1]*p[1]);
        auto d = (r - major_radius)*(r - major_radius) + p[2]*p[2];
        return abs(d - minor_radius*minor_radius) < 1e-8;
    }

    // MetricSpace
    ScalarType distance(const PointType& a, const PointType& b) const {
        return (a - b).norm();  // ambient approximation
    }

    // Manifold
    PointType exp_map(const PointType& p, const TangentVector& v, ScalarType t) const {
        return project(p + v * t);
    }
    TangentVector log_map(const PointType& p, const PointType& q) const {
        auto diff = q - p;
        auto n = normal(p);
        return diff - n * n.dot(diff);
    }

    // RiemannianManifold
    ScalarType metric_at(const PointType&, const TangentVector& u, const TangentVector& v) const {
        return u.dot(v);
    }

    // Surface
    PointType project(const PointType& p) const {
        using std::sqrt; using std::atan2; using std::cos; using std::sin;
        auto theta = atan2(p[1], p[0]);
        auto r = sqrt(p[0]*p[0] + p[1]*p[1]);
        auto phi = atan2(p[2], r - major_radius);
        return PointType{
            (major_radius + minor_radius * cos(phi)) * cos(theta),
            (major_radius + minor_radius * cos(phi)) * sin(theta),
            minor_radius * sin(phi)
        };
    }
    TangentVector normal(const PointType& p) const {
        using std::sqrt; using std::atan2; using std::cos; using std::sin;
        auto theta = atan2(p[1], p[0]);
        auto r = sqrt(p[0]*p[0] + p[1]*p[1]);
        auto phi = atan2(p[2], r - major_radius);
        return TangentVector{cos(phi)*cos(theta), cos(phi)*sin(theta), sin(phi)};
    }
};

static_assert(RiemannianManifold<FlatTorus>);
static_assert(Surface<FlatTorus>);
static_assert(Complete<FlatTorus>);
```

This gives you:
- `Mesh<FlatTorus>`, `subdivide_once()`, `subdivide()`, `LodChain<FlatTorus>`
- `geodesic_distances()`, `shortest_path()`, `geodesic_voronoi()`
- `parallel_transport()`
- `verify_metric()`, `verify_exp_log()`
- `App::add_mesh(mesh, color)` in the Vulkan viewer

### Alternatively: Use ParametricSurface

Instead of implementing all methods manually, use `ParametricSurface<T>` (in `spaces/parametric.hpp`):

```cpp
auto my_surface = parametric(
    [](double u, double v) -> Vec3 {
        return Vec3{u * std::cos(v), u * std::sin(v), u * u};
    },
    Domain{0.0, 2.0, 0.0, 2 * std::numbers::pi},
    true,   // periodic in u?
    true    // periodic in v?
);

auto mesh = tessellate(my_surface, 32, 32);
```

ParametricSurface auto-computes: `project()`, `normal()`, `exp_map()`, `log_map()`, `metric_at()`, `distance()` via finite differences and Newton iteration (`find_params()` in `spaces/parametric.hpp`).

### Letting your space into the scene DSL: `chart_of`

Satisfying `RiemannianManifold` or `Surface` is not what the scene DSL
(`io::build`) requires. Its operations — `offset_surface`,
`sample_surface_uniform`, tessellation — call `evaluate`, `normal_at`,
`area_element`, `domain` and `periodic_u`/`periodic_v`, and none of them
calls `project`, `normal` or `exp_map`. That requirement is the `Chart`
concept, in `spaces/chart.hpp`. It is not part of the
`Set → Manifold → Surface` hierarchy on purpose: a chart is a choice of
coordinates *on* a space, and a space can carry many or none.

A space joins the DSL by having a `chart_of` overload written beside it,
in its own namespace. Nothing in Spatium is edited, and there is no
registry — the call is resolved by ADL at compile time:

```cpp
namespace my_lib {

struct Ribbon {
    using ScalarType = double;   // Chart reads the scalar from here
    double radius = 1.0, height = 2.0;
};

spatium::ParametricSurface<double> chart_of(const Ribbon& r) {
    return spatium::ParametricSurface<double>(
        [r](double u, double v) -> spatium::Vec<double, 3> {
            return {r.radius * std::cos(u), r.radius * std::sin(u), v};
        },
        {0.0, 2 * std::numbers::pi, 0.0, r.height},
        /*periodic_u=*/true, /*periodic_v=*/false);
}

} // namespace my_lib

spatium::io::build::Trace<double> trace;
auto base  = trace.space(my_lib::Ribbon{.radius = 2.0, .height = 1.0});
auto shell = trace.offset_shell(base, 0.25, EdgeRule::ZeroThickness);
auto dots  = trace.scatter(trace.cube(Vec3{0.05, 0.05, 0.05}), base, 32);
```

Every DSL operation works on it, because entering as a chart is entering
as the thing those operations were always written against.

`Sphere<2, T>` ships with a `chart_of` — the standard θ/φ chart, which is
its first parametrization rather than a generalization of an existing
one. (`Sphere<N, T>` is the N-sphere in R^{N+1}, so the ordinary sphere
in R³ is `Sphere<2, T>`; `Sphere<3, T>` lives in R⁴ and has no chart of
this shape, and `Chartable<Sphere<3, double>, double>` is `false`.)

Two things to know when writing one:

- **Where your chart degenerates is your business, and it is normal.**
  The sphere chart's `v` edges collapse to the poles, where
  `area_element` (which is `r² sin v`) reaches zero and the
  finite-difference normal goes with it. `periodic_v()` is `false` there
  and the surface is still closed — a flag describes the chart, not the
  topology, which is why `is_closed()` tests geometry instead.
- **A chart is not an exact analytic form.** Giving a space a chart does
  not make a renderer able to hit it without triangles; that is the
  separate `exact` slot on a trace node, and a node with a chart but no
  exact form tessellates. See `RenderLevel` in `io/build.hpp`.

### Verification

After defining your space, verify axioms:

```cpp
FlatTorus torus;
std::array samples = { /* points on the torus */ };

auto metric_ok = verify_metric(torus, std::span{samples});       // core/verify.hpp
auto explog_ok = verify_exp_log(torus, std::span{samples}, 1e-4); // core/verify.hpp
```

---

## Defining a Custom Geometry Primitive

Geometry primitives satisfy concepts from `geometry/concepts.hpp`.

### Minimal: Shape (in `geometry/concepts.hpp`)

```cpp
struct Ellipse2D {
    using ScalarType = double;
    using PointType = Vec<double, 2>;
    static constexpr std::size_t ambient_dimension = 2;

    PointType center;
    double semi_major, semi_minor;

    PointType centroid() const { return center; }
};

static_assert(geometry::Shape<Ellipse2D>);
```

### Full: all geometry concepts

```cpp
struct Ellipse2D {
    // ... (as above) ...

    bool contains(const PointType& p) const { /* x²/a² + y²/b² ≤ 1 */ }
    double measure() const { return std::numbers::pi * semi_major * semi_minor; }
    geometry::Box<2> bounding_box() const { /* ... */ }
    double distance(const PointType& p) const { /* ... */ }
    PointType project(const PointType& p) const { /* ... */ }
};

static_assert(geometry::ClosedShape<Ellipse2D>);     // contains
static_assert(geometry::Measurable<Ellipse2D>);       // measure
static_assert(geometry::Bounded<Ellipse2D>);           // bounding_box
static_assert(geometry::DistanceQueryable<Ellipse2D>); // distance + project
```

With `Bounded`: works with `BVH<Ellipse2D>` for spatial queries.
With `Measurable`: works with `difference_area()`, `symmetric_difference_area()`.

**Naming convention for `measure()`.**  The `Measurable` concept is
deliberately dimension-agnostic — `measure()` is whatever the natural
Hausdorff k-measure is for the shape's intrinsic dimension k (length for
`k = 1`, area for `k = 2`, volume for `k = 3`, …).  This is what lets a
single boolean-ops or integration kernel work uniformly across `Segment`,
`Triangle`, `Polygon`, `Box`, `Simplex<N, K>` for any `K`.

For the common cases the built-in shapes additionally expose a friendlier
synonym so user code reads as it would on paper:

| Shape                              | k | Idiomatic alias        |
|------------------------------------|:-:|------------------------|
| `Segment`, `Line`                  | 1 | `length()`             |
| `Triangle`, `Polygon`, `Disk`, `Circle`, 2-D `Box` | 2 | `area()`               |
| 3-D `Box`, eventually 3-D simplex  | 3 | `volume()` *(planned)* |
| `Simplex<N, K>` for general `K`    | k | `measure()` only       |

When you add a new shape, give it `measure()` first (so the concept is
satisfied and your shape composes with the generic algorithms) and then,
if the dimension is one of the table rows above, add the matching alias as
a one-liner that **forwards** to `measure()` rather than duplicating the
formula.  That way the concept-side and the user-facing names are
guaranteed to stay in agreement — the only place the math lives is inside
`measure()`.

```cpp
double measure() const { /* the actual formula */ }
double area()    const { return measure(); }   // 2-D shape: alias only
```

---

## Adding Intersection Overloads

Define a free function in `spatium::geometry`. The `operator|` in `geometry/make.hpp` auto-dispatches to `intersect()` via a requires clause:

```cpp
namespace spatium::geometry {

template<Scalar T>
Result<Vec<T, 2>> intersect(const Ray<2, T>& ray, const Ellipse2D& ellipse) {
    // Substitute ray into ellipse equation → quadratic
    // Use solve_quadratic() from algebra/polynomial.hpp
    auto roots = solve_quadratic(a, b, c);
    // ... filter real positive roots ...
}

} // namespace spatium::geometry
```

Now `ray | ellipse` works automatically.

---

## Adding Distance Overloads

Add free functions in `spatium::geometry`:

```cpp
namespace spatium::geometry {

template<Scalar T>
T distance(const Vec<T, 2>& p, const Ellipse2D& e) {
    return e.distance(p);
}

template<Scalar T>
T distance(const Ellipse2D& a, const Ellipse2D& b) {
    // ...
}

} // namespace spatium::geometry
```

---

## Defining a Custom Algebraic Structure

Satisfy concepts from `algebra/concepts.hpp`.

### Group (in `algebra/concepts.hpp`)

```cpp
struct ZnGroup {
    using ElementType = int;
    int n;

    ElementType identity() const { return 0; }
    ElementType compose(const ElementType& a, const ElementType& b) const {
        return (a + b) % n;
    }
    ElementType inverse(const ElementType& a) const {
        return (n - a) % n;
    }
};

static_assert(algebra::Group<ZnGroup>);
```

This gives you: `algebra::power(zn, elem, k)`, `algebra::commutator(zn, a, b)`.

### LieGroup (in `algebra/concepts.hpp`)

Add `AlgebraType`, `exp(v)`, `log(a)`:

```cpp
struct MyLieGroup {
    using ElementType = Mat3;
    using AlgebraType = Vec3;
    // ... Group methods ...
    ElementType exp(const AlgebraType& v) const { /* ... */ }
    AlgebraType log(const ElementType& a) const { /* ... */ }
};

static_assert(algebra::LieGroup<MyLieGroup>);
```

This gives you: `algebra::adjoint(g, elem, v)`.

### Ring (in `algebra/concepts.hpp`)

```cpp
struct DoubleRing {
    using ElementType = double;
    ElementType add(const ElementType& a, const ElementType& b) const { return a + b; }
    ElementType negate(const ElementType& a) const { return -a; }
    ElementType zero() const { return 0.0; }
    ElementType multiply(const ElementType& a, const ElementType& b) const { return a * b; }
    ElementType one() const { return 1.0; }
};

static_assert(algebra::Ring<DoubleRing>);
```

This gives you: `algebra::poly_eval(ring, coeffs, x)` — Horner's method polynomial evaluation.

---

## Using Custom Scalar Types

Any type satisfying `Scalar` (in `core/concepts.hpp`) works with all Spatium types:

- `std::regular<T>` — default constructible, copyable, equality comparable
- `std::totally_ordered<T>` — all comparison operators
- Arithmetic: `+`, `-`, `*`, `/`, unary `-`
- Literals: `T{0}`, `T{1}`

For math functions, ensure ADL finds them:

```cpp
struct MyFloat {
    // ... arithmetic operators ...
    friend MyFloat sqrt(const MyFloat& x) { /* ... */ }
    friend MyFloat sin(const MyFloat& x) { /* ... */ }
};
```

Spatium uses `using std::sqrt; sqrt(x)` pattern, so ADL-discovered functions take precedence.

Built-in: `Real50`, `Real100` from Boost.Multiprecision (`core/precision.hpp`, built with `-DSPATIUM_BOOST=ON`).

---

## Adding a Quadric Surface

Use `Quadric<T>` (in `geometry/ray_surface.hpp`) with a custom 4×4 symmetric matrix:

```cpp
// Hyperboloid of one sheet: x²/a² + y²/b² - z²/c² = 1
template<Scalar T>
Quadric<T> hyperboloid(T a, T b, T c) {
    auto Q = Matrix<T, 4, 4>{};
    Q(0, 0) = T{1} / (a * a);
    Q(1, 1) = T{1} / (b * b);
    Q(2, 2) = T{-1} / (c * c);
    Q(3, 3) = T{-1};
    return {Q};
}
```

Then `ray_quadric(ray, hyperboloid(2, 3, 1))` works automatically — the polynomial solver handles the quadratic equation.

---

## Surface Adapter: Geometry Shape → Surface

`ShapeSurface` (in `geometry/surface_adapter.hpp`) converts geometry shapes into 2D manifolds for mesh operations:

```cpp
// Requires: DistanceQueryable + (HasPointNormal | HasConstNormal)
auto triangle_surface = ShapeSurface(my_triangle);
auto mesh = mesh::Mesh<decltype(triangle_surface)>{/* ... */};
```

The adapter auto-generates: `contains()`, `distance()`, `project()`, `normal()`, `exp_map()`, `log_map()`, `metric_at()`.
