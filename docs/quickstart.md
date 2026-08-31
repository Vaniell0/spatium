# Quick Start Guide

## Installation

### Nix (recommended)

```bash
git clone https://github.com/Vaniell0/spatium.git
cd spatium
nix develop       # enters dev shell with all deps
```

### Manual

Requirements: GCC 15+ or Clang 19+, CMake 3.28+, Catch2 v3, Boost (for multiprecision).

```bash
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
ninja -C build
ctest --test-dir build
```

## First Steps

### Include the umbrella header

```cpp
#include <spatium/spatium.hpp>
using namespace spatium;
using namespace spatium::literals;
```

This gives you everything: Vec, spaces, geometry, mesh, morphisms, Complex, polynomials, formatting.

### Vectors

```cpp
Vec3 a{1, 2, 3};
Vec3 b{4, 5, 6};

auto c = a + b;              // (5, 7, 9)
auto d = a * 2.0;            // (2, 4, 6)
auto dot = a.dot(b);         // 32
auto n = a.norm();            // 3.741...
auto u = a.normalized();      // unit vector
auto cross = a.cross(b);     // 3D only

// UDL syntax (algebra/literals.hpp)
auto v = 3.0_x + 2.0_y + 1.0_z;  // Vec3{3, 2, 1}

std::println("{}", a);        // (1, 2, 3)
```

### Geometric Primitives

```cpp
// Factory syntax (preferred) — see geometry/make.hpp
auto t = tri(Vec3{0,0,0}, Vec3{1,0,0}, Vec3{0,1,0});
auto s = seg(Vec3{0,0,0}, Vec3{1,1,1});
auto b = box(Vec3{-1,-1,-1}, Vec3{1,1,1});
auto r = *ray(Vec3{0,0,5}, Vec3{0,0,-1});    // unwrap Result
auto p = *plane(Vec3{0,0,1}, Vec3{0,0,0});

t.area();                     // 0.5
t.normal();                   // (0, 0, 1)
t.contains(Vec3{0.1, 0.1, 0}); // true
t.bounding_box();             // Box[(0,0,0) <-> (1,1,0)]
```

### Intersection

Two ways: free function or pipe operator.

```cpp
auto r = *ray(Vec3{0.25, 0.25, 5}, Vec3{0, 0, -1});
auto t = tri(Vec3{0,0,0}, Vec3{1,0,0}, Vec3{0,1,0});

auto hit1 = intersect(r, t);    // Result<Vec3>
auto hit2 = r | t;               // same — operator| in make.hpp:92

if (hit2)
    std::println("hit at {}", *hit2);
```

### Chaining with Result pipe-unwrap

```cpp
// line(...) | plane(...) | transform — auto-unwraps Result at each step
auto result = line(Vec3{0,0,0}, Vec3{1,0,0})
            | *plane(Vec3{1,0,0}, Vec3{3,0,0})
            | translate(Vec3{10,0,0});
// Result<Vec3> → intersect at (3,0,0) → translate to (13,0,0)
```

### Distance

```cpp
// Point-to-shape
auto d = distance(Vec3{5, 0, 0}, t);  // distance.hpp

// Shape-to-shape
auto d2 = distance(box1, box2);
auto d3 = distance(tri1, tri2);
auto d4 = distance(poly1, poly2);
```

### Polygon Boolean Operations

```cpp
auto a = poly<2, double>({{0,0}, {2,0}, {2,2}, {0,2}});
auto b = poly<2, double>({{1,1}, {3,1}, {3,3}, {1,3}});

auto inter = a & b;    // Result<Polygon> — intersection
auto diff  = a - b;    // Result<Polygon> — difference
auto uni   = a + b;    // Result<Polygon> — union (convex hull, 2D)
                       // (`|` is reserved for intersect/pipe, not union)
```

### Spaces and Points

```cpp
E3 space;
auto d = space.distance(Vec3{0,0,0}, Vec3{3,4,0});  // 5.0

// Typed points prevent mixing spaces
auto p = pt<E3>(Vec3{1, 2, 3});
auto q = pt<E3>(Vec3{4, 5, 6});
auto dist = p.distance_to(q, space);

// Sphere — great-circle geodesics
S2 sphere;
auto north = pt<S2>(Vec3{0, 0, 1});
auto east  = pt<S2>(Vec3{1, 0, 0});
auto geodesic_dist = north.distance_to(east, sphere);  // π/2
auto tangent = north.log(east, sphere);
auto midpoint = north.exp(tangent, 0.5, sphere);
```

### Morphisms

```cpp
auto scale = morph<E3, E3>([](const Vec3& p) { return p * 2.0; });
auto shift = morph<E3, E3>([](const Vec3& p) { return p + Vec3{10,0,0}; });
auto proj  = morph<E3, E2>([](const Vec3& p) -> Vec2 { return {p[0], p[1]}; });

// Pipeline: point | f | g
auto result = pt<E3>(Vec3{1, 2, 3}) | scale | shift | proj;  // P(12, 4)
```

### Transforms

```cpp
// Eager composition
auto t = translate(Vec3{1,0,0}) * rotate_z(45_deg) * scale(2.0);
auto p = t(Vec3{1, 0, 0});

// Lazy chains — no matrix multiply until needed
auto chain = lazy(translate(Vec3{1,0,0})) * lazy(rotate_z(45_deg)) * lazy(scale(2.0));
auto p2 = chain.apply(Vec3{1, 0, 0});     // sequential apply
auto mat = chain.collapse();               // combine to single matrix
```

### Complex Numbers and Polynomials

```cpp
Complex64 z{3, 4};           // 3 + 4i
z.magnitude();                // 5
z.conjugate();                // 3 - 4i
auto w = sqrt(Complex64{-1}); // i

// Solve x² + 1 = 0
auto roots = solve_quadratic(1.0, 0.0, 1.0);  // → {i, -i}

// Solve x³ - 6x² + 11x - 6 = 0
auto reals = real_roots_cubic(1.0, -6.0, 11.0, -6.0);  // → {1, 2, 3}
```

### Ray-Quadric Intersection

```cpp
auto sphere = Quadric<>::sphere(3.0);
auto r = *ray(Vec3{-10, 0, 0}, Vec3{1, 0, 0});

// Intersection: returns sorted hits
auto hits = ray_quadric(r, sphere);
for (auto& h : hits)
    std::println("t={:.1f} point={} normal={}", h.t, h.point, h.normal);

// Miss: proximity metric via complex roots
auto r2 = *ray(Vec3{-10, 5, 0}, Vec3{1, 0, 0});
auto prox = ray_quadric_proximity(r2, sphere);
if (prox)
    std::println("miss distance: {:.3f}", prox->miss_distance);
```

### Mesh Subdivision

```cpp
S2 sphere;
auto base = mesh::icosahedron(sphere);
auto refined = mesh::subdivide(base, sphere, 3);  // 642 vertices, 1280 faces

// All vertices stay on sphere
for (auto& v : refined.vertices)
    assert(sphere.contains(v));

// LOD chain
auto lod = mesh::LodChain<S2>::build(base, sphere, 4);
auto coarse = lod.at(0);   // 20 faces
auto fine   = lod.at(4);   // 5120 faces
```

### Parametric and Implicit Surfaces

```cpp
// Parametric: f(u,v) → R³
auto torus = make_torus(2.0, 0.5);
auto mesh = tessellate(torus, 32, 16);

// Implicit: F(x,y,z) = 0
auto gyroid = make_gyroid();
auto mesh2 = marching_cubes(gyroid, 64);
```

### Generic Algebra (Group/LieGroup/Ring)

```cpp
SO3 so3;
auto R = so3.rz(45_deg);

// Concept-constrained generic functions
auto R4 = algebra::power(so3, R, 4);              // R⁴ = 180° rotation
auto comm = algebra::commutator(so3, so3.rx(0.3), so3.ry(0.3));  // [Rx, Ry]
auto ad = algebra::adjoint(so3, R, Vec3{1, 0, 0}); // Ad_R(v)
```

### BVH Ray Casting

```cpp
std::vector<Triangle3> triangles = /* ... */;
auto bvh = BVH<Triangle3>::build(triangles);

auto r = *ray(Vec3{0, 0, 10}, Vec3{0, 0, -1});
auto hit = bvh.ray_cast(r);
if (hit)
    std::println("hit triangle {} at t={:.3f}", hit->index, hit->t);
```

### IO

```cpp
// OBJ
auto mesh = *load_obj("model.obj");
save_obj(mesh, "output.obj");

// STL
auto mesh2 = *load_stl("model.stl");
save_stl(mesh2, "output.stl");

// SVG (2D projection)
auto svg = mesh_to_svg(mesh, OrthoProjection{});
svg.save("view.svg");
```

### Arbitrary Precision

```cpp
#include <spatium/core/precision.hpp>

Euclidean<3, Real50> space;
Vec<Real50, 3> a{Real50{0}, Real50{0}, Real50{0}};
Vec<Real50, 3> b{Real50{3}, Real50{4}, Real50{0}};
auto d = space.distance(a, b);  // 5.000...000 (50 digits)
```

### Axiom Verification

```cpp
S2 sphere;
std::array samples = {Vec3{0,0,1}, Vec3{1,0,0}, Vec3{0,1,0}};

auto r = verify_metric(sphere, std::span{samples});
if (!r) std::println("FAILED: {}", r.failure);

auto r2 = verify_exp_log(sphere, std::span{samples});
if (!r2) std::println("FAILED: {}", r2.failure);
```
