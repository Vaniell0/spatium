// Raycast throughput benchmarks across BVH, brute force and analytical
// ray-quadric. Uses a single tessellated unit sphere to give BVH vs
// brute-force a fair comparison; quadric path uses the analytical sphere.

#include <benchmark/benchmark.h>
#include <spatium/spatial/bvh.hpp>
#include <spatium/geometry/triangle.hpp>
#include <spatium/geometry/line.hpp>
#include <spatium/geometry/ray_surface.hpp>
#include <spatium/geometry/ray_parametric.hpp>
#include <spatium/geometry/make.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/mesh/subdivision.hpp>
#include <spatium/spaces/sphere.hpp>
#include <spatium/spaces/parametric.hpp>
#include <spatium/io/build.hpp>
#include <random>
#include <vector>

using namespace spatium;
namespace bd = spatium::io::build;
using namespace spatium::geometry;
using namespace spatium::mesh;
using namespace spatium::spatial;

namespace {

std::vector<Triangle3> unit_sphere_triangles(std::size_t subdivisions) {
    Sphere<2> s;
    auto mesh = subdivide(icosahedron(s), s, subdivisions);
    std::vector<Triangle3> tris;
    tris.reserve(mesh.face_count());
    for (auto [a, b, c] : mesh.triangles())
        tris.push_back(Triangle3(a, b, c));
    return tris;
}

std::vector<Ray<3>> make_rays(std::size_t count, double spread = 0.6) {
    std::mt19937 rng(2026);
    std::uniform_real_distribution<double> j(-spread, spread);
    std::vector<Ray<3>> out;
    out.reserve(count);
    Vec3 origin{3, 3, 3};
    Vec3 aim{0, 0, 0};
    auto forward = (aim - origin).normalized();
    Vec3 up{0, 0, 1};
    auto right = forward.cross(up).normalized();
    up = right.cross(forward).normalized();
    while (out.size() < count) {
        auto d = Vec3{forward + right * j(rng) + up * j(rng)}.normalized();
        if (auto r = ray(origin, d)) out.push_back(*r);
    }
    return out;
}

} // namespace

// ── BVH ray_cast on tessellated sphere ────────────────────────

static void BM_RayCast_BVH_Sphere(benchmark::State& state) {
    auto tris = unit_sphere_triangles(static_cast<std::size_t>(state.range(0)));
    auto bvh = BVH<Triangle3>::build(tris);
    auto rays = make_rays(256);
    std::size_t i = 0;
    for (auto _ : state) {
        auto hit = bvh.ray_cast(rays[i++ & 255]);
        benchmark::DoNotOptimize(hit);
    }
    state.counters["tris"] = static_cast<double>(tris.size());
}
BENCHMARK(BM_RayCast_BVH_Sphere)->Arg(2)->Arg(4)->Arg(5);

static void BM_RayCast_Brute_Sphere(benchmark::State& state) {
    auto tris = unit_sphere_triangles(static_cast<std::size_t>(state.range(0)));
    auto rays = make_rays(256);
    std::size_t i = 0;
    for (auto _ : state) {
        auto& r = rays[i++ & 255];
        double best_t = 1e30;
        bool any = false;
        for (auto& tri : tris) {
            if (auto it = intersect(r, tri)) {
                double t = (*it - r.origin).dot(r.direction);
                if (t >= 0 && t < best_t) { best_t = t; any = true; }
            }
        }
        benchmark::DoNotOptimize(any);
    }
    state.counters["tris"] = static_cast<double>(tris.size());
}
BENCHMARK(BM_RayCast_Brute_Sphere)->Arg(2)->Arg(4)->Arg(5);

// ── Analytical ray-quadric ────────────────────────────────────

static void BM_RayQuadric_Sphere(benchmark::State& state) {
    auto q = Quadric<double>::sphere(1.0);
    auto rays = make_rays(256);
    std::size_t i = 0;
    for (auto _ : state) {
        auto hits = ray_quadric(rays[i++ & 255], q);
        benchmark::DoNotOptimize(hits);
    }
}
BENCHMARK(BM_RayQuadric_Sphere);

// The same call, but with every ray actually hitting. The default spread
// of 0.6 rad against a unit sphere seen from 5.2 away is mostly *misses*,
// which measures the cheap path and hides the expensive one: a hit is
// where the result container has to hold something. Measuring only the
// miss path is how "the allocation is the cost" stays unfalsified.
// 0.1 rad keeps the worst-case lateral offset inside the radius, so this
// is 100% hits by construction rather than by luck.
static void BM_RayQuadric_SphereAllHit(benchmark::State& state) {
    auto q = Quadric<double>::sphere(1.0);
    auto rays = make_rays(256, 0.1);
    std::size_t i = 0;
    for (auto _ : state) {
        auto hits = ray_quadric(rays[i++ & 255], q);
        benchmark::DoNotOptimize(hits);
    }
}
BENCHMARK(BM_RayQuadric_SphereAllHit);

static void BM_RayQuadric_Ellipsoid(benchmark::State& state) {
    auto q = Quadric<double>::ellipsoid(1.0, 0.5, 0.3);
    auto rays = make_rays(256);
    std::size_t i = 0;
    for (auto _ : state) {
        auto hits = ray_quadric(rays[i++ & 255], q);
        benchmark::DoNotOptimize(hits);
    }
}
BENCHMARK(BM_RayQuadric_Ellipsoid);

static void BM_RayQuadric_Cylinder(benchmark::State& state) {
    auto q = Quadric<double>::cylinder_z(0.5);
    auto rays = make_rays(256);
    std::size_t i = 0;
    for (auto _ : state) {
        auto hits = ray_quadric(rays[i++ & 255], q);
        benchmark::DoNotOptimize(hits);
    }
}
BENCHMARK(BM_RayQuadric_Cylinder);

static void BM_RayQuadric_Proximity(benchmark::State& state) {
    auto q = Quadric<double>::sphere(1.0);
    // Rays aimed away from sphere to hit the miss-proximity path
    auto rays = make_rays(256, 1.5);
    std::size_t i = 0;
    for (auto _ : state) {
        auto prox = ray_quadric_proximity(rays[i++ & 255], q);
        benchmark::DoNotOptimize(prox);
    }
}
BENCHMARK(BM_RayQuadric_Proximity);

// ── Analytical ray-torus ─────────────────────────────────────

static void BM_RayTorus(benchmark::State& state) {
    Torus<double> torus{.major_radius = 1.0, .minor_radius = 0.3};
    auto rays = make_rays(256);
    std::size_t i = 0;
    for (auto _ : state) {
        auto hits = ray_torus(rays[i++ & 255], torus);
        benchmark::DoNotOptimize(hits);
    }
}
BENCHMARK(BM_RayTorus);

static void BM_RayTorus_Proximity(benchmark::State& state) {
    Torus<double> torus{.major_radius = 1.0, .minor_radius = 0.3};
    auto rays = make_rays(256, 2.5);
    std::size_t i = 0;
    for (auto _ : state) {
        auto prox = ray_torus_proximity(rays[i++ & 255], torus);
        benchmark::DoNotOptimize(prox);
    }
}
BENCHMARK(BM_RayTorus_Proximity);

// ── Ray-ParametricSurface via Newton UV ──────────────────────
// Compares against BM_RayTorus (closed-form quartic) on the same torus.

static void BM_RayParametric_Torus(benchmark::State& state) {
    auto surf = make_torus<double>(1.0, 0.3);
    auto rays = make_rays(256);
    std::size_t i = 0;
    for (auto _ : state) {
        auto hits = ray_parametric(rays[i++ & 255], surf);
        benchmark::DoNotOptimize(hits);
    }
}
BENCHMARK(BM_RayParametric_Torus);

static void BM_RayParametric_Mobius(benchmark::State& state) {
    auto surf = make_mobius<double>(1.0, 0.3);
    auto rays = make_rays(256);
    std::size_t i = 0;
    for (auto _ : state) {
        auto hits = ray_parametric(rays[i++ & 255], surf);
        benchmark::DoNotOptimize(hits);
    }
}
BENCHMARK(BM_RayParametric_Mobius);

static void BM_RayParametric_TorusFirst(benchmark::State& state) {
    auto surf = make_torus<double>(1.0, 0.3);
    auto rays = make_rays(256);
    std::size_t i = 0;
    for (auto _ : state) {
        auto hit = ray_parametric_first(rays[i++ & 255], surf);
        benchmark::DoNotOptimize(hit);
    }
}
BENCHMARK(BM_RayParametric_TorusFirst);

// ── A scene of many small round objects, two ways ─────────────
//
// The shape examples/donut_demo.cpp actually builds: ~19 800 dust specks,
// each an icosahedron of 12 vertices and 20 faces. Rendered as triangles
// that is ~396 000 leaves; rendered as exact spheres it is 19 800.
//
// Both trees do the same culling work -- the question is only what sits
// at the leaf and how many of them there are. A single exact quadric hit
// costs about what a single triangle test costs (see BM_RayQuadric_Sphere
// vs BM_RayTriangleHit above), so any win here comes from leaf count and
// from never building the tessellation at all, not from the intersection
// being cheaper per call.
//
// Measured (12-core, loaded machine, so read the ratios):
//
//                     build      per ray    leaves     memory
//   triangles        291 ms       121 ns      396k     ~28.5 MB
//   exact spheres   11.8 ms       143 ns     19.8k      ~3.5 MB
//
// Per ray the analytic path is *slower* by about a fifth, and that is
// not a defect to fix: a sphere's AABB overlaps its neighbours more than
// a triangle's does, so traversal visits more nodes, and the leaf test
// itself is two matrix-vector products and two square roots against
// Moller-Trumbore's cross products and one division.
//
// The decision is made elsewhere. At 960x720 the whole frame costs
// 291 + 84 = 375 ms through triangles and 12 + 99 = 111 ms through
// spheres -- 3.4x, driven entirely by tree construction, which an
// animated scene pays again on every single frame while the per-ray cost
// is paid once per pixel. The tessellation also never exists, which is
// the eightfold memory difference and the ~20% of frame time
// benchmarks/bench_trace.cpp attributes to assembling those meshes.

namespace {

constexpr std::size_t kSpecks = 19800;
constexpr double kSpeckRadius = 0.014;

std::vector<Vec3> speck_centers() {
    std::mt19937 rng(4242);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::vector<Vec3> c;
    c.reserve(kSpecks);
    for (std::size_t i = 0; i < kSpecks; ++i) c.push_back(Vec3{u(rng), u(rng), u(rng)});
    return c;
}

std::vector<BoundedQuadric<double>> speck_quadrics() {
    std::vector<BoundedQuadric<double>> out;
    out.reserve(kSpecks);
    for (const auto& c : speck_centers())
        out.push_back(BoundedQuadric<double>::sphere(c, kSpeckRadius));
    return out;
}

std::vector<Triangle3> speck_triangles() {
    Sphere<2> unit;
    auto ico = icosahedron(unit);
    std::vector<Triangle3> out;
    out.reserve(kSpecks * ico.face_count());
    for (const auto& c : speck_centers())
        for (auto [a, b, d] : ico.triangles())
            out.push_back(Triangle3(Vec3{c + a * kSpeckRadius},
                                    Vec3{c + b * kSpeckRadius},
                                    Vec3{c + d * kSpeckRadius}));
    return out;
}

} // namespace

static void BM_SpeckScene_Build_Triangles(benchmark::State& state) {
    auto tris = speck_triangles();
    for (auto _ : state) benchmark::DoNotOptimize(BVH<Triangle3>::build(tris));
    state.counters["leaves"] = static_cast<double>(tris.size());
}
BENCHMARK(BM_SpeckScene_Build_Triangles)->Unit(benchmark::kMillisecond);

static void BM_SpeckScene_Build_Quadrics(benchmark::State& state) {
    auto qs = speck_quadrics();
    for (auto _ : state) benchmark::DoNotOptimize(BVH<BoundedQuadric<double>>::build(qs));
    state.counters["leaves"] = static_cast<double>(qs.size());
}
BENCHMARK(BM_SpeckScene_Build_Quadrics)->Unit(benchmark::kMillisecond);

static void BM_SpeckScene_RayCast_Triangles(benchmark::State& state) {
    auto tris = speck_triangles();
    auto bvh = BVH<Triangle3>::build(tris);
    auto rays = make_rays(256, 0.9);
    std::size_t i = 0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(bvh.ray_cast(rays[i++ % rays.size()]));
    }
    state.counters["leaves"] = static_cast<double>(tris.size());
}
BENCHMARK(BM_SpeckScene_RayCast_Triangles);

static void BM_SpeckScene_RayCast_Quadrics(benchmark::State& state) {
    auto qs = speck_quadrics();
    auto bvh = BVH<BoundedQuadric<double>>::build(qs);
    auto rays = make_rays(256, 0.9);
    std::size_t i = 0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(bvh.ray_cast(rays[i++ % rays.size()]));
    }
    state.counters["leaves"] = static_cast<double>(qs.size());
}
BENCHMARK(BM_SpeckScene_RayCast_Quadrics);

// ── The dough: one torus, or the 25 600 triangles it becomes ─────
//
// examples/donut_demo.cpp tessellates its dough at 160x80, which is
// 12 800 quads and so 25 600 triangles -- for one shape that has an
// exact closed form. Torus now satisfies Bounded (bounded by
// construction, unlike a general quadric, so no clip wrapper), which is
// what lets the same BVH hold it directly.
//
// This is a different question from the speck scene above. There the win
// came from leaf *count* at equal per-leaf cost, and the per-ray number
// went the wrong way. Here the leaf count collapses from 25 600 to 1
// while the leaf itself gets much more expensive -- a quartic solve
// rather than a Moller-Trumbore test -- so the two effects pull against
// each other and the ratio had to be measured rather than assumed.
//
// Measured (12-core, loaded machine, same ray set for both):
//
//                      build      per ray    leaves
//   25 600 triangles   10.8 ms     40.8 ns    25.6k
//   one exact torus       ~0       20.9 ns        1
//
// The per-ray cost went the right way here -- half, not a fifth worse.
// Corrected 2026-09-15 after BM_DSL_Tori below: that is an artefact of
// this case's extreme ratio, one leaf against 25 600, and not a general
// property. With 64 tori against 147 456 triangles the per-ray numbers
// come back to a wash (44 ns against 49). Read this 2x as "one leaf
// beats a deep tree", not as "exact hits are faster".
//
// The claim that survives all three measurements -- specks, this, and
// the DSL tori -- is about build and memory, not about rays: the
// tessellation never exists. Here that is ~1.8 MB of triangles replaced
// by a 40-byte struct and 10.8 ms a frame that is not spent.

namespace {

constexpr std::size_t kDoughU = 160, kDoughV = 80;

std::vector<Triangle3> dough_triangles() {
    auto surf = make_torus<double>(2.0, 1.0);
    auto m = tessellate(surf, kDoughU, kDoughV);
    std::vector<Triangle3> tris;
    tris.reserve(m.faces.size());
    for (const auto& f : m.faces)
        tris.push_back(Triangle3{m.vertices[f[0]], m.vertices[f[1]], m.vertices[f[2]]});
    return tris;
}

} // namespace

static void BM_Dough_Build_Triangles(benchmark::State& state) {
    auto tris = dough_triangles();
    for (auto _ : state) {
        auto bvh = BVH<Triangle3>::build(tris);
        benchmark::DoNotOptimize(bvh);
    }
    state.counters["leaves"] = static_cast<double>(tris.size());
}
BENCHMARK(BM_Dough_Build_Triangles)->Unit(benchmark::kMillisecond);

static void BM_Dough_RayCast_Triangles(benchmark::State& state) {
    auto bvh = BVH<Triangle3>::build(dough_triangles());
    auto rays = make_rays(256, 3.0);
    std::size_t i = 0;
    for (auto _ : state) benchmark::DoNotOptimize(bvh.ray_cast(rays[i++ % rays.size()]));
}
BENCHMARK(BM_Dough_RayCast_Triangles);

static void BM_Dough_RayCast_ExactTorus(benchmark::State& state) {
    std::vector<Torus<double>> one{Torus<double>{.major_radius = 2.0, .minor_radius = 1.0}};
    auto bvh = BVH<Torus<double>>::build(one);
    auto rays = make_rays(256, 3.0);
    std::size_t i = 0;
    for (auto _ : state) benchmark::DoNotOptimize(bvh.ray_cast(rays[i++ % rays.size()]));
    state.counters["leaves"] = 1.0;
}
BENCHMARK(BM_Dough_RayCast_ExactTorus);

// ── The same scene, through the DSL, at two levels ───────────────
//
// Every measurement above hand-builds its shapes. This one goes through
// io::build: a Trace of 64 tori, materialized, then bucketed by what
// Placed::render_level() says. That matters because the DSL is where the
// choice actually gets made -- a number measured on hand-built shapes
// says the exact path *could* be faster, and this one says the path a
// caller reaches through `.rendered_as()` is.
//
// Both buckets come from the same trace. The difference is one call.
//
// Measured (12-core, loaded machine):
//
//                    build       per ray    leaves
//   exact          0.011 ms      44.2 ns        64
//   tessellated      78.8 ms     48.5 ns   147 456
//
// Seven thousand times the build cost, for a tenth of the per-ray cost
// back. That ratio is the whole argument for the exact path, and it is
// the same argument the speck scene made -- the per-ray numbers are a
// wash in both, and what actually differs is that one of them never
// builds the triangles at all. An animated scene pays the build on every
// frame while the per-ray cost is paid once per pixel.

namespace {

constexpr std::size_t kTori = 64;

bd::Trace<double> tori_trace(bool force_tessellation) {
    bd::Trace<double> scene;
    std::vector<bd::Handle<double>> handles;
    for (std::size_t i = 0; i < kTori; ++i) {
        auto h = scene.torus(0.22, 0.07);
        if (force_tessellation) h.rendered_as(bd::RenderLevel::Tessellated);
        handles.push_back(h);
    }
    scene.compose(handles);
    return scene;
}

// The tori are built at the origin, so they are placed here instead --
// materialize() has no per-node placement yet (that is instancing, still
// open), and overlapping leaves would measure BVH degeneracy rather than
// the level choice.
//
// Placed as a cloud around the origin at the scale the speck scene uses,
// deliberately: make_rays() fires from (3,3,3) toward the origin, so a
// scene spread wider than that gets a ray set that mostly misses it, and
// the comparison becomes one of hit rates rather than of leaf cost. A
// first version of this benchmark laid the tori out on a 42-unit grid
// and reported the exact path as slower; that number was measuring the
// layout, not the level.
std::vector<Vec3> tori_centers() {
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::vector<Vec3> c;
    c.reserve(kTori);
    while (c.size() < kTori) {
        Vec3 p{u(rng), u(rng), u(rng)};
        if (p.norm() <= 1.0) c.push_back(p);
    }
    return c;
}

} // namespace

static void BM_DSL_Tori_Exact(benchmark::State& state) {
    auto scene = tori_trace(false);
    auto placed = bd::materialize(scene, scene.size() - 1);

    auto centers = tori_centers();
    std::vector<Torus<double>> tori;
    std::size_t i = 0;
    for (const auto& obj : placed) {
        if (obj.render_level() != bd::RenderLevel::Exact) continue;
        auto t = *obj.exact_as<Torus<double>>();
        t.center = centers[i++];
        tori.push_back(t);
    }
    auto bvh = BVH<Torus<double>>::build(tori);
    auto rays = make_rays(256, 0.9);
    std::size_t r = 0;
    for (auto _ : state) benchmark::DoNotOptimize(bvh.ray_cast(rays[r++ % rays.size()]));
    state.counters["leaves"] = static_cast<double>(tori.size());
}
BENCHMARK(BM_DSL_Tori_Exact);

static void BM_DSL_Tori_Tessellated(benchmark::State& state) {
    auto scene = tori_trace(true);
    auto placed = bd::materialize(scene, scene.size() - 1);

    auto centers = tori_centers();
    std::vector<Triangle3> tris;
    std::size_t i = 0;
    for (const auto& obj : placed) {
        auto shift = centers[i++];
        auto m = obj.mesh();
        for (const auto& f : m.faces)
            tris.push_back(Triangle3{Vec3{m.vertices[f[0]] + shift},
                                     Vec3{m.vertices[f[1]] + shift},
                                     Vec3{m.vertices[f[2]] + shift}});
    }
    auto bvh = BVH<Triangle3>::build(tris);
    auto rays = make_rays(256, 0.9);
    std::size_t r = 0;
    for (auto _ : state) benchmark::DoNotOptimize(bvh.ray_cast(rays[r++ % rays.size()]));
    state.counters["leaves"] = static_cast<double>(tris.size());
}
BENCHMARK(BM_DSL_Tori_Tessellated);

static void BM_DSL_Tori_Build_Exact(benchmark::State& state) {
    auto scene = tori_trace(false);
    auto centers = tori_centers();
    for (auto _ : state) {
        auto placed = bd::materialize(scene, scene.size() - 1);
        std::vector<Torus<double>> tori;
        std::size_t i = 0;
        for (const auto& obj : placed) {
            if (obj.render_level() != bd::RenderLevel::Exact) continue;
            auto t = *obj.exact_as<Torus<double>>();
            t.center = centers[i++];
            tori.push_back(t);
        }
        auto bvh = BVH<Torus<double>>::build(std::move(tori));
        benchmark::DoNotOptimize(bvh);
    }
}
BENCHMARK(BM_DSL_Tori_Build_Exact)->Unit(benchmark::kMillisecond);

static void BM_DSL_Tori_Build_Tessellated(benchmark::State& state) {
    auto scene = tori_trace(true);
    auto centers = tori_centers();
    for (auto _ : state) {
        auto placed = bd::materialize(scene, scene.size() - 1);
        std::vector<Triangle3> tris;
        std::size_t i = 0;
        for (const auto& obj : placed) {
            auto shift = centers[i++];
            auto m = obj.mesh();
            for (const auto& f : m.faces)
                tris.push_back(Triangle3{Vec3{m.vertices[f[0]] + shift},
                                         Vec3{m.vertices[f[1]] + shift},
                                         Vec3{m.vertices[f[2]] + shift}});
        }
        auto bvh = BVH<Triangle3>::build(std::move(tris));
        benchmark::DoNotOptimize(bvh);
    }
}
BENCHMARK(BM_DSL_Tori_Build_Tessellated)->Unit(benchmark::kMillisecond);
