#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/spatial/bvh.hpp>
#include <spatium/geometry/triangle.hpp>
#include <spatium/geometry/line.hpp>
#include <spatium/geometry/ray_surface.hpp>
#include <cmath>

using namespace spatium;
using namespace spatium::geometry;
using namespace spatium::spatial;
using Catch::Matchers::WithinAbs;

// ── Helpers ───────────────────────────────────────────────────

static std::vector<Triangle3> make_grid_triangles(int nx, int ny) {
    std::vector<Triangle3> tris;
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            double fx = x, fy = y;
            tris.push_back(Triangle3({fx, fy, 0}, {fx + 1, fy, 0}, {fx, fy + 1, 0}));
            tris.push_back(Triangle3({fx + 1, fy, 0}, {fx + 1, fy + 1, 0}, {fx, fy + 1, 0}));
        }
    }
    return tris;
}

// ── Build ─────────────────────────────────────────────────────

TEST_CASE("BVH build from triangles", "[bvh]") {
    auto tris = make_grid_triangles(5, 5);
    REQUIRE(tris.size() == 50);

    auto bvh = BVH<Triangle3>::build(tris);
    CHECK(bvh.node_count() > 0);
    CHECK(bvh.shapes().size() == 50);
}

TEST_CASE("BVH build empty", "[bvh]") {
    auto bvh = BVH<Triangle3>::build({});
    CHECK(bvh.node_count() == 0);
    CHECK(bvh.shapes().empty());
}

TEST_CASE("BVH build single element", "[bvh]") {
    std::vector<Triangle3> tris = {Triangle3({0, 0, 0}, {1, 0, 0}, {0, 1, 0})};
    auto bvh = BVH<Triangle3>::build(tris);
    CHECK(bvh.node_count() == 1);
    CHECK(bvh.shapes().size() == 1);
}

// ── ray_cast ──────────────────────────────────────────────────

TEST_CASE("BVH ray_cast hit", "[bvh]") {
    auto tris = make_grid_triangles(5, 5);
    auto bvh = BVH<Triangle3>::build(tris);

    // Shoot ray down at (2.5, 2.5)
    auto ray = *Ray3::from(Vec3{2.5, 2.5, 10.0}, Vec3{0.0, 0.0, -1.0});
    auto hit = bvh.ray_cast(ray);
    REQUIRE(hit.has_value());
    CHECK_THAT(hit->point[2], WithinAbs(0.0, 1e-10));
    CHECK_THAT(hit->t, WithinAbs(10.0, 1e-10));
}

TEST_CASE("BVH ray_cast miss", "[bvh]") {
    auto tris = make_grid_triangles(5, 5);
    auto bvh = BVH<Triangle3>::build(tris);

    // Shoot ray above and parallel
    auto ray = *Ray3::from(Vec3{2.5, 2.5, 10.0}, Vec3{1.0, 0.0, 0.0});
    auto hit = bvh.ray_cast(ray);
    CHECK_FALSE(hit.has_value());
}

TEST_CASE("BVH ray_cast empty", "[bvh]") {
    auto bvh = BVH<Triangle3>::build({});
    auto ray = *Ray3::from(Vec3{0, 0, 1}, Vec3{0, 0, -1});
    CHECK_FALSE(bvh.ray_cast(ray).has_value());
}

TEST_CASE("BVH ray_cast finds closest hit", "[bvh]") {
    // Two triangles at different Z: z=0 and z=-5
    std::vector<Triangle3> tris = {
        Triangle3({-1, -1, 0}, {1, -1, 0}, {0, 1, 0}),
        Triangle3({-1, -1, -5}, {1, -1, -5}, {0, 1, -5}),
    };
    auto bvh = BVH<Triangle3>::build(tris);
    auto ray = *Ray3::from(Vec3{0, 0, 10}, Vec3{0, 0, -1});
    auto hit = bvh.ray_cast(ray);
    REQUIRE(hit.has_value());
    CHECK_THAT(hit->point[2], WithinAbs(0.0, 1e-10));  // closest
}

// ── ray_test ──────────────────────────────────────────────────

TEST_CASE("BVH ray_test hit", "[bvh]") {
    auto tris = make_grid_triangles(5, 5);
    auto bvh = BVH<Triangle3>::build(tris);
    auto ray = *Ray3::from(Vec3{2.5, 2.5, 10.0}, Vec3{0.0, 0.0, -1.0});
    CHECK(bvh.ray_test(ray));
}

TEST_CASE("BVH ray_test miss", "[bvh]") {
    auto tris = make_grid_triangles(5, 5);
    auto bvh = BVH<Triangle3>::build(tris);
    auto ray = *Ray3::from(Vec3{100.0, 100.0, 10.0}, Vec3{0.0, 0.0, -1.0});
    CHECK_FALSE(bvh.ray_test(ray));
}

// ── nearest ───────────────────────────────────────────────────

TEST_CASE("BVH nearest point", "[bvh]") {
    auto tris = make_grid_triangles(5, 5);
    auto bvh = BVH<Triangle3>::build(tris);

    // Point directly above (2.5, 2.5, 3)
    auto result = bvh.nearest(Vec3{2.5, 2.5, 3.0});
    REQUIRE(result.has_value());
    CHECK_THAT(result->distance, WithinAbs(3.0, 1e-10));
    CHECK_THAT(result->point[2], WithinAbs(0.0, 1e-10));
}

TEST_CASE("BVH nearest point on surface", "[bvh]") {
    auto tris = make_grid_triangles(5, 5);
    auto bvh = BVH<Triangle3>::build(tris);

    auto result = bvh.nearest(Vec3{2.5, 2.5, 0.0});
    REQUIRE(result.has_value());
    CHECK_THAT(result->distance, WithinAbs(0.0, 1e-10));
}

TEST_CASE("BVH nearest empty", "[bvh]") {
    auto bvh = BVH<Triangle3>::build({});
    CHECK_FALSE(bvh.nearest(Vec3{0, 0, 0}).has_value());
}

// ── query_box ─────────────────────────────────────────────────

TEST_CASE("BVH query_box all", "[bvh]") {
    auto tris = make_grid_triangles(5, 5);
    auto bvh = BVH<Triangle3>::build(tris);

    Box3 query{Vec3{-1, -1, -1}, Vec3{6, 6, 1}};
    auto result = bvh.query_box(query);
    CHECK(result.size() == 50);
}

TEST_CASE("BVH query_box partial", "[bvh]") {
    auto tris = make_grid_triangles(5, 5);
    auto bvh = BVH<Triangle3>::build(tris);

    // Query a small region — should return subset
    Box3 query{Vec3{0.5, 0.5, -0.5}, Vec3{1.5, 1.5, 0.5}};
    auto result = bvh.query_box(query);
    CHECK(result.size() > 0);
    CHECK(result.size() < 50);
}

TEST_CASE("BVH query_box miss", "[bvh]") {
    auto tris = make_grid_triangles(5, 5);
    auto bvh = BVH<Triangle3>::build(tris);

    Box3 query{Vec3{100, 100, 100}, Vec3{200, 200, 200}};
    auto result = bvh.query_box(query);
    CHECK(result.empty());
}

TEST_CASE("BVH query_box empty tree", "[bvh]") {
    auto bvh = BVH<Triangle3>::build({});
    Box3 query{Vec3{0, 0, 0}, Vec3{1, 1, 1}};
    CHECK(bvh.query_box(query).empty());
}

// ── Hit barycentric + normal ─────────────────────────────────

TEST_CASE("BVH ray_cast Triangle3 Hit barycentric + normal", "[bvh]") {
    std::vector<Triangle3> tris{Triangle3({0,0,0}, {2,0,0}, {0,2,0})};
    auto bvh = BVH<Triangle3>::build(tris);

    // Centroid ray: u ≈ 1/3, v ≈ 1/3, normal = +Z
    Ray3 ray{Vec3{2.0/3, 2.0/3, 1}, Vec3{0, 0, -1}};
    auto hit = bvh.ray_cast(ray);
    REQUIRE(hit);
    CHECK_THAT(hit->u, WithinAbs(1.0/3, 1e-9));
    CHECK_THAT(hit->v, WithinAbs(1.0/3, 1e-9));
    CHECK_THAT(hit->normal[0], WithinAbs(0.0, 1e-9));
    CHECK_THAT(hit->normal[1], WithinAbs(0.0, 1e-9));
    CHECK_THAT(hit->normal[2], WithinAbs(1.0, 1e-9));

    // Corner at vertex 1 → u=1, v=0
    Ray3 r1{Vec3{2, 0, 1}, Vec3{0, 0, -1}};
    auto h1 = bvh.ray_cast(r1);
    REQUIRE(h1);
    CHECK_THAT(h1->u, WithinAbs(1.0, 1e-9));
    CHECK_THAT(h1->v, WithinAbs(0.0, 1e-9));

    // Corner at vertex 2 → u=0, v=1
    Ray3 r2{Vec3{0, 2, 1}, Vec3{0, 0, -1}};
    auto h2 = bvh.ray_cast(r2);
    REQUIRE(h2);
    CHECK_THAT(h2->u, WithinAbs(0.0, 1e-9));
    CHECK_THAT(h2->v, WithinAbs(1.0, 1e-9));
}

// ── Rays that graze a bound exactly ───────────────────────────
//
// Seventeen tests above and not one of them touches a bound. That is the
// gap this closes, and it is worth closing *before* anything narrows the
// bounds rather than after.
//
// A BVH box is an accelerator, so its only obligation is to be
// conservative: it may be looser than the shape and must never be
// tighter. A box that has shrunk by one ulp rejects a ray before the leaf
// test ever runs, and the shape then has a hole in it that looks like a
// shading bug rather than like a rejected ray -- the same shape of
// failure as `miss = 0.0`, a value that reads as an answer and means
// something else.
//
// The cases below all sit exactly on a bound, where "exactly" is
// representable: integers and halves, so the arithmetic is not itself in
// question.

TEST_CASE("BVH: a ray in the plane of a bound's face still reaches the shape",
          "[bvh]") {
    // A unit sphere at the origin bounds to [-1, 1]^3. A ray along +x at
    // y = 1 lies in the plane of the box's own face and is tangent to the
    // sphere: a genuine, single, touching hit. Narrow the bound inward by
    // any amount and the slab test discards the ray before the sphere is
    // ever asked.
    std::vector<BoundedQuadric<double>> shapes{BoundedQuadric<double>::sphere(1.0)};
    auto bvh = BVH<BoundedQuadric<double>>::build(shapes);

    Ray<3, double> grazing{Vec<double, 3>{-4.0, 1.0, 0.0}, Vec<double, 3>{1.0, 0.0, 0.0}};
    auto hit = bvh.ray_cast(grazing);
    REQUIRE(hit.has_value());
    CHECK_THAT(hit->point[1], WithinAbs(1.0, 1e-9));

    // And one ulp outside really does miss, so the case above is not
    // passing because everything passes.
    Ray<3, double> just_past{Vec<double, 3>{-4.0, std::nextafter(1.0, 2.0), 0.0},
                             Vec<double, 3>{1.0, 0.0, 0.0}};
    CHECK_FALSE(bvh.ray_cast(just_past).has_value());
}

TEST_CASE("BVH: a ray aimed at a bound's corner still reaches the shape", "[bvh]") {
    // The corner is the worst case for a narrowed bound, because all three
    // slabs are on their limit at once. This triangle's AABB corner is its
    // own vertex, so hitting the corner is hitting the shape.
    std::vector<Triangle3> tris{Triangle3({0, 0, 0}, {2, 0, 0}, {0, 2, 0})};
    auto bvh = BVH<Triangle3>::build(tris);

    Ray<3, double> at_corner{Vec<double, 3>{0.0, 0.0, 5.0}, Vec<double, 3>{0.0, 0.0, -1.0}};
    auto hit = bvh.ray_cast(at_corner);
    REQUIRE(hit.has_value());
    CHECK_THAT(hit->point[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(hit->point[1], WithinAbs(0.0, 1e-12));
}

TEST_CASE("BVH: a flat shape has a zero-thickness bound and is still found", "[bvh]") {
    // Every triangle in a plane gives an AABB with no extent in one axis,
    // so the slab test on that axis is entirely a question of whether the
    // comparison is inclusive. This is the degenerate bound that any
    // change to how bounds are stored will meet first.
    auto tris = make_grid_triangles(4, 4);
    auto bvh = BVH<Triangle3>::build(tris);

    Ray<3, double> down{Vec<double, 3>{1.25, 1.25, 3.0}, Vec<double, 3>{0.0, 0.0, -1.0}};
    auto hit = bvh.ray_cast(down);
    REQUIRE(hit.has_value());
    CHECK_THAT(hit->point[2], WithinAbs(0.0, 1e-12));

    // A ray *in* the plane is refused, and by the leaf rather than by the
    // box. Asserted the other way round first, which was wrong: a
    // coplanar ray meets a triangle in nothing or in a segment, never in
    // a point, so Möller-Trumbore's vanishing determinant is the correct
    // answer and not a numerical accident. Pinned because it was
    // untested, and because the next person to narrow a bound will want
    // to know that this `false` is the leaf's and not the tree's.
    Ray<3, double> along{Vec<double, 3>{-1.0, 1.0, 0.0}, Vec<double, 3>{1.0, 0.0, 0.0}};
    CHECK_FALSE(bvh.ray_test(along));

    // What the zero-width slab must not do is reject a ray that arrives
    // at a shallow angle and genuinely crosses the plane. "Shallow" has
    // to mean shallow-but-real: a first attempt used a slope of 1e-9,
    // which puts Möller-Trumbore's determinant under its own epsilon, so
    // the leaf refused it as parallel and the test was measuring the leaf
    // again rather than the box. Worth recording, because that is two
    // ways in a row for a "grazing" test to accidentally examine the
    // wrong half.
    //
    // This one crosses z = 0 at (1, 1.25), well inside a triangle rather
    // than on a vertex, with a determinant around 5e-3.
    Ray<3, double> shallow{Vec<double, 3>{-1.0, 1.25, 0.01},
                           Vec<double, 3>{1.0, 0.0, -0.005}};
    CHECK(bvh.ray_test(shallow));
}

// ── A tree deeper than the traversal stack ──────────────────────
//
// Triangles across the x axis at x = 1.5^i: the centroids spread over a
// hundred orders of magnitude, so every binned split peels a few of the
// largest off and the tree is a long spine. A ray down the axis crosses
// every box, and near-first traversal leaves the far child of each level
// on the stack. The stack used to be a fixed array of 64 with no bound
// check, which such a tree overruns; run under AddressSanitizer this test
// reported the write past its end.
namespace {

std::vector<Triangle3> spine(int n) {
    std::vector<Triangle3> tris;
    for (int i = 0; i < n; ++i) {
        const double x = std::pow(1.5, i);
        tris.push_back(Triangle3({x, -1, -1}, {x, 2, -1}, {x, -1, 2}));
    }
    return tris;
}

template<typename Tree>
std::size_t depth_of(const Tree& bvh, std::size_t node = 0) {
    const auto& n = bvh.nodes()[node];
    if (n.count > 0) return 1;
    return 1 + std::max(depth_of(bvh, node + 1), depth_of(bvh, n.first));
}

}  // namespace

TEST_CASE("BVH queries agree with brute force on a tree deeper than 64", "[bvh]") {
    const auto tris = spine(600);
    const auto bvh = BVH<Triangle3>::build(tris);
    INFO("depth " << depth_of(bvh));
    REQUIRE(depth_of(bvh) > 64);
    CHECK(bvh.depth() == depth_of(bvh));   // what the walk sizes its stack by

    const Ray<3, double> down_axis{{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
    const auto hit = bvh.ray_cast(down_axis);
    REQUIRE(hit);
    CHECK(hit->index == 0);
    CHECK_THAT(hit->t, WithinAbs(2.0, 1e-12));
    CHECK(bvh.ray_test(down_axis));

    // From beyond the far end looking back: the nearest is the last one.
    const double far = std::pow(1.5, 600);
    const Ray<3, double> back{{far * 2.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}};
    const auto back_hit = bvh.ray_cast(back);
    REQUIRE(back_hit);
    CHECK(back_hit->index == 599);

    const auto near = bvh.nearest(Vec<double, 3>{0.0, 0.0, 0.0});
    REQUIRE(near);
    CHECK(near->index == 0);

    const auto all = bvh.query_box(Box<3, double>{{0.0, -2.0, -2.0}, {far * 2.0, 3.0, 3.0}});
    CHECK(all.size() == 600);
}

// ── A ray that ends ─────────────────────────────────────────────
//
// A shadow ray stops at the light and a swept query at the end of its
// step; an occluder beyond either is not an occluder. Without a t_max
// the tree could only answer for a ray that never ends.
TEST_CASE("BVH ray queries stop at t_max", "[bvh]") {
    std::vector<Triangle3> tris{Triangle3({1, -1, -1}, {1, 2, -1}, {1, -1, 2}),
                                Triangle3({3, -1, -1}, {3, 2, -1}, {3, -1, 2})};
    const auto bvh = BVH<Triangle3>::build(tris);
    const Ray<3, double> ray{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};

    CHECK(bvh.ray_cast(ray, 2.0)->index == 0);
    CHECK_FALSE(bvh.ray_cast(ray, 0.5));
    CHECK(bvh.ray_test(ray, 2.0));
    CHECK_FALSE(bvh.ray_test(ray, 0.5));

    // Between the two, looking at the far one through a light at 2.5.
    const Ray<3, double> shadow{{2.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
    CHECK_FALSE(bvh.ray_test(shadow, 0.5));
    CHECK(bvh.ray_test(shadow, 1.5));
}

// ── Any-hit over the shapes ray_cast already handles ────────────
//
// ray_cast reaches a shape through ray_hit, ray_test through intersect(),
// so a tree of exact quadrics could answer "what does this ray hit first"
// and not "does it hit anything": BVH<BoundedQuadric>::ray_test did not
// compile. It goes through ray_hit now, the same as ray_cast.
TEST_CASE("BVH any-hit works for every shape first-hit works for", "[bvh]") {
    std::vector<BoundedQuadric<double>> balls{BoundedQuadric<double>::sphere(1.0)};
    const auto bvh = BVH<BoundedQuadric<double>>::build(balls);
    CHECK(bvh.ray_test(Ray<3, double>{{-5.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}));
    CHECK_FALSE(bvh.ray_test(Ray<3, double>{{-5.0, 3.0, 0.0}, {1.0, 0.0, 0.0}}));
    CHECK_FALSE(bvh.ray_test(Ray<3, double>{{-5.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}, 3.0));
    CHECK(bvh.ray_cast(Ray<3, double>{{-5.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}));
}

// ── The same tree over balls ────────────────────────────────────
//
// One hierarchy, the bound swapped: a ball tree must give the answers the
// box tree gives, which are the brute-force answers. First hit by ray and
// nearest point, over random rays and points, on a scene with both
// scattered triangles and a flat grid.
TEST_CASE("A BVH over balls answers as the BVH over boxes does", "[bvh][bound]") {
    std::vector<Triangle3> tris = make_grid_triangles(12, 12);
    std::uint64_t state = 12345;
    auto uniform = [&state] {
        state = state * 6364136223846793005ull + 1442695040888963407ull;
        return static_cast<double>(state >> 11) / static_cast<double>(1ull << 53);
    };
    for (int i = 0; i < 400; ++i) {
        const double x = uniform() * 12, y = uniform() * 12, z = uniform() * 4;
        tris.push_back(Triangle3({x, y, z}, {x + 0.4, y, z + 0.1}, {x, y + 0.4, z + 0.3}));
    }
    const auto boxes = BVH<Triangle3>::build(tris);
    const auto balls = BVH<Triangle3, Ball<3, double>>::build(tris);

    for (int i = 0; i < 2000; ++i) {
        const Vec<double, 3> o{uniform() * 12, uniform() * 12, 6.0};
        const Vec<double, 3> d{uniform() - 0.5, uniform() - 0.5, -1.0};
        const auto ray = *Ray<3, double>::from(o, d);
        const auto a = boxes.ray_cast(ray);
        const auto b = balls.ray_cast(ray);
        REQUIRE(a.has_value() == b.has_value());
        if (a) {
            CHECK(a->index == b->index);
            CHECK(a->t == b->t);
        }
        CHECK(boxes.ray_test(ray, 3.0) == balls.ray_test(ray, 3.0));

        const Vec<double, 3> p{uniform() * 14 - 1, uniform() * 14 - 1, uniform() * 6 - 1};
        const auto na = boxes.nearest(p);
        const auto nb = balls.nearest(p);
        REQUIRE(na);
        REQUIRE(nb);
        CHECK_THAT(nb->distance, WithinAbs(na->distance, 1e-12));
    }
}

TEST_CASE("A merged ball contains both of its children", "[bvh][bound]") {
    using B = Ball<3, double>;
    const B a{{0.0, 0.0, 0.0}, 1.0}, b{{3.0, 0.0, 0.0}, 0.5}, inner{{0.2, 0.0, 0.0}, 0.1};
    const auto m = merge(a, b);
    CHECK(Vec<double, 3>{a.c - m.c}.norm() + a.r <= m.r);
    CHECK(Vec<double, 3>{b.c - m.c}.norm() + b.r <= m.r);
    CHECK(merge(a, inner).r == a.r);   // one inside the other: the outer one
    CHECK_THAT(m.r, WithinAbs(2.25, 1e-12));
}
