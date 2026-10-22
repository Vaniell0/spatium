#include <catch2/catch_test_macros.hpp>
#include <spatium/io/build.hpp>
#include <spatium/render/cooked_scene.hpp>

#include <bit>
#include <cstdint>
#include <vector>

using namespace spatium;
namespace bd = spatium::io::build;

namespace {

std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }

// A torus with no closed form (so it is flattened to triangles), with
// small cylinders scattered over it.
bd::Cooked<double> torus_scene(double minor_r, std::size_t u, std::size_t v, std::size_t count,
                               std::uint32_t seed) {
    bd::Trace<double> tr;
    auto torus = tr.torus(2.0, minor_r, u, v);
    auto item = tr.cylinder(0.011, 0.062, 8, 2);
    auto sc = tr.scatter(item, torus, count, seed, 0.45, bd::SeatAxis::X);
    auto root = tr.compose({torus, sc});
    return bd::cook(tr, root.index, 0.0);
}

bool same_tree(const render::CookedScene<double>& a, const render::CookedScene<double>& b) {
    const auto &na = a.triangles.nodes(), &nb = b.triangles.nodes();
    if (na.size() != nb.size() || a.triangles.prim_indices() != b.triangles.prim_indices() ||
        a.prims.size() != b.prims.size() || a.triangles.depth() != b.triangles.depth())
        return false;
    for (std::size_t i = 0; i < na.size(); ++i) {
        if (na[i].first != nb[i].first || na[i].count != nb[i].count) return false;
        for (int k = 0; k < 3; ++k)
            if (bits(na[i].bounds.min_corner[k]) != bits(nb[i].bounds.min_corner[k]) ||
                bits(na[i].bounds.max_corner[k]) != bits(nb[i].bounds.max_corner[k]))
                return false;
    }
    return true;
}

// A handful of rays across the torus, from the side and from above.
void same_hits(const render::CookedScene<double>& a, const render::CookedScene<double>& b) {
    for (int i = -4; i <= 4; ++i) {
        const double x = 0.37 * i;
        geometry::Ray<3, double> side{Vec<double, 3>{x, -6.0, 0.1 * i}, Vec<double, 3>{0, 1, 0}};
        geometry::Ray<3, double> top{Vec<double, 3>{2.0 + 0.2 * x, 0.3 * x, 5.0},
                                     Vec<double, 3>{0, 0, -1}};
        for (const auto& r : {side, top}) {
            const auto ha = a.triangles.ray_cast(r), hb = b.triangles.ray_cast(r);
            REQUIRE(ha.has_value() == hb.has_value());
            if (ha) {
                CHECK(ha->index == hb->index);
                CHECK(bits(ha->t) == bits(hb->t));
            }
        }
    }
}

}  // namespace

TEST_CASE("a cache reuses the triangle tree when only the scatter changes", "[render][cache]") {
    render::StaticCache<double> cache;
    const auto first = torus_scene(1.0, 48, 24, 200, 42);
    const auto second = torus_scene(1.0, 48, 24, 350, 43);

    auto a = render::lay_out(first, {}, nullptr, &cache);
    CHECK(cache.hits() == 0);
    CHECK(cache.misses() == 1);
    auto b = render::lay_out(second, {}, nullptr, &cache);
    CHECK(cache.hits() == 1);
    CHECK(cache.misses() == 1);

    // Same scene as one built with no cache, tree and instances alike.
    const auto plain = render::lay_out(second);
    CHECK(same_tree(b, plain));
    CHECK(same_tree(a, render::lay_out(first)));
    same_hits(b, plain);
    CHECK(b.looks.size() == plain.looks.size());
    CHECK(b.looks.size() == 350);
    CHECK(b.instances.shapes().size() == plain.instances.shapes().size());
    for (std::size_t i = 0; i < plain.prims.size(); i += 97)
        CHECK(bits(b.prims[i].tri[0][0]) == bits(plain.prims[i].tri[0][0]));
}

TEST_CASE("changing the static geometry misses the cache", "[render][cache]") {
    render::StaticCache<double> cache;
    (void)render::lay_out(torus_scene(1.0, 48, 24, 50, 1), {}, nullptr, &cache);
    CHECK(cache.misses() == 1);

    // Radius, then each step count: each changes the triangles.
    (void)render::lay_out(torus_scene(0.9, 48, 24, 50, 1), {}, nullptr, &cache);
    (void)render::lay_out(torus_scene(1.0, 40, 24, 50, 1), {}, nullptr, &cache);
    (void)render::lay_out(torus_scene(1.0, 48, 20, 50, 1), {}, nullptr, &cache);
    CHECK(cache.hits() == 0);
    CHECK(cache.misses() == 4);

    // And the changed torus is its own correct tree, then hits.
    const auto changed = torus_scene(0.9, 48, 24, 50, 2);
    auto c = render::lay_out(changed, {}, nullptr, &cache);
    CHECK(cache.hits() == 1);
    const auto plain = render::lay_out(changed);
    CHECK(same_tree(c, plain));
    same_hits(c, plain);
}

TEST_CASE("a cache keeps at most its capacity", "[render][cache]") {
    render::StaticCache<double> cache(2);
    for (double r : {0.8, 0.9, 1.0})
        (void)render::lay_out(torus_scene(r, 24, 12, 10, 1), {}, nullptr, &cache);
    CHECK(cache.size() == 2);
    (void)render::lay_out(torus_scene(0.8, 24, 12, 10, 1), {}, nullptr, &cache);  // evicted
    CHECK(cache.hits() == 0);
    (void)render::lay_out(torus_scene(1.0, 24, 12, 10, 1), {}, nullptr, &cache);
    CHECK(cache.hits() == 1);
}
