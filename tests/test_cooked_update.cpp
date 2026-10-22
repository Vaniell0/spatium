#include <catch2/catch_test_macros.hpp>
#include <spatium/io/build.hpp>
#include <spatium/render/cooked_scene.hpp>

#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

using namespace spatium;
namespace bd = spatium::io::build;

namespace {

bd::Cooked<double> scatter_scene(std::size_t count, std::uint32_t seed) {
    bd::Trace<double> tr;
    auto torus = tr.torus(2.0, 1.0, 48, 24);
    auto item = tr.cylinder(0.011, 0.062, 8, 2);
    auto sc = tr.scatter(item, torus, count, seed, 0.45, bd::SeatAxis::X);
    auto root = tr.compose({torus, sc});
    return bd::cook(tr, root.index, 0.0);
}

// Rays aimed at the torus from all around, so a fair share reach an
// instance.
std::vector<geometry::Ray<3, double>> rays(std::size_t n) {
    std::mt19937_64 g(7);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::vector<geometry::Ray<3, double>> out;
    for (std::size_t i = 0; i < n; ++i) {
        const double a = u(g) * 6.2831853, r = 5.0 + u(g);
        const Vec<double, 3> o{r * std::cos(a), r * std::sin(a), (u(g) - 0.5) * 4.0};
        const double b = u(g) * 6.2831853;
        const Vec<double, 3> target{2.0 * std::cos(b), 2.0 * std::sin(b), (u(g) - 0.5) * 1.5};
        out.push_back(*geometry::Ray<3, double>::from(o, Vec<double, 3>{target - o}));
    }
    return out;
}

// The two scenes must agree on every ray: same distance, same instance
// (placement and look), same triangle hit.
std::size_t same_hits(const render::CookedScene<double>& a, const render::CookedScene<double>& b) {
    REQUIRE(a.instances.size() == b.instances.size());
    REQUIRE(a.looks.size() == a.instances.size());
    REQUIRE(b.looks.size() == b.instances.size());
    std::size_t instance_hits = 0;
    for (const auto& r : rays(1500)) {
        const auto ha = a.instances.ray_cast(r), hb = b.instances.ray_cast(r);
        REQUIRE(ha.has_value() == hb.has_value());
        if (ha) {
            ++instance_hits;
            REQUIRE(ha->t == hb->t);
            const auto &ia = a.instances.shapes()[ha->index], &ib = b.instances.shapes()[hb->index];
            REQUIRE(ia.translation[0] == ib.translation[0]);
            REQUIRE(ia.translation[1] == ib.translation[1]);
            REQUIRE(ia.translation[2] == ib.translation[2]);
            REQUIRE(a.looks[ha->index].color[0] == b.looks[hb->index].color[0]);
            REQUIRE(a.looks[ha->index].opacity == b.looks[hb->index].opacity);
        }
        const auto ta = a.triangles.ray_cast(r), tb = b.triangles.ray_cast(r);
        REQUIRE(ta.has_value() == tb.has_value());
        if (ta) REQUIRE(ta->t == tb->t);
        REQUIRE(a.instances.ray_test(r) == b.instances.ray_test(r));
    }
    return instance_hits;
}

}  // namespace

TEST_CASE("updating the instance tree gives the hits of a fresh lay_out", "[render][update]") {
    render::StaticCache<double> cache;
    const std::size_t N = 3000;
    render::InstanceUpdate u;

    SECTION("one more instance") {
        auto prev = render::lay_out(scatter_scene(N, 42), {}, nullptr, &cache);
        auto next = render::lay_out(std::move(prev), scatter_scene(N + 1, 42), {}, nullptr, &cache, &u);
        CHECK(u.reused);
        CHECK_FALSE(u.rebuilt);
        CHECK(u.inserted == 1);
        CHECK(u.removed == 0);
        CHECK(u.kept == N);
        CHECK(same_hits(next, render::lay_out(scatter_scene(N + 1, 42), {}, nullptr, &cache)) > 100);
    }

    SECTION("one fewer instance") {
        auto prev = render::lay_out(scatter_scene(N, 42), {}, nullptr, &cache);
        auto next = render::lay_out(std::move(prev), scatter_scene(N - 1, 42), {}, nullptr, &cache, &u);
        CHECK_FALSE(u.rebuilt);
        CHECK(u.removed == 1);
        CHECK(u.inserted == 0);
        CHECK(same_hits(next, render::lay_out(scatter_scene(N - 1, 42), {}, nullptr, &cache)) > 100);
    }

    SECTION("nothing changed") {
        auto prev = render::lay_out(scatter_scene(N, 42), {}, nullptr, &cache);
        auto next = render::lay_out(std::move(prev), scatter_scene(N, 42), {}, nullptr, &cache, &u);
        CHECK(u.kept == N);
        CHECK(u.removed + u.inserted == 0);
        CHECK(same_hits(next, render::lay_out(scatter_scene(N, 42), {}, nullptr, &cache)) > 100);
    }

    SECTION("a new seed changes everything: built again") {
        auto prev = render::lay_out(scatter_scene(N, 42), {}, nullptr, &cache);
        auto next = render::lay_out(std::move(prev), scatter_scene(N, 43), {}, nullptr, &cache, &u);
        CHECK(u.rebuilt);
        CHECK(same_hits(next, render::lay_out(scatter_scene(N, 43), {}, nullptr, &cache)) > 100);
    }

    SECTION("half the instances removed") {
        auto prev = render::lay_out(scatter_scene(N, 42), {}, nullptr, &cache);
        auto next = render::lay_out(std::move(prev), scatter_scene(N / 2, 42), {}, nullptr, &cache, &u);
        CHECK(u.removed == N - N / 2);
        CHECK(same_hits(next, render::lay_out(scatter_scene(N / 2, 42), {}, nullptr, &cache)) > 50);
    }

    SECTION("a run of small changes in both directions") {
        auto scene = render::lay_out(scatter_scene(N, 42), {}, nullptr, &cache);
        std::size_t count = N;
        const int steps[] = {+1, +5, -3, +40, -12, -1, +7, -60, +2, +30};
        for (int s : steps) {
            count = std::size_t(int(count) + s);
            scene = render::lay_out(std::move(scene), scatter_scene(count, 42), {}, nullptr, &cache, &u);
            CHECK(u.reused);
            CHECK_FALSE(u.rebuilt);
            CHECK(same_hits(scene, render::lay_out(scatter_scene(count, 42), {}, nullptr, &cache)) > 50);
        }
    }
}

TEST_CASE("a different quadric table is not reused", "[render][update]") {
    // A different cylinder: same count, different closed form.
    bd::Trace<double> tr;
    auto torus = tr.torus(2.0, 1.0, 48, 24);
    auto item = tr.cylinder(0.02, 0.1, 8, 2);
    auto sc = tr.scatter(item, torus, 500, 42, 0.45, bd::SeatAxis::X);
    auto root = tr.compose({torus, sc});
    auto other = bd::cook(tr, root.index, 0.0);

    render::InstanceUpdate u;
    auto prev = render::lay_out(scatter_scene(500, 42));
    auto next = render::lay_out<double>(std::move(prev), other, {}, nullptr, nullptr, &u);
    CHECK(u.rebuilt);
    CHECK(same_hits(next, render::lay_out(other)) > 20);
}
