#include <catch2/catch_test_macros.hpp>
#include <spatium/io/build.hpp>
#include <spatium/render/cooked_scene.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
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

// ── The same scene from a fresh cook, after each kind of edit ────────
//
// A scene regenerated from the previous one must be the scene a fresh cook
// and lay-out make. The edits are the ones an editor makes: a count moved,
// a seed, the surface under the scatter, the scattered shape, the scatter's
// motion, a scatter removed or added.

namespace {

struct Params {
    std::size_t count = 3000;
    std::uint32_t seed = 42;
    double major_r = 2.0;
    double item_r = 0.011;
    bool scatter = true;
    std::size_t second = 0;          // a second scatter of this many, same item
    double lift = 0.0;               // the scatter's motion: a translation up
};

bd::Cooked<double> scene_of(const Params& p) {
    bd::Trace<double> tr;
    auto torus = tr.torus(p.major_r, 1.0, 48, 24);
    auto item = tr.cylinder(p.item_r, 0.062, 8, 2);
    std::vector<bd::Handle<double>> parts{torus};
    if (p.scatter) {
        auto sc = tr.scatter(item, torus, p.count, p.seed, 0.45, bd::SeatAxis::X);
        if (p.lift != 0.0) sc.moving(bd::VecField<double>::translation({0.0, 0.0, p.lift}));
        parts.push_back(sc);
    }
    if (p.second) parts.push_back(tr.scatter(item, torus, p.second, 7, 0.45, bd::SeatAxis::X));
    auto root = tr.compose(parts);
    return bd::cook(tr, root.index, 0.0);
}

// An instance as plain numbers: where it is, how it is turned and scaled,
// what it looks like, and which shape of the table it points at.
using Row = std::array<double, 22>;

std::vector<Row> instance_list(const render::CookedScene<double>& s) {
    std::vector<Row> rows;
    const auto& insts = s.instances.shapes();
    for (std::size_t i = 0; i < insts.size(); ++i) {
        const auto& a = insts[i];
        const auto& l = s.looks[i];
        Row r{};
        std::size_t k = 0;
        for (std::size_t j = 0; j < 3; ++j) r[k++] = a.translation[j];
        r[k++] = a.scale;
        for (std::size_t m = 0; m < 3; ++m)
            for (std::size_t n = 0; n < 3; ++n) r[k++] = a.rotation(m, n);
        for (std::size_t j = 0; j < 3; ++j) r[k++] = l.color[j];
        r[k++] = l.roughness;
        for (std::size_t j = 0; j < 3; ++j) r[k++] = l.emissive[j];
        r[k++] = l.opacity;
        r[k++] = static_cast<double>(a.shape - s.quadrics.data());
        rows.push_back(r);
    }
    std::sort(rows.begin(), rows.end());
    return rows;
}

}  // namespace

TEST_CASE("a regenerated scene is the scene a fresh cook makes", "[render][update]") {
    render::StaticCache<double> cache;
    const Params base;

    struct Edit {
        const char* name;
        Params to;
        bool rebuilt;
        bool update_path;   // the instance tree of the previous scene was the starting point
    };
    auto with = [&](auto f) { Params p = base; f(p); return p; };
    const Edit edits[] = {
        {"count + 1", with([](Params& p) { p.count += 1; }), false, true},
        {"count - 1", with([](Params& p) { p.count -= 1; }), false, true},
        {"count + 100", with([](Params& p) { p.count += 100; }), false, true},
        {"seed", with([](Params& p) { p.seed = 43; }), true, true},
        {"torus radius", with([](Params& p) { p.major_r = 2.5; }), true, true},
        {"scattered shape", with([](Params& p) { p.item_r = 0.02; }), true, false},
        {"motion field", with([](Params& p) { p.lift = 0.05; }), true, true},
        {"scatter removed", with([](Params& p) { p.scatter = false; }), true, false},
        {"second scatter added", with([](Params& p) { p.second = 100; }), false, true},
    };

    for (const auto& e : edits) {
        DYNAMIC_SECTION(e.name) {
            auto prev = render::lay_out(scene_of(base), {}, nullptr, &cache);
            render::InstanceUpdate u;
            const auto cooked = scene_of(e.to);
            auto next = render::lay_out(std::move(prev), cooked, {}, nullptr, &cache, &u);
            auto fresh = render::lay_out(cooked, {}, nullptr, &cache);

            CHECK(u.rebuilt == e.rebuilt);
            CHECK(u.reused == e.update_path);
            if (e.update_path) CHECK(u.hashed == 0);   // no instance was compared with another
            CHECK(next.instances.size() == fresh.instances.size());
            CHECK(next.prims.size() == fresh.prims.size());
            CHECK(instance_list(next) == instance_list(fresh));
            if (fresh.instances.size() > 100) CHECK(same_hits(next, fresh) > 50);
        }
    }
}

TEST_CASE("an update decided by generators touches only the tail", "[render][update]") {
    render::StaticCache<double> cache;
    Params p;
    auto scene = render::lay_out(scene_of(p), {}, nullptr, &cache);

    SECTION("count + 1") {
        p.count += 1;
        render::InstanceUpdate u;
        scene = render::lay_out(std::move(scene), scene_of(p), {}, nullptr, &cache, &u);
        CHECK(u.generators_matched == 1);
        CHECK(u.hashed == 0);
        CHECK(u.kept == 3000);
        CHECK(u.removed == 0);
        CHECK(u.inserted == 1);
    }
    SECTION("count - 100") {
        p.count -= 100;
        render::InstanceUpdate u;
        scene = render::lay_out(std::move(scene), scene_of(p), {}, nullptr, &cache, &u);
        CHECK(u.kept == 2900);
        CHECK(u.removed == 100);
        CHECK(u.inserted == 0);
    }
    SECTION("a seed change is not matched") {
        p.seed = 43;
        render::InstanceUpdate u;
        scene = render::lay_out(std::move(scene), scene_of(p), {}, nullptr, &cache, &u);
        CHECK(u.generators_matched == 0);
        CHECK(u.kept == 0);
        CHECK(u.rebuilt);
    }
    SECTION("a run of changes keeps following the scene") {
        const int steps[] = {+1, +50, -7, +100, -100, -1, +3};
        for (int d : steps) {
            p.count = std::size_t(int(p.count) + d);
            render::InstanceUpdate u;
            scene = render::lay_out(std::move(scene), scene_of(p), {}, nullptr, &cache, &u);
            CHECK(u.hashed == 0);
            CHECK(instance_list(scene) == instance_list(render::lay_out(scene_of(p), {}, nullptr, &cache)));
        }
    }
}

TEST_CASE("instances with no origin are still compared one by one", "[render][update]") {
    // A scatter over a target handed to space() has no version: the chart is
    // a (u,v) map and cannot be compared with another. Its instances are
    // compared by hashing, as before.
    auto build = [](std::size_t count) {
        bd::Trace<double> tr;
        auto chart = tr.space(make_torus<double>(2.0, 1.0), 48, 24);
        auto item = tr.cylinder(0.011, 0.062, 8, 2);
        auto sc = tr.scatter(item, chart, count, 42, 0.45, bd::SeatAxis::X);
        auto root = tr.compose({chart, sc});
        return bd::cook(tr, root.index, 0.0);
    };
    CHECK(build(500).generators().empty());

    render::StaticCache<double> cache;
    auto prev = render::lay_out(build(2000), {}, nullptr, &cache);
    render::InstanceUpdate u;
    auto next = render::lay_out(std::move(prev), build(2001), {}, nullptr, &cache, &u);
    CHECK(u.hashed == 4001);
    CHECK(u.kept == 2000);
    CHECK(u.inserted == 1);
    CHECK(u.generators_matched == 0);
    CHECK(instance_list(next) == instance_list(render::lay_out(build(2001), {}, nullptr, &cache)));
}

// ── Cooking from the previous cook ───────────────────────────────────

namespace {

bool same_object(const bd::Object<double>& a, const bd::Object<double>& b) {
    return std::memcmp(&a.translation, &b.translation, sizeof a.translation) == 0 && a.scale == b.scale &&
           std::memcmp(&a.rotation_q, &b.rotation_q, sizeof a.rotation_q) == 0 &&
           std::memcmp(&a.material, &b.material, sizeof a.material) == 0 && a.shape == b.shape &&
           a.source_node == b.source_node && a.site == b.site && a.instanceable == b.instanceable;
}

bd::Trace<double> trace_of(const Params& p, std::size_t& root) {
    bd::Trace<double> tr;
    auto torus = tr.torus(p.major_r, 1.0, 48, 24);
    auto item = tr.cylinder(p.item_r, 0.062, 8, 2);
    std::vector<bd::Handle<double>> parts{torus};
    if (p.scatter) {
        auto sc = tr.scatter(item, torus, p.count, p.seed, 0.45, bd::SeatAxis::X);
        if (p.lift != 0.0) sc.moving(bd::VecField<double>::translation({0.0, 0.0, p.lift}));
        parts.push_back(sc);
    }
    if (p.second) parts.push_back(tr.scatter(item, torus, p.second, 7, 0.45, bd::SeatAxis::X));
    root = tr.compose(parts).index;
    return tr;
}

}  // namespace

TEST_CASE("a cook that continues the previous one is the cook that starts over", "[build][update]") {
    // Counts either side of the sampler's checkpoint spacing, so a resume
    // starts from a state and not from the beginning.
    const std::size_t every = bd::SiteCheckpoints<double>::every;
    struct Step { const char* name; Params from, to; std::size_t reused; };
    auto p = [](std::size_t count) { Params q; q.count = count; return q; };
    auto edit = [](Params q, auto f) { f(q); return q; };
    const Step steps[] = {
        {"count + 1", p(3000), p(3001), 3000},
        {"count + 1 across a checkpoint", p(every - 1), p(every + 1), every - 1},
        {"count + 100 from many checkpoints", p(5 * every + 3), p(5 * every + 103), 5 * every + 3},
        {"count - 1", p(3000), p(2999), 2999},
        {"count - 100 below a checkpoint", p(3 * every + 50), p(3 * every - 50), 3 * every - 50},
        {"same", p(3000), p(3000), 3000},
        {"seed", p(3000), edit(p(3000), [](Params& q) { q.seed = 43; }), 0},
        {"torus radius", p(3000), edit(p(3000), [](Params& q) { q.major_r = 2.5; }), 0},
        {"motion field", p(3000), edit(p(3000), [](Params& q) { q.lift = 0.05; }), 0},
        {"scattered shape", p(3000), edit(p(3000), [](Params& q) { q.item_r = 0.02; }), 0},
        {"second scatter added", p(3000), edit(p(3000), [](Params& q) { q.second = 100; }), 3000},
        {"scatter removed", p(3000), edit(p(3000), [](Params& q) { q.scatter = false; }), 0},
    };
    for (const auto& st : steps) {
        DYNAMIC_SECTION(st.name) {
            std::size_t r0, r1;
            auto t0 = trace_of(st.from, r0);
            auto t1 = trace_of(st.to, r1);
            const auto before = bd::cook(t0, r0, 0.0);
            const auto fresh = bd::cook(t1, r1, 0.0);
            const auto resumed = bd::cook(t1, r1, 0.0, &before);

            REQUIRE(resumed.object_count() == fresh.object_count());
            REQUIRE(resumed.shape_count() == fresh.shape_count());
            for (std::size_t i = 0; i < fresh.object_count(); ++i)
                REQUIRE(same_object(resumed.objects()[i], fresh.objects()[i]));
            CHECK(resumed.shared_objects() == fresh.shared_objects());
            for (std::size_t i = 0; i < fresh.shapes().size(); ++i)
                CHECK(resumed.shapes()[i].instances == fresh.shapes()[i].instances);

            REQUIRE(resumed.generators().size() == fresh.generators().size());
            std::size_t reused = 0;
            for (std::size_t i = 0; i < fresh.generators().size(); ++i) {
                const auto &a = resumed.generators()[i], &b = fresh.generators()[i];
                CHECK(a.node == b.node);
                CHECK(a.version.same_as(b.version));
                CHECK(a.first_object == b.first_object);
                CHECK(a.count == b.count);
                CHECK(a.shape == b.shape);
                CHECK(b.reused == 0);
                reused += a.reused;
            }
            CHECK(reused == st.reused);
            CHECK(fresh.generators().empty() == !st.to.scatter);

            // And it can be continued again: the states it leaves behind are
            // the ones a fresh cook leaves.
            auto t2 = trace_of(edit(st.to, [](Params& q) { q.count += 7; }), r1);
            const auto again = bd::cook(t2, r1, 0.0, &resumed);
            const auto from_scratch = bd::cook(t2, r1, 0.0);
            REQUIRE(again.object_count() == from_scratch.object_count());
            for (std::size_t i = 0; i < from_scratch.object_count(); ++i)
                REQUIRE(same_object(again.objects()[i], from_scratch.objects()[i]));
        }
    }
}

TEST_CASE("a regenerated scene cooked from the previous cook is laid out like a fresh one", "[render][update]") {
    render::StaticCache<double> cache;
    std::size_t r0, r1;
    Params a, b;
    b.count = a.count + 100;
    auto ta = trace_of(a, r0);
    auto tb = trace_of(b, r1);
    const auto cooked_a = bd::cook(ta, r0, 0.0);
    auto scene = render::lay_out(cooked_a, {}, nullptr, &cache);
    const auto cooked_b = bd::cook(tb, r1, 0.0, &cooked_a);
    render::InstanceUpdate u;
    scene = render::lay_out(std::move(scene), cooked_b, {}, nullptr, &cache, &u);
    CHECK(u.hashed == 0);
    CHECK(u.inserted == 100);
    CHECK(instance_list(scene) == instance_list(render::lay_out(bd::cook(tb, r1, 0.0), {}, nullptr, &cache)));
}

TEST_CASE("a laid-out scene keeps where its instances came from", "[render][update]") {
    Params p;
    p.second = 40;
    const auto scene = render::lay_out(scene_of(p));
    REQUIRE(scene.origins.size() == scene.instances.size());
    REQUIRE(scene.generators.size() == 2);
    std::vector<std::size_t> seen(scene.generators.size(), 0);
    for (const render::InstanceOrigin& o : scene.origins) {
        REQUIRE(o.generator != render::InstanceOrigin::none);
        REQUIRE(o.site < scene.generators[o.generator].count);
        ++seen[o.generator];
    }
    for (std::size_t i = 0; i < seen.size(); ++i) {
        const render::LaidGenerator& g = scene.generators[i];
        CHECK(seen[i] == g.count);
        CHECK_FALSE(g.key.empty());
    }
}

TEST_CASE("sites from the k-th on are the tail of all the sites", "[build][update]") {
    std::size_t root;
    const auto tr = trace_of(Params{.count = 3000}, root);
    std::size_t scatter = 0;
    for (std::size_t i = 0; i < tr.size(); ++i)
        if (tr.node(i).kind == bd::Kind::Scatter) scatter = i;

    const bd::ScatterSites<double> all = bd::scatter_spots_from(tr, scatter, 0.0, 0);
    REQUIRE(all.spots.size() == 3000);
    const auto plain = bd::scatter_spots(tr, scatter, 0.0);
    REQUIRE(plain.size() == 3000);
    CHECK(std::memcmp(plain.data(), all.spots.data(), plain.size() * sizeof(plain[0])) == 0);

    for (std::size_t first : {1ul, 1023ul, 1024ul, 1500ul, 2999ul}) {
        // Resumed from the states the first call left, and from nothing.
        for (const auto* known : {&all.checkpoints, static_cast<const bd::SiteCheckpoints<double>*>(nullptr)}) {
            const auto tail = bd::scatter_spots_from(tr, scatter, 0.0, first, known);
            REQUIRE(tail.spots.size() == 3000 - first);
            CHECK(std::memcmp(tail.spots.data(), all.spots.data() + first,
                              tail.spots.size() * sizeof(plain[0])) == 0);
        }
    }
}

TEST_CASE("a scatter's generator words are what its version is made of", "[build][identity]") {
    std::size_t root;
    const auto tr = trace_of(Params{}, root);
    std::size_t scatter = 0;
    for (std::size_t i = 0; i < tr.size(); ++i)
        if (tr.node(i).kind == bd::Kind::Scatter) scatter = i;
    std::vector<std::uint64_t> words;
    REQUIRE(bd::scatter_generator_words(tr, tr.node(scatter), words));
    const auto version = *bd::node_version(tr, scatter, 0.0);
    // The version's key starts with the node kind, then these words.
    REQUIRE(version.key.size() > words.size());
    CHECK(std::equal(words.begin(), words.end(), version.key.begin() + 1));
}

// ── A new time replaces only what reads the time ─────────────────────

namespace {

// Two scatters over one torus: `first` moves by a constant, `second` by an
// amount that grows with the time. Either can be left out.
bd::Cooked<double> timed_scene(double t, bool fixed, bool timed, const bd::Cooked<double>* previous = nullptr) {
    using F = bd::Field<double>;
    using V = bd::VecField<double>;
    bd::Trace<double> tr;
    auto torus = tr.torus(2.0, 1.0, 48, 24);
    auto item = tr.cylinder(0.011, 0.062, 8, 2);
    std::vector<bd::Handle<double>> parts{torus};
    if (fixed) {
        auto a = tr.scatter(item, torus, 2000, 42, 0.45, bd::SeatAxis::X);
        a.moving(V::translation({0.0, 0.0, 0.05}));
        parts.push_back(a);
    }
    if (timed) {
        auto b = tr.scatter(item, torus, 1500, 7, 0.45, bd::SeatAxis::X);
        b.moving(V::point() + V::make(F{0.0}, F{0.0}, F::t() * F{0.1}));
        parts.push_back(b);
    }
    const auto root = tr.compose(parts).index;
    return previous ? bd::cook(tr, root, t, previous) : bd::cook(tr, root, t);
}

}  // namespace

TEST_CASE("a new time replaces the scatters that read it and no others", "[render][update]") {
    struct Case { const char* name; bool fixed, timed; std::size_t reused; };
    const Case cases[] = {
        {"a scatter that does not read the time", true, false, 2000},
        {"a scatter that reads the time", false, true, 0},
        {"one of each", true, true, 2000},
    };
    for (const auto& c : cases) {
        DYNAMIC_SECTION(c.name) {
            render::StaticCache<double> cache;
            const auto before = timed_scene(0.0, c.fixed, c.timed);
            const auto fresh = timed_scene(2.5, c.fixed, c.timed);
            const auto resumed = timed_scene(2.5, c.fixed, c.timed, &before);

            REQUIRE(resumed.object_count() == fresh.object_count());
            for (std::size_t i = 0; i < fresh.object_count(); ++i)
                REQUIRE(same_object(resumed.objects()[i], fresh.objects()[i]));
            std::size_t reused = 0;
            for (const auto& g : resumed.generators()) reused += g.reused;
            CHECK(reused == c.reused);

            auto scene = render::lay_out(before, {}, nullptr, &cache);
            render::InstanceUpdate u;
            scene = render::lay_out(std::move(scene), resumed, {}, nullptr, &cache, &u);
            const auto expected = render::lay_out(fresh, {}, nullptr, &cache);
            CHECK(instance_list(scene) == instance_list(expected));
            CHECK(scene.prims.size() == expected.prims.size());
            CHECK(same_hits(scene, expected) > 20);
            if (c.reused) CHECK(u.reused);
        }
    }
}

TEST_CASE("a static part laid out before is found by its description, not made again", "[render][update]") {
    render::StaticCache<double> cache;
    Params a;
    const auto first = render::lay_out(scene_of(a), {}, nullptr, &cache);
    CHECK(cache.misses() == 1);

    // Same static part, a different scatter on top: found without emitting.
    Params more = a;
    more.count += 50;
    const auto second = render::lay_out(scene_of(more), {}, nullptr, &cache);
    CHECK(cache.hits() == 1);
    CHECK(cache.misses() == 1);
    REQUIRE(second.prims.size() == first.prims.size());
    for (std::size_t i = 0; i < first.prims.size(); i += 37) {
        REQUIRE(second.prims[i].vertex_normals[0][0] == first.prims[i].vertex_normals[0][0]);
        REQUIRE(second.prims[i].tri.vertices[0][1] == first.prims[i].tri.vertices[0][1]);
    }
    const auto fresh = render::lay_out(scene_of(more));   // no cache: the baseline
    CHECK(same_hits(second, fresh) > 50);

    // A changed static part is a miss, and its triangles are its own.
    Params wider = a;
    wider.major_r = 2.5;
    const auto third = render::lay_out(scene_of(wider), {}, nullptr, &cache);
    CHECK(cache.misses() == 2);
    CHECK(same_hits(third, render::lay_out(scene_of(wider))) > 50);
}
