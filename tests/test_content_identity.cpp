// Content identity: what a field, a motion and a trace node ARE, as a value
// that survives rebuilding them, and that regeneration can compare instead
// of looking at the output.
#include <catch2/catch_test_macros.hpp>
#include <spatium/algebra/noise.hpp>
#include <spatium/io/build.hpp>
#include <spatium/render/cooked_scene.hpp>

#include <cstdint>
#include <memory>
#include <vector>

using namespace spatium;
namespace bd = spatium::io::build;

namespace {

using F = bd::Field<double>;
using V = bd::VecField<double>;

// Every id equal for the life of this object, to show that equality is not
// decided by the hash.
struct CollidingIds {
    CollidingIds() { bd::detail::id_mask = 0; }
    ~CollidingIds() { bd::detail::id_mask = ~std::uint64_t{0}; }
};

F expression(double c) {
    // (u * c + v) * sin(u) + v * 2.5
    return (F::u() * F{c} + F::v()) * sin(F::u()) + F::v() * F{2.5};
}

F stateless_leaf() {
    return F::opaque([](double u, double) { return u * 2.0; });
}

}  // namespace

TEST_CASE("a field's content id does not depend on how it was built", "[field][identity]") {
    SECTION("separately built equal fields") {
        const F a = expression(3.0), b = expression(3.0);
        CHECK(a.content_id() == b.content_id());
        CHECK(bd::same_value(a, b));
    }

    SECTION("equal ops merged or not: a different layout of the same expression") {
        const F e = sin(F::u()) + F::v();
        const F merged = e * e;
        const F spread = e.without_consing() * e;
        REQUIRE(merged.size() < spread.size());
        CHECK(merged.content_id() == spread.content_id());
        CHECK(bd::same_value(merged, spread));
    }

    SECTION("a copy, and a field with the same ops appended to it") {
        F acc = F::u();
        acc += F{2.0};
        const F copy = acc;
        CHECK(copy.content_id() == acc.content_id());
        F longer = acc;
        longer += F::v();
        CHECK(longer.content_id() != acc.content_id());
        // The id of an op already there is what it was.
        CHECK(longer.content_id(acc.size() - 1) == acc.content_id());
    }

    SECTION("different expressions differ") {
        CHECK(expression(3.0).content_id() != expression(4.0).content_id());
        CHECK_FALSE(bd::same_value(expression(3.0), expression(4.0)));
        CHECK_FALSE(bd::same_value(F::u() - F::v(), F::v() - F::u()));
        CHECK_FALSE(bd::same_value(F{0.0}, F{-0.0}));   // equal as numbers, not as bits
    }

    SECTION("a noise table and a gather table count in full") {
        algebra::PerlinNoise n1(1), n2(2);
        auto pts = [](double shift) {
            auto t = std::make_shared<std::vector<Vec<double, 3>>>();
            for (int i = 0; i < 50; ++i) t->push_back(Vec<double, 3>{double(i), double(i) * 2, double(i) + shift});
            return std::shared_ptr<const std::vector<Vec<double, 3>>>(t);
        };
        const F a = noise(n1, F::u(), F::v(), F{0.5}), b = noise(n1, F::u(), F::v(), F{0.5});
        CHECK(bd::same_value(a, b));
        CHECK_FALSE(bd::same_value(a, noise(n2, F::u(), F::v(), F{0.5})));

        // Two tables equal in their first point and size and different in the last.
        const F g1 = gather(pts(0.0), F::u(), 2);
        auto changed = std::make_shared<std::vector<Vec<double, 3>>>(*pts(0.0));
        changed->back()[2] += 1.0;
        const F g2 = gather(changed, F::u(), 2);
        CHECK(g1.content_id() != g2.content_id());
        CHECK_FALSE(bd::same_value(g1, g2));
        CHECK(bd::same_value(g1, gather(pts(0.0), F::u(), 2)));   // equal contents, another table
    }
}

TEST_CASE("changing one leaf changes the ids on the path to the root and nowhere else",
          "[field][identity]") {
    const F a = expression(3.0), b = expression(4.0);
    REQUIRE(a.size() == b.size());

    // Which ops read the constant that differs: it, and everything above it.
    std::size_t leaf = a.size();
    for (std::size_t i = 0; i < a.size(); ++i)
        if (a.op(i).op == bd::Op::Const && a.op(i).value == 3.0) leaf = i;
    REQUIRE(leaf < a.size());
    std::vector<char> above(a.size(), 0);
    above[leaf] = 1;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const auto k = bd::arity(a.op(i).op);
        if ((k > 0 && above[a.op(i).a]) || (k > 1 && above[a.op(i).b]) || (k > 2 && above[a.op(i).c]))
            above[i] = 1;
    }
    REQUIRE(above[a.size() - 1]);   // the root reads it

    std::size_t differ = 0, same = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        CHECK((a.content_id(i) != b.content_id(i)) == bool(above[i]));
        (above[i] ? differ : same) += 1;
    }
    CHECK(differ >= 3);
    CHECK(same >= 3);   // u, v, sin(u), 2.5, v*2.5 are untouched
}

TEST_CASE("an opaque leaf keeps its identity through copies and has a new one when written again",
          "[field][identity]") {
    const double k = 2.0;
    const F a = F::opaque([k](double u, double) { return u * k; });
    const F copy = a;
    const F fresh = F::opaque([k](double u, double) { return u * k; });   // the same text, written again

    CHECK(copy.content_id() == a.content_id());
    CHECK(bd::same_value(copy, a));
    CHECK(a.content_id() != fresh.content_id());
    CHECK_FALSE(bd::same_value(a, fresh));

    // The identity travels into the fields built on it.
    CHECK(bd::same_value(a + F::v(), copy + F::v()));
    CHECK_FALSE(bd::same_value(a + F::v(), fresh + F::v()));

    // A stateless callable is its type: written once in a function, made
    // twice, it is one function; another lambda is another.
    CHECK(bd::same_value(stateless_leaf(), stateless_leaf()));
    CHECK_FALSE(bd::same_value(stateless_leaf(), F::opaque([](double u, double) { return u * 2.0; })));
}

TEST_CASE("a hash collision is a candidate, never a verdict", "[field][identity]") {
    SECTION("fields") {
        CollidingIds all_equal;
        const F a = expression(3.0), b = expression(4.0), a2 = expression(3.0);
        REQUIRE(a.content_id() == b.content_id());   // the ids no longer tell them apart
        CHECK_FALSE(bd::same_value(a, b));
        CHECK(bd::same_value(a, a2));
        CHECK_FALSE(bd::same_value(stateless_leaf(), F::opaque([k = 1.0](double u, double) { return u * k; })));
    }
    SECTION("motions") {
        CollidingIds all_equal;
        const V a = V::translation({0, 0, 1}), b = V::translation({0, 0, 2});
        REQUIRE(a.content_id() == b.content_id());
        CHECK_FALSE(bd::same_value(a, b));
        CHECK(bd::same_value(a, V::translation({0, 0, 1})));
    }
    SECTION("trace nodes") {
        CollidingIds all_equal;
        auto version = [](std::uint32_t seed, std::size_t count) {
            bd::Trace<double> tr;
            auto torus = tr.torus(2.0, 1.0, 48, 24);
            auto item = tr.cylinder(0.011, 0.062, 8, 2);
            auto sc = tr.scatter(item, torus, count, seed, 0.45, bd::SeatAxis::X);
            return *bd::node_version(tr, sc.index, 0.0);
        };
        const auto a = version(42, 100), b = version(43, 100);
        CHECK(a.generator == b.generator);   // collided
        CHECK_FALSE(a.same_generator(b));
        CHECK(a.same_generator(version(42, 150)));
        CHECK_FALSE(a.same_as(version(42, 150)));
    }
    SECTION("regeneration") {
        // A new seed moves every copy. With every id equal the generators
        // look the same to the hash; the scene must still come out as a
        // fresh one does.
        CollidingIds all_equal;
        auto scene = [](std::uint32_t seed) {
            bd::Trace<double> tr;
            auto torus = tr.torus(2.0, 1.0, 48, 24);
            auto item = tr.cylinder(0.011, 0.062, 8, 2);
            auto sc = tr.scatter(item, torus, 400, seed, 0.45, bd::SeatAxis::X);
            auto root = tr.compose({torus, sc});
            return bd::cook(tr, root.index, 0.0);
        };
        auto prev = render::lay_out(scene(42));
        render::InstanceUpdate u;
        auto next = render::lay_out<double>(std::move(prev), scene(43), {}, nullptr, nullptr, &u);
        auto fresh = render::lay_out(scene(43));
        CHECK(u.generators_matched == 0);
        CHECK(u.kept == 0);
        REQUIRE(next.instances.size() == fresh.instances.size());
        for (std::size_t i = 0; i < 20; ++i) {
            // Every one of the new scene's instances is in the result.
            bool found = false;
            for (const auto& s : next.instances.shapes())
                found = found || (s.translation[0] == fresh.instances.shapes()[i].translation[0] &&
                                  s.translation[1] == fresh.instances.shapes()[i].translation[1] &&
                                  s.translation[2] == fresh.instances.shapes()[i].translation[2]);
            CHECK(found);
        }
    }
}

TEST_CASE("a motion's id is its structure", "[field][identity]") {
    SECTION("built twice") {
        auto build = [](double up) {
            return rotated(scaled(V::translation({0, 0, up}), smoothstep(F::t())), F{0.0}, F{0.0}, F::t());
        };
        CHECK(bd::same_value(build(1.0), build(1.0)));
        CHECK_FALSE(bd::same_value(build(1.0), build(2.0)));
        CHECK(build(1.0).content_id() == build(1.0).content_id());
        CHECK(V{}.content_id() == V::point().content_id());   // the identity motion
        CHECK_FALSE(bd::same_value(V{}, V::translation({0, 0, 1})));
    }
    SECTION("a component field is part of it") {
        auto make = [](double c) { return V::make(F{c} * F::t(), F{0.0}, F{0.0}); };
        CHECK(bd::same_value(make(1.0), make(1.0)));
        CHECK_FALSE(bd::same_value(make(1.0), make(2.0)));
    }
    SECTION("callables are opaque leaves") {
        const double k = 1.0;
        auto cap = [k](double t) { return Vec<double, 3>{0, 0, t * k}; };
        const V a = V::opaque_of_time(cap), b = V::opaque_of_time(cap);
        CHECK_FALSE(bd::same_value(a, b));   // two captured states: not known to be one
        auto stateless = [] { return V::opaque_of_time([](double t) { return Vec<double, 3>{0, 0, t}; }); };
        CHECK(bd::same_value(stateless(), stateless()));
    }
}

TEST_CASE("a trace node's version separates the generator from the count", "[build][identity]") {
    struct Made {
        bd::Trace<double> tr;
        std::size_t scatter = 0, torus = 0, chart = 0, group = 0;
    };
    auto make = [](std::size_t count, std::uint32_t seed, double radius, double lift) {
        auto m = std::make_unique<Made>();
        auto torus = m->tr.torus(radius, 1.0, 48, 24);
        auto chart = m->tr.space(make_torus<double>(2.0, 1.0));
        auto item = m->tr.cylinder(0.011, 0.062, 8, 2);
        auto sc = m->tr.scatter(item, torus, count, seed, 0.45, bd::SeatAxis::X);
        if (lift != 0.0) sc.moving(V::translation({0, 0, lift}));
        m->scatter = sc.index;
        m->torus = torus.index;
        m->chart = chart.index;
        m->group = m->tr.compose({torus, sc}).index;
        return m;
    };
    auto version = [](const Made& m, std::size_t node, double t = 0.0) { return bd::node_version(m.tr, node, t); };

    const auto base = make(100, 42, 2.0, 0.0);
    const auto v0 = *version(*base, base->scatter);

    SECTION("equal when built again") {
        const auto again = make(100, 42, 2.0, 0.0);
        CHECK(v0.same_as(*version(*again, again->scatter)));
        CHECK(v0.id == version(*again, again->scatter)->id);
    }
    SECTION("the count is not part of the generator") {
        const auto more = make(101, 42, 2.0, 0.0);
        const auto v = *version(*more, more->scatter);
        CHECK(v0.same_generator(v));
        CHECK(v0.generator == v.generator);
        CHECK(v0.id != v.id);
        CHECK(v.count == 101);
    }
    SECTION("everything else is") {
        for (const auto& other : {make(100, 43, 2.0, 0.0), make(100, 42, 2.5, 0.0), make(100, 42, 2.0, 0.1)}) {
            const auto v = *version(*other, other->scatter);
            CHECK_FALSE(v0.same_generator(v));
            CHECK(v0.generator != v.generator);
        }
        // A motion changes it, and so does the time it is asked at.
        CHECK(version(*base, base->scatter, 0.5)->generator != v0.generator);
    }
    SECTION("a node with nothing to compare it by has no version") {
        CHECK(version(*base, base->chart) == std::nullopt);   // a chart handed to space()
        CHECK(version(*base, base->group) == std::nullopt);   // a group is its children
        CHECK(version(*base, base->torus).has_value());
    }
    SECTION("a plain node's version follows its fields") {
        auto torus_with = [](auto&& decorate) {
            bd::Trace<double> tr;
            auto h = tr.torus(2.0, 1.0, 48, 24);
            decorate(h);
            return *bd::node_version(tr, h.index, 0.0);
        };
        const auto plain = torus_with([](auto&) {});
        CHECK(plain.same_as(torus_with([](auto&) {})));
        CHECK(plain.id != torus_with([](auto& h) { h.moving(V::translation({0, 0, 1})); }).id);
        CHECK(plain.id != torus_with([](auto& h) { h.colored(V::constant({1, 0, 0})); }).id);
        CHECK(plain.id != torus_with([](auto& h) { h.glowing(V::constant({0, 1, 0})); }).id);
        CHECK(plain.id != torus_with([](auto& h) { h.glowing(Vec<double, 3>{0, 1, 0}); }).id);
    }
}

TEST_CASE("cooked objects say where they came from", "[build][identity]") {
    bd::Trace<double> tr;
    auto torus = tr.torus(2.0, 1.0, 48, 24);
    auto item = tr.cylinder(0.011, 0.062, 8, 2);
    auto sc = tr.scatter(item, torus, 50, 42, 0.45, bd::SeatAxis::X);
    auto second = tr.scatter(item, torus, 20, 7, 0.45, bd::SeatAxis::X);
    auto root = tr.compose({torus, sc, second});
    const auto cooked = bd::cook(tr, root.index, 0.0);

    REQUIRE(cooked.generators().size() == 2);
    std::size_t covered = 0;
    for (const auto& g : cooked.generators()) {
        CHECK(g.version.same_as(*bd::node_version(tr, g.node, 0.0)));
        CHECK(g.count == tr.node(g.node).count);
        for (std::size_t s = 0; s < g.count; ++s) {
            const auto& o = cooked.objects()[g.first_object + s];
            CHECK(o.source_node == g.node);
            CHECK(o.site == s);
            CHECK(o.shape == g.shape);
        }
        covered += g.count;
    }
    // The torus is one object of its own, with no site.
    std::size_t sited = 0;
    for (const auto& o : cooked.objects()) sited += o.site != bd::Object<double>::no_site;
    CHECK(sited == covered);
    CHECK(cooked.object_count() == covered + 1);
}
