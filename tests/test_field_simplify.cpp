// Rewrites of a field that change no value.
//
// The claim is held against the field it came from: on every input, special values included
// (+0, -0, infinities, NaN), the simplified field gives the same bits, and it never has more
// ops. Each rewrite that is NOT made is a case it must fail: x + 0, x * 0, x - x and x / x keep
// their ops, because each changes a value somewhere, and x + x keeps its two, since 2 * x is three.

#include <catch2/catch_test_macros.hpp>
#include <spatium/io/build.hpp>
#include <spatium/io/field.hpp>

#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

namespace bd = spatium::io::build;
using F = bd::Field<double>;

namespace {

const std::vector<double> kInputs{0.0, -0.0, 1.0, -1.0, 0.5, 2.0, -3.25, 1e300, -1e-300,
                                  std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(),
                                  std::numeric_limits<double>::quiet_NaN()};

bool same(double x, double y) { return std::bit_cast<std::uint64_t>(x) == std::bit_cast<std::uint64_t>(y); }

void holds(const F& f) {
    const F g = f.simplified();
    CHECK(g.size() <= f.reachable_from_root().size());
    CHECK(g.topologically_ordered());
    for (double u : kInputs)
        for (double v : kInputs) {
            INFO("u = " << u << ", v = " << v);
            CHECK(same(g(u, v), f(u, v)));
        }
}

// a random field over u, v and the constants the rules look for
F grow(std::mt19937& rng, int depth) {
    std::uniform_int_distribution<int> pick(0, depth > 0 ? 13 : 3);
    switch (pick(rng)) {
        case 0: return F::u();
        case 1: return F::v();
        case 2: { const double c[] = {0.0, -0.0, 1.0, 2.0, -1.0, 0.5}; return F(c[std::uniform_int_distribution<int>(0, 5)(rng)]); }
        case 3: return F(1.0);
        case 4: return grow(rng, depth - 1) + grow(rng, depth - 1);
        case 5: return grow(rng, depth - 1) - grow(rng, depth - 1);
        case 6: return grow(rng, depth - 1) * grow(rng, depth - 1);
        case 7: return grow(rng, depth - 1) / grow(rng, depth - 1);
        case 8: { F x = grow(rng, depth - 1); return x + x; }
        case 9: return sqrt(grow(rng, depth - 1));
        case 10: return sin(grow(rng, depth - 1));
        case 11: return min(grow(rng, depth - 1), grow(rng, depth - 1));
        case 12: return less(grow(rng, depth - 1), grow(rng, depth - 1));
        default: return exp(grow(rng, depth - 1)) * F(1.0);
    }
}

} // namespace

TEST_CASE("field simplify: each rewrite gives the same bits, special values included", "[field][simplify]") {
    const F u = F::u(), v = F::v();
    CHECK((F(2.0) * F(2.0) * u).simplified().size() == 3);   // 4, u, the product
    CHECK((u * F(1.0)).simplified().size() == 1);
    CHECK((F(1.0) * u / F(1.0)).simplified().size() == 1);
    CHECK((u - F(0.0)).simplified().size() == 1);
    CHECK((u + F(-0.0)).simplified().size() == 1);
    CHECK((sqrt(F(4.0)) + log(F(1.0))).simplified().size() == 1);
    for (const F& f : {F(2.0) * F(2.0) * u, u * F(1.0), F(1.0) * u / F(1.0), u - F(0.0), u + F(-0.0),
                       sqrt(F(4.0)) + log(F(1.0)) * v, sin(F(0.5)) * cos(v) + exp(F(-1.0))})
        holds(f);

    // the cases it must fail: each of these would change a value, so none is rewritten
    for (const F& f : {u + F(0.0), u * F(0.0), u - u, u / u, u + u}) {
        CHECK(f.simplified().size() == f.reachable_from_root().size());
        holds(f);
    }
}

TEST_CASE("field simplify: random fields keep their values to the bit", "[field][simplify]") {
    std::mt19937 rng(7);
    std::size_t before = 0, after = 0;
    for (int i = 0; i < 300; ++i) {
        const F f = grow(rng, 4);
        holds(f);
        before += f.reachable_from_root().size();
        after += f.simplified().size();
    }
    INFO(before << " ops before, " << after << " after");
    CHECK(after < before);
}

TEST_CASE("field simplify: a name bound to a number in a scene is one constant", "[field][simplify][build_dsl]") {
    bd::Trace<double> trace;
    const auto r = trace.let("r", 0.5);
    const auto torus = trace.torus(2.0, r);
    const F written = F(2.0) * r.field() * r.field();
    const auto shell = trace.offset(torus, written);
    const auto& node = trace.node(shell.index);
    CHECK(node.thickness.size() == 1);
    CHECK(node.thickness(0.0, 0.0) == 0.5);
    CHECK(node.thickness_written.size() == written.size());
}
