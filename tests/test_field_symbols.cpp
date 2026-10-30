// Named symbols in a field, and a field written as text and as LaTeX.
//
// Each claim is held against a second route: a symbol bound to a number or to an input is the
// same expression, by canonical key and bit for bit, as the field written with that number or
// input in its place. Each check gets a case it must fail: an unbound symbol evaluates to NaN
// and does not lower, and naming a symbol the field does not read changes nothing.

#include <catch2/catch_test_macros.hpp>
#include <spatium/io/field.hpp>
#include <spatium/io/field_glsl.hpp>
#include <spatium/io/field_pod.hpp>
#include <spatium/io/field_text.hpp>

#include <bit>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace bd = spatium::io::build;
using F = bd::Field<double>;

namespace {

std::vector<std::uint64_t> key(const F& f) {
    std::vector<std::uint64_t> k;
    REQUIRE(f.canonical_key(k));
    return k;
}

std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }

} // namespace

TEST_CASE("field symbols: a symbol is open until bound, and bound it is the field written out", "[field][symbol]") {
    const F r = F::symbol("r"), x = F::symbol("x");
    const F sphere = x * x + F::u() * F::u() - r * r;

    CHECK(sphere.free_symbols() == std::vector<std::string>{"r", "x"});
    CHECK(std::isnan(sphere(0.5, 0.0)));
    const auto refused = bd::lower(sphere);
    REQUIRE_FALSE(refused.has_value());
    CHECK(refused.error().message.find("symbol 'x'") != std::string::npos);   // the first one it meets

    // bound to a number and to an input: the same expression as the one written with them
    const F bound = sphere.bind("r", F(2.0)).bind("x", F::v());
    const F written = F::v() * F::v() + F::u() * F::u() - F(2.0) * F(2.0);
    CHECK(bound.free_symbols().empty());
    CHECK(key(bound) == key(written));
    CHECK(bound.content_id() == written.content_id());
    for (double u : {-1.5, 0.0, 0.25, 3.0})
        for (double v : {-2.0, 0.5, 1.0})
            CHECK(bits(bound(u, v)) == bits(written(u, v)));
    REQUIRE(bd::lower(bound).has_value());

    // binding merges: x * u with x bound to u holds u once
    const F square = (x * F::u()).bind("x", F::u());
    CHECK(square.size() == (F::u() * F::u()).size());

    // two symbols of one name are one op; of two names, two
    CHECK((r + r).size() == 2);
    CHECK((r + x).size() == 3);
    CHECK(r.content_id() != x.content_id());

    // the case it must fail: a name the field does not read changes nothing
    CHECK(key(sphere.bind("q", F(1.0))) == key(sphere));
    CHECK(sphere.bind("q", F(1.0)).free_symbols() == sphere.free_symbols());
}

TEST_CASE("field symbols: a symbol at the root, bound to a field already in the pool, stays the root", "[field][symbol]") {
    // bind rebuilds by appending, and a root that lands on an earlier op is copied to the end
    const F f = F::symbol("a");
    const F g = f.bind("a", F::u() + F(1.0));
    CHECK(g(2.0, 0.0) == 3.0);
    CHECK(g.topologically_ordered());
    CHECK(F::symbol("a").bind("a", F::u())(5.0, 0.0) == 5.0);
}

TEST_CASE("field text: precedence decides the parentheses, and symbols keep their names", "[field][text]") {
    const F u = F::u(), v = F::v();
    CHECK(*bd::to_text(sqrt(u) + F(2.0) * log(v), "x", "y") == "sqrt(x) + 2 * log(y)");
    // the case the old printer got wrong: u / v * 2 reads as (u / v) * 2
    CHECK(*bd::to_text(u / (v * F(2.0))) == "u / (v * 2)");
    CHECK(*bd::to_text(u / v * F(2.0)) == "u / v * 2");
    CHECK(*bd::to_text((u + v) * (u - v)) == "(u + v) * (u - v)");
    CHECK(*bd::to_text(u - (v - F(1.0))) == "u - (v - 1)");
    CHECK(*bd::to_text(u + F(-0.5)) == "u - 0.5");
    CHECK(*bd::to_text(u * F(-2.0)) == "u * (-2)");
    // a factor of exactly 1 is left out, which changes no value: x * 1 is x to the last bit
    CHECK(*bd::to_text(F(1.0) * sqrt(u)) == "sqrt(u)");
    CHECK(*bd::to_text(u * F(1.0) + v) == "u + v");
    CHECK(*bd::to_latex(F(0.5) * u + F(1.25) * F(1.0)) == "0.5 \\cdot u + 1.25");
    // the case it must fail: a factor near 1 is a number, not a 1
    CHECK(*bd::to_text(F(1.0000000000000002) * u) == "1.0000000000000002 * u");
    CHECK(*bd::to_text(F::symbol("mass") * F(2.0) / F::symbol("r")) == "mass * 2 / r");
    CHECK_FALSE(bd::to_text(F([](double, double) { return 0.0; })).has_value());
}

TEST_CASE("field latex: the same walk, spelled for a page", "[field][latex]") {
    const F u = F::u(), v = F::v();
    CHECK(*bd::to_latex(sqrt(u) + F(2.0) * log(v), "x", "y") == "\\sqrt{x} + 2 \\cdot \\ln\\left(y\\right)");
    CHECK(*bd::to_latex(u / (v * F(2.0))) == "\\frac{u}{v \\cdot 2}");
    CHECK(*bd::to_latex((u + v) * F(3.0)) == "\\left(u + v\\right) \\cdot 3");
    CHECK(*bd::to_latex(F::symbol("theta") + F::symbol("r_0") + F::symbol("mass")) ==
          "\\theta + r_{0} + \\mathrm{mass}");
    CHECK(*bd::to_latex(F(1e-26) * u) == "10^{-26} \\cdot u");
    CHECK(*bd::to_latex(F(1.8e-5) * u) == "1.8 \\cdot 10^{-5} \\cdot u");
    CHECK(*bd::to_latex(F(-2.5e20) + u) == "-2.5 \\cdot 10^{20} + u");
    CHECK(*bd::to_latex(less(u, v)) == "\\left[u < v\\right]");
    CHECK(*bd::to_latex(F::coord(1) * F::coord(1)) == "x^{1} \\cdot x^{1}");
    CHECK_FALSE(bd::to_latex(F([](double, double) { return 0.0; })).has_value());
}
