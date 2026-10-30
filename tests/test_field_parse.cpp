// A field read from text.
//
// The claim is held against the printer it inverts: for a random field, what `to_text` writes
// reads back to the same values bit for bit (special inputs included), prints to the same text
// again, and, where the printer shortened nothing, is the same ops by canonical key. A text
// written by hand is held against the field built in code. Each check gets a case it must fail:
// a text that is not a field is refused at the column where it stops being one.

#include <catch2/catch_test_macros.hpp>
#include <spatium/io/field.hpp>
#include <spatium/io/field_parse.hpp>
#include <spatium/io/field_text.hpp>

#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <random>
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

bool same(double x, double y) { return std::bit_cast<std::uint64_t>(x) == std::bit_cast<std::uint64_t>(y); }

F read(std::string_view text, std::string u = "u", std::string v = "v") {
    const auto f = bd::parse(text, std::move(u), std::move(v));
    INFO(text);
    REQUIRE(f.has_value());
    return *f;
}

std::string refusal(std::string_view text) {
    const auto f = bd::parse(text);
    INFO(text);
    REQUIRE_FALSE(f.has_value());
    return f.error().message;
}

const std::vector<double> kInputs{0.0, -0.0, 1.0, -1.0, 0.5, 2.0, -3.25, 1e300, -1e-300,
                                  std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()};

// A random field over u, v, the symbol a and constants. With `shortened`, it may hold what the
// printer writes more simply -- a factor of exactly 1, a negative constant added -- so only its
// values and its text come back, not its ops.
F grow(std::mt19937& rng, int depth, bool shortened) {
    std::uniform_int_distribution<int> pick(0, depth > 0 ? 15 : 3);
    const auto next = [&] { return grow(rng, depth - 1, shortened); };
    switch (pick(rng)) {
        case 0: return F::u();
        case 1: return F::v();
        case 2: return F::symbol("a");
        case 3: {
            const double plain[] = {2.0, 0.5, 3.25, 1e-5, 1.8e20, 0.1, -0.0, -2.0};
            const double more[] = {1.0, -1.0, -1e-300, 1e300, 0.0, -2.5, 7.0, 1.0 / 3};
            return F((shortened ? more : plain)[std::uniform_int_distribution<int>(0, 7)(rng)]);
        }
        case 4: return next() + next();
        case 5: return next() - next();
        case 6: return next() * next();
        case 7: return next() / next();
        case 8: return sqrt(next());
        case 9: return sin(next());
        case 10: return cos(next());
        case 11: return log(next());
        case 12: return exp(next());
        case 13: return min(next(), next());
        case 14: return max(next(), next());
        default: return less(next(), next());
    }
}

// u + (-2) prints as u - 2, which reads back as a subtraction: the same bits, other ops
bool adds_negative_constant(const F& f) {
    const F g = f.reachable_from_root();
    for (std::size_t i = 0; i < g.size(); ++i) {
        const auto& n = g.op(i);
        if (n.op == bd::Op::Add && g.op(n.b).op == bd::Op::Const && std::signbit(g.op(n.b).value)) return true;
    }
    return false;
}

} // namespace

TEST_CASE("field parse: a text written by hand is the field built in code", "[field][parse]") {
    const F u = F::u(), v = F::v(), x = F::symbol("x"), p = F::symbol("p");
    CHECK(key(read("u + v * 2")) == key(u + v * F(2.0)));
    CHECK(key(read("(u + v) * 2")) == key((u + v) * F(2.0)));
    CHECK(key(read("u - v - 1")) == key(u - v - F(1.0)));
    CHECK(key(read("u - (v - 1)")) == key(u - (v - F(1.0))));
    CHECK(key(read("u / v / 2")) == key(u / v / F(2.0)));
    CHECK(key(read("x^2 - p")) == key(x * x - p));
    CHECK(key(read("x^3")) == key(x * x * x));
    CHECK(key(read("-x^2")) == key(F(-1.0) * (x * x)));
    CHECK(key(read("-2^2")) == key(F(-1.0) * (F(2.0) * F(2.0))));
    CHECK(key(read("2 * -x")) == key(F(2.0) * (F(-1.0) * x)));
    CHECK(key(read("sqrt(p) + 2 * log(x)")) == key(sqrt(p) + F(2.0) * log(x)));
    CHECK(key(read("min(u, max(v, 0.5)) * less(u, 1e-5)")) == key(min(u, max(v, F(0.5))) * less(u, F(1e-5))));
    CHECK(key(read("exp(p * log(x))")) == key(exp(p * log(x))));
    // a product after a number or a closing parenthesis may leave out its '*'
    CHECK(key(read("2x^2 - 3x")) == key(F(2.0) * (x * x) - F(3.0) * x));
    CHECK(key(read("(x - 1)(x + 2)")) == key((x - F(1.0)) * (x + F(2.0))));
    CHECK(key(read("2 (u + v)")) == key(F(2.0) * (u + v)));
    CHECK(key(read("sqrt(x)x")) == key(sqrt(x) * x));
    // the inputs under the names the caller gives; any other name is a symbol
    CHECK(key(read("x * y", "x", "y")) == key(u * v));
    CHECK(key(read("u * w")) == key(u * F::symbol("w")));
    CHECK(read("u * w").free_symbols() == std::vector<std::string>{"w"});

    // a minus before a number is the number: negative zero, -inf, a number in exponent form
    CHECK(read("-2").size() == 1);
    CHECK(same(read("-0")(0.0, 0.0), -0.0));
    CHECK(same(read("-inf")(0.0, 0.0), -std::numeric_limits<double>::infinity()));
    CHECK(std::isnan(read("nan")(0.0, 0.0)));
    CHECK(same(read("1.8e+20")(0.0, 0.0), 1.8e20));
    CHECK(same(read(".5")(0.0, 0.0), 0.5));
    // pi is the double nearest pi, and prints as its name, which reads back to the same bits
    CHECK(same(read("pi")(0.0, 0.0), std::numbers::pi));
    CHECK(*bd::to_text(read("2 * pi - u")) == "2 * pi - u");
    CHECK(*bd::to_latex(read("u - pi")) == "u - \\pi");
}

TEST_CASE("field parse: what to_text writes reads back to the same field", "[field][parse]") {
    std::mt19937 rng(17);
    int same_ops = 0;
    for (int i = 0; i < 400; ++i) {
        const F f = grow(rng, 5, i % 2 == 1);
        const auto text = bd::to_text(f);
        REQUIRE(text.has_value());
        INFO(*text);
        const F g = read(*text);
        CHECK(*bd::to_text(g) == *text);
        // the symbol is read back as the symbol; bound, the two give the same bits
        const F fb = f.bind("a", F(0.75)), gb = g.bind("a", F(0.75));
        for (double a : kInputs)
            for (double b : kInputs) {
                INFO("u = " << a << ", v = " << b);
                CHECK(same(gb(a, b), fb(a, b)));
            }
        bool plain = !adds_negative_constant(f);
        const F r = f.reachable_from_root();
        for (std::size_t j = 0; j < r.size(); ++j)
            if (r.op(j).op == bd::Op::Const && r.op(j).value == 1.0) plain = false;
        if (plain) {
            CHECK(key(g) == key(f));
            ++same_ops;
        }
    }
    CHECK(same_ops >= 200);

    // the printer keeps the parentheses a sum and a product need to read back in the same order:
    // floating-point addition does not reassociate, (1e300 + -1e300) + 1 is 1 and 1e300 + (-1e300 + 1) is 0
    const F u = F::u(), v = F::v();
    const F right = u + (v + F(1.0));
    CHECK(*bd::to_text(right) == "u + (v + 1)");
    CHECK(*bd::to_text(u * (v * F(2.0))) == "u * (v * 2)");
    CHECK(*bd::to_text(u * (v / F(2.0))) == "u * (v / 2)");
    CHECK(*bd::to_latex(right) == "u + v + 1");
    CHECK(read("u + (v + 1)")(1e300, -1e300) == 0.0);
    CHECK(read("u + v + 1")(1e300, -1e300) == 1.0);
    // and writes a subtraction only for a negative constant added, not for any sum starting with '-'
    CHECK(*bd::to_text(u + (F(-2.0) - v)) == "u + (-2 - v)");
    CHECK(*bd::to_text(u + F(-2.0)) == "u - 2");
}

TEST_CASE("field parse: a text that is not a field is refused where it stops being one", "[field][parse]") {
    CHECK(refusal("u +") == "parse: at column 4, expected a number, a name, a function or '(', found the end");
    CHECK(refusal("u v") == "parse: at column 3, expected an operator or the end, found 'v'");
    CHECK(refusal("(u + v") == "parse: at column 7, expected ')', found the end");
    CHECK(refusal("sin u").find("'sin' is a function: write sin(...)") != std::string::npos);
    CHECK(refusal("tan(u)") == "parse: at column 1, expected sin, cos, log, exp, sqrt, min, max or less "
                               "before '(', found 'tan'");
    CHECK(refusal("min(u)").find("',' between the two arguments of min") != std::string::npos);
    CHECK(refusal("sqrt(u, v)").find("')' after the one argument of sqrt") != std::string::npos);
    CHECK(refusal("u^2.5").find("a whole power from 1 to 64") != std::string::npos);
    CHECK(refusal("u^0").find("column 3") != std::string::npos);
    CHECK(refusal("u^-1").find("a whole power") != std::string::npos);
    CHECK(refusal("u^2^2").find("found '^'") != std::string::npos);
    CHECK(refusal("") == "parse: at column 1, expected a number, a name, a function or '(', found the end");
}
