// The derivative of a field, as a field: each operation by its rule, a function of the language by
// the rule its row of the table writes. Each rule is held against a dual number, which computes the
// same derivative from the function by its own rules (algebra/dual.hpp), both in intervals; and the
// check is given a wrong derivative, which it must refuse. In a document, a derivative as a
// function, at a point after `at`, and at the value chosen for its parameter are the same number.

#include <catch2/catch_test_macros.hpp>
#include <spatium/answer/evaluate.hpp>
#include <spatium/answer/respond.hpp>
#include <spatium/io/document.hpp>
#include <spatium/io/field_derivative.hpp>
#include <spatium/io/field_parse.hpp>
#include <spatium/io/field_text.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace {

namespace bd = spatium::io::build;
using spatium::answer::Expr;

Expr read(std::string_view text) {
    auto f = bd::parse(text, "", "");
    REQUIRE(f);
    return *f;
}

std::vector<double> grid(double lo, double hi) {
    std::vector<double> out;
    for (int k = 0; k < 16; ++k) out.push_back(lo + (hi - lo) * (k + 0.5) / 16);
    return out;
}

bd::Document doc(std::string_view text) {
    auto d = bd::read_document(text);
    INFO(text);
    REQUIRE(d);
    return *d;
}

} // namespace

TEST_CASE("derivative: every function of the language, by its row's rule, against a dual number", "[io][derivative]") {
    for (const bd::Function& fn : bd::functions) {
        INFO(fn.name);
        for (int k = 0; k < fn.arity; ++k) REQUIRE(bd::parse(fn.derivative[static_cast<std::size_t>(k)], "", ""));
        // each argument a function of x, so the chain rule runs through every argument
        const std::string text = fn.arity == 1 ? std::format("{}(x^2 + 1)", fn.name) : std::format("{}(x, 2x - 1)", fn.name);
        const Expr f = read(text);
        const auto df = bd::derivative(f, "x");
        REQUIRE(df);
        const auto agreed = spatium::answer::detail::agree_with_dual(f, *df, "x", grid(-3, 3));
        INFO(*bd::to_text(*df) << " " << agreed.refused);
        CHECK(agreed.compared >= 15);   // `less`, min and max may tie at one point
        CHECK_FALSE(agreed.disagrees);
    }
}

TEST_CASE("derivative: sums, products, quotients and powers, against a dual number", "[io][derivative]") {
    for (const char* text : {"x^3 - 2x + 1", "sin(x)^2 * exp(-x)", "(x^2 + 1) / (x - 4)", "x^0.5 * log(x)", "2^x", "x^x",
                             "sqrt(abs(x)) - cos(3x)", "1 / (1 + exp(-x))", "max(sin(x), cos(x)) * x", "log(x^2 + 1) / x"}) {
        INFO(text);
        const Expr f = read(text);
        const auto df = bd::derivative(f, "x");
        REQUIRE(df);
        const auto agreed = spatium::answer::detail::agree_with_dual(f, *df, "x", grid(0.1, 3.7));
        INFO(*bd::to_text(*df));
        CHECK(agreed.compared == 16);
        CHECK_FALSE(agreed.disagrees);
    }
}

TEST_CASE("derivative: the check refuses a wrong derivative", "[io][derivative]") {
    // the case it must fail: sin's derivative written as itself, and one off by a constant factor
    const Expr f = read("sin(x) * x");
    for (const char* wrong : {"sin(x) * x", "2 * (cos(x) * x + sin(x))"}) {
        INFO(wrong);
        const auto agreed = spatium::answer::detail::agree_with_dual(f, read(wrong), "x", grid(-2, 2));
        CHECK(agreed.disagrees);
    }
}

TEST_CASE("derivative: what does not read x adds no term, and constants stay as written", "[io][derivative]") {
    CHECK(*bd::to_text(*bd::derivative(read("2 * x"), "x")) == "2");
    CHECK(*bd::to_text(*bd::derivative(read("a * y + 3"), "x")) == "0");
    CHECK(*bd::to_text(*bd::derivative(read("a * sin(x)"), "x")) == "a * cos(x)");
    // 1/3 is a quotient, not the double nearest it: in fractions the derivative at 1 is exactly 1
    const auto e = spatium::answer::evaluate(doc("g = d/dx x^3 / 3 at 1"), "g");
    REQUIRE(e.exact);
    CHECK(e.exact->text() == "1");
}

TEST_CASE("derivative: noise and a net that read x are refused by name", "[io][derivative]") {
    auto net = Expr::net_of({0.5, -0.25, 1.0, 0.3}, 1, {Expr::symbol("x")});
    REQUIRE(net);
    const auto d = bd::derivative(*net, "x");
    REQUIRE_FALSE(d);
    CHECK(d.error().message.find("Net") != std::string::npos);
    // a net that does not read x is a constant in it
    CHECK(bd::derivative(*net, "y"));
}

TEST_CASE("derivative: as a function, at a point, and at its parameter's chosen value, one number", "[io][derivative][evaluate]") {
    using spatium::answer::Route;
    const auto function = spatium::answer::evaluate(doc("g = d/dx sin(x)^2 * x"), "g");
    REQUIRE(function.taken);
    CHECK(function.route() == Route::Derivative);
    CHECK(function.values.empty());
    REQUIRE_FALSE(function.candidates.front().checks.empty());
    const double at = 0.7;
    const double by_function = function.written.bind("x", Expr(at))(0.0, 0.0);
    const auto point = spatium::answer::evaluate(doc("g = d/dx sin(x)^2 * x at 0.7"), "g");
    const auto chosen = spatium::answer::evaluate(doc("x in [0, 2]\ng = d/dx sin(x)^2 * x"), "g", {{"x", at}});
    REQUIRE(point.taken);
    REQUIRE(chosen.taken);
    CHECK(point.values.front() == by_function);
    CHECK(chosen.values.front() == by_function);
    CHECK(std::abs(by_function - (std::sin(2 * at) * at + std::sin(at) * std::sin(at))) <= 1e-15);
    // a value read by the derivative is put in before it is taken
    const auto read_value = spatium::answer::evaluate(doc("a = 3\ng = d/dx a * x^2 at 2"), "g");
    REQUIRE(read_value.taken);
    CHECK(read_value.values.front() == 12);
}

TEST_CASE("derivative: what is not asked well is refused, by name", "[io][derivative][document]") {
    // a point for a parameter, a variable that is a value, no variable, no value at the point
    CHECK_FALSE(bd::read_document("p in [0, 1]\ng = d/dp p^2 at 0.5"));
    CHECK_FALSE(bd::read_document("a = 2\ng = d/da a^2"));
    CHECK_FALSE(bd::read_document("at = 2"));
    const auto none = spatium::answer::evaluate(doc("g = d/dx log(x) at 0"), "g");
    CHECK_FALSE(none.taken);
    // a name d divided by another is still a value, not a derivative
    const auto value = bd::read_document("d = 2\nq = 1\nr = d/q");
    REQUIRE(value);
    CHECK(value->statements.back().kind == bd::StatementKind::Value);
}

TEST_CASE("derivative: the rows of a formula editor, as a function and at a point", "[io][derivative][mathjson]") {
    const auto function = spatium::answer::respond_rows(R"({"rows": [["D", ["Power", "x", 3], "x"]]})");
    CHECK(function.string_or("route", "") == "derivative");
    const auto point = spatium::answer::respond_rows(
        R"({"rows": [["Delimiter", ["Sequence", ["D", ["Power", "x", 3], "x"], ["Equal", "x", 2]], {"str": ","}]]})");
    INFO(point.dump());
    CHECK(point.number_or("value", 0) == 12);
    // the point given for another name than the derivative's variable
    const auto other = spatium::answer::respond_rows(
        R"({"rows": [["Delimiter", ["Sequence", ["D", ["Power", "x", 3], "x"], ["Equal", "y", 2]], {"str": ","}]]})");
    CHECK(other.find("refused"));
    // a root's row, which also starts as a Sequence, is still a root
    const auto root = spatium::answer::respond_rows(
        R"({"rows": [["Delimiter", ["Sequence", ["Equal", ["Power", "x", 2], 2], ["Element", "x", ["Interval", 0, 2]]], {"str": ","}]]})");
    CHECK(root.string_or("route", "") != "derivative");
    CHECK_FALSE(root.find("refused"));
}
