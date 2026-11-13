// A field in its shortest form: like terms collected, one denominator, roots out of their square
// factors, exact values of sines -- each kept only when a page shows no more of it and its values
// hold. The case it must fail: a rewrite that cancels x / x to 1, which has a value at 0 where x / x
// has none, and a candidate that is not the same function.

#include <catch2/catch_test_macros.hpp>
#include <spatium/io/field_derivative.hpp>
#include <spatium/io/field_parse.hpp>
#include <spatium/io/field_shortest.hpp>
#include <spatium/io/field_text.hpp>

#include <string>

namespace {

namespace bd = spatium::io::build;

bd::Field<double> read(std::string_view text) {
    auto f = bd::parse(text, "", "");
    REQUIRE(f);
    return *f;
}

std::string shortest(std::string_view text) { return *bd::to_text(bd::shortest(read(text))); }

std::string shortest_derivative(std::string_view text) { return *bd::to_text(bd::shortest(*bd::derivative(read(text), "x"))); }

} // namespace

TEST_CASE("shortest: like terms collected, one denominator, as a derivative leaves them", "[io][shortest]") {
    CHECK(shortest_derivative("x^3") == "3 * x^2");
    CHECK(shortest_derivative("x^3 / 3") == "x^2");
    CHECK(shortest_derivative("sin(x)^2") == "2 * cos(x) * sin(x)");
    CHECK(shortest_derivative("log(x) / x") == "(1 - log(x)) / x^2");
    CHECK(shortest_derivative("sqrt(x^2 + 1)") == "x / sqrt(x^2 + 1)");
    CHECK(shortest_derivative("1 / (1 + exp(-x))") == "exp(-x) / (exp(-x) + 1)^2");
    CHECK(shortest("x * y - y * x + 3") == "3");
    CHECK(shortest("(x + 1) * (x - 1)") == "x^2 - 1");
    CHECK(shortest("pi / 6 + pi / 3") == "pi / 2");
}

TEST_CASE("shortest: roots out of their square factors, and the exact values of sines", "[io][shortest]") {
    CHECK(shortest("sqrt(24) + sqrt(150)") == "7 * sqrt(6)");
    CHECK(shortest("sqrt(50)") == "sqrt(50)");   // 5 sqrt(2) shows as much and its text is longer
    CHECK(shortest("1 / sqrt(2)") == "sqrt(2) / 2");
    CHECK(shortest("sin(pi / 6)") == "1 / 2");
    CHECK(shortest("cos(pi / 4)") == "sqrt(2) / 2");
    CHECK(*bd::to_latex(bd::shortest(read("sin(4 * pi / 3)"))) == "-\\frac{\\sqrt{3}}{2}");
    CHECK(shortest("8^(1/3)") == "2");
    CHECK(shortest("2^0.5 + 8^(1/3)") == "sqrt(2) + 2");
    CHECK(shortest("x^0.5") == "x^0.5");   // x^(1 / 2) shows as much, and its text is longer
    CHECK(shortest("12^(1/2)") == "2 * sqrt(3)");
}

TEST_CASE("shortest: what a page would show more of, or a value that does not hold, is left as written", "[io][shortest]") {
    CHECK(shortest("(x + 1)^2") == "(x + 1)^2");   // x^2 + 2x + 1 shows more
    CHECK(shortest("x / x") == "x / x");           // 1 has a value at 0, x / x none
    CHECK(shortest("2^x") == "2^x");
    // the check itself refuses a different function, and one with a value where the first has none
    CHECK_FALSE(bd::same_values(read("x + 1"), read("x")));
    CHECK_FALSE(bd::same_values(read("x / x"), read("1")));
    CHECK(bd::same_values(read("(x + 1) * (x - 1)"), read("x^2 - 1")));
    CHECK(bd::shown_length("\\left(x + 1\\right)^{2}") == 6);
    CHECK(bd::shown_length("x^{2} + 2x + 1") == 7);
}
