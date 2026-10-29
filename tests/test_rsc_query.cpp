// Questions with a free symbol, answered by synthesis.
//
// Each claim is held against a second route: a found answer against the closed form the
// question is known to have, at parameters it was not fitted on. Each check gets a case it must
// fail: a root asked for on a bracket where the function does not change sign is refused, and an
// integral whose closed form needs a function the operations do not have takes pieces.

#include <rsc/spatium/query.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <numbers>
#include <random>

using Catch::Matchers::WithinAbs;
using rsc::Expr;

namespace {

// The answer's value at a parameter: the piece that holds it.
double at(const rsc::Answer& a, double p) {
    for (const auto& s : a.pieces)
        if (p >= s.lo && p <= s.hi)
            return s.expr(p, 0.0);
    return std::numeric_limits<double>::quiet_NaN();
}

void holds(const rsc::Answer& a, double lo, double hi, auto&& truth, double tol) {
    std::mt19937 rng(3);
    std::uniform_real_distribution<double> u(lo, hi);
    for (int i = 0; i < 200; ++i) {
        const double p = u(rng);
        REQUIRE_THAT(at(a, p), WithinAbs(truth(p), tol * std::abs(truth(p))));
    }
}

} // namespace

TEST_CASE("query: the root of x^2 - p is sqrt(p)", "[rsc][query]") {
    const Expr x = Expr::u(), p = Expr::v();
    const rsc::Query q{.ask = rsc::Ask::Solve, .f = x * x - p, .x_lo = 0, .x_hi = 10, .p_lo = 0.5, .p_hi = 4,
                       .tolerance = 1e-10};
    const auto a = rsc::answer(q);
    REQUIRE(a.answered);
    CHECK(a.pieces.size() == 1);
    CHECK(a.pieces.front().terms == 1);
    holds(a, 0.5, 4, [](double v) { return std::sqrt(v); }, 1e-10);
}

TEST_CASE("query: the integral of x^p over [0, 1] is 1 / (p + 1)", "[rsc][query]") {
    // x^p written as exp(p log x): the question needs the logarithm and the exponential of the IR,
    // the answer does not
    const Expr x = Expr::u(), p = Expr::v();
    const rsc::Query q{.ask = rsc::Ask::Integrate, .f = exp(p * log(x)), .x_lo = 0, .x_hi = 1, .p_lo = 0,
                       .p_hi = 3, .tolerance = 1e-9};
    const auto a = rsc::answer(q);
    REQUIRE(a.answered);
    CHECK(a.pieces.size() == 1);
    holds(a, 0, 3, [](double v) { return 1 / (v + 1); }, 1e-9);
}

TEST_CASE("query: the integral of exp(-p x^2) over [0, 1] needs erf, and takes pieces", "[rsc][query]") {
    // sqrt(pi / p) erf(sqrt p) / 2: erf is not among the operations; the answer is still certified
    const Expr x = Expr::u(), p = Expr::v();
    const rsc::Query q{.ask = rsc::Ask::Integrate, .f = exp(Expr(-1.0) * p * x * x), .x_lo = 0, .x_hi = 1,
                       .p_lo = 0.1, .p_hi = 4, .tolerance = 1e-6};
    const auto a = rsc::answer(q);
    INFO(a.pieces.size() << " pieces");
    REQUIRE(a.answered);
    CHECK(a.pieces.size() > 1);
    holds(a, 0.1, 4, [](double v) { return std::sqrt(std::numbers::pi / v) * std::erf(std::sqrt(v)) / 2; }, 1e-6);
}

TEST_CASE("query: a root on a bracket with no change of sign is refused, not fitted", "[rsc][query]") {
    const Expr x = Expr::u(), p = Expr::v();
    // x^2 + p > 0 on [0, 10] for every p in [0.5, 4]
    const auto a = rsc::answer({.ask = rsc::Ask::Solve, .f = x * x + p, .x_lo = 0, .x_hi = 10, .p_lo = 0.5, .p_hi = 4});
    CHECK_FALSE(a.answered);
    CHECK(a.pieces.empty());
    CHECK(a.refused.find("no value") != std::string::npos);
}
