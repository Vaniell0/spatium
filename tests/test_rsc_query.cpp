// Questions with a free symbol, answered by synthesis.
//
// Each claim is held against a second route: a found answer against the closed form the
// question is known to have, at parameters it was not fitted on. Each check gets a case it must
// fail: a root asked for on a bracket where the function does not change sign is refused, and an
// integral whose closed form needs a function the operations do not have takes pieces.

#include <rsc/spatium/query.hpp>
#include <spatium/io/field_text.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <bit>
#include <cmath>
#include <cstdint>
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

TEST_CASE("query: the integral of exp(-p x^2) over [0, 1] needs erf, and says it found no closed form", "[rsc][query]") {
    // sqrt(pi / p) erf(sqrt p) / 2: erf is not among the operations; the answer is still certified,
    // as a Chebyshev polynomial, and marked so -- or, without that fallback, as pieces of short terms
    const Expr x = Expr::u(), p = Expr::v();
    rsc::Query q{.ask = rsc::Ask::Integrate, .f = exp(Expr(-1.0) * p * x * x), .x_lo = 0, .x_hi = 1,
                 .p_lo = 0.1, .p_hi = 4, .tolerance = 1e-6};
    const auto a = rsc::answer(q);
    INFO(a.pieces.size() << " pieces");
    REQUIRE(a.answered);
    CHECK(std::ranges::all_of(a.pieces, [](const rsc::Solution& s) { return s.form == rsc::Form::Chebyshev; }));
    holds(a, 0.1, 4, [](double v) { return std::sqrt(std::numbers::pi / v) * std::erf(std::sqrt(v)) / 2; }, 1e-6);
    q.polynomial_degree = 0;
    const auto short_terms = rsc::answer(q);
    REQUIRE(short_terms.answered);
    CHECK(short_terms.pieces.size() > a.pieces.size());
    holds(short_terms, 0.1, 4, [](double v) { return std::sqrt(std::numbers::pi / v) * std::erf(std::sqrt(v)) / 2; }, 1e-6);
}

TEST_CASE("query: a root on a bracket with no change of sign is refused, not fitted", "[rsc][query]") {
    const Expr x = Expr::u(), p = Expr::v();
    // x^2 + p > 0 on [0, 10] for every p in [0.5, 4]
    const auto a = rsc::answer({.ask = rsc::Ask::Solve, .f = x * x + p, .x_lo = 0, .x_hi = 10, .p_lo = 0.5, .p_hi = 4});
    CHECK_FALSE(a.answered);
    CHECK(a.pieces.empty());
    CHECK(a.refused.find("no value") != std::string::npos);
}

TEST_CASE("query: written in named symbols, the question is the one written in U and V", "[rsc][query]") {
    // the same root asked twice; the symbols are bound to U and V before anything is evaluated,
    // so the two answers are one answer to the last bit
    const Expr x = Expr::symbol("x"), p = Expr::symbol("p");
    const rsc::Query named{.ask = rsc::Ask::Solve, .f = x * x - p, .unknown = "x", .parameter = "p", .x_lo = 0,
                           .x_hi = 10, .p_lo = 0.5, .p_hi = 4, .tolerance = 1e-10};
    const rsc::Query positional{.ask = rsc::Ask::Solve, .f = Expr::u() * Expr::u() - Expr::v(), .x_lo = 0,
                                .x_hi = 10, .p_lo = 0.5, .p_hi = 4, .tolerance = 1e-10};
    const auto a = rsc::answer(named), b = rsc::answer(positional);
    REQUIRE(a.answered);
    REQUIRE(a.pieces.size() == b.pieces.size());
    for (double v : {0.5, 1.0, 2.25, 3.7})
        CHECK(std::bit_cast<std::uint64_t>(at(a, v)) == std::bit_cast<std::uint64_t>(at(b, v)));
    // and read under the parameter's name, in code and on a page
    CHECK(spatium::io::build::to_text(a.pieces.front().expr, named.parameter)->find("sqrt(p)") != std::string::npos);
    CHECK(spatium::io::build::to_latex(a.pieces.front().expr, named.parameter)->find("\\sqrt{p}") != std::string::npos);
}

TEST_CASE("query: a symbol the question reads but does not name is refused by name", "[rsc][query]") {
    const Expr x = Expr::symbol("x"), p = Expr::symbol("p"), a = Expr::symbol("a");
    const auto r = rsc::answer({.ask = rsc::Ask::Solve, .f = x * x - a * p, .unknown = "x", .parameter = "p",
                                .x_lo = 0, .x_hi = 10, .p_lo = 0.5, .p_hi = 4});
    CHECK_FALSE(r.answered);
    CHECK(r.pieces.empty());
    CHECK(r.refused.find("'a'") != std::string::npos);
}

TEST_CASE("query: a question typed as text is the question built in code", "[rsc][query]") {
    // the text names the unknown and the parameter; the answer is the one asked in code with the
    // same symbols, to the last bit
    auto asked = rsc::ask("solve x^2 - p = 0 for x in [0, 10], p in [0.5, 4]");
    REQUIRE(asked.has_value());
    CHECK(asked->ask == rsc::Ask::Solve);
    CHECK(asked->unknown == "x");
    CHECK(asked->parameter == "p");
    CHECK((asked->x_lo == 0 && asked->x_hi == 10 && asked->p_lo == 0.5 && asked->p_hi == 4));
    CHECK(*spatium::io::build::to_text(asked->f) == "x^2 - p");
    asked->tolerance = 1e-10;
    const Expr x = Expr::symbol("x"), p = Expr::symbol("p");
    const rsc::Query built{.ask = rsc::Ask::Solve, .f = x * x - p, .unknown = "x", .parameter = "p", .x_lo = 0,
                           .x_hi = 10, .p_lo = 0.5, .p_hi = 4, .tolerance = 1e-10};
    const auto a = rsc::answer(*asked), b = rsc::answer(built);
    REQUIRE(a.answered);
    REQUIRE(a.pieces.size() == b.pieces.size());
    for (double v : {0.5, 1.0, 2.25, 3.7})
        CHECK(std::bit_cast<std::uint64_t>(at(a, v)) == std::bit_cast<std::uint64_t>(at(b, v)));

    // an integral, an equation with a right side, and a bound that is an expression in no name
    const auto in = rsc::ask("integrate exp(p * log(x)) for x in [0, 1], p in [0, 4]");
    REQUIRE(in.has_value());
    CHECK(in->ask == rsc::Ask::Integrate);
    CHECK((in->unknown == "x" && in->parameter == "p"));
    const auto eq = rsc::ask("solve x * x = p for x in [0, sqrt(2) * 10], p in [1 / 3, 4]");
    REQUIRE(eq.has_value());
    CHECK(*spatium::io::build::to_text(eq->f) == "x^2 - p");
    CHECK(eq->x_hi == std::sqrt(2.0) * 10);
    CHECK(eq->p_lo == 1.0 / 3);
}

TEST_CASE("query: a text that is not a question is refused, and says where", "[rsc][query]") {
    const auto refused = [](std::string_view text) {
        const auto q = rsc::ask(text);
        INFO(text);
        REQUIRE_FALSE(q.has_value());
        return q.error().message;
    };
    CHECK(refused("find x^2 - p for x in [0, 1], p in [0, 1]").find("'solve' or 'integrate'") != std::string::npos);
    CHECK(refused("solve x^2 - p over x in [0, 1], p in [0, 1]").find("'=' or 'for', found 'over'") != std::string::npos);
    CHECK(refused("integrate x = p for x in [0, 1], p in [0, 1]").find("'for', found '='") != std::string::npos);
    CHECK(refused("solve x - p for x in [0, 1]").find("the parameter the answer is a function of") != std::string::npos);
    CHECK(refused("solve x - p for x in [1, 0], p in [0, 1]") == "ask: x in [1, 0] is empty");
    CHECK(refused("solve x - p for x in [0, a], p in [0, 1]").find("x reads 'a', which no line above defines") != std::string::npos);
    CHECK(refused("solve x - p for x in [0, 1], x in [0, 1]").find("named twice") != std::string::npos);
    CHECK(refused("solve x - p for x in [0, 1], p in [0, 1] please").find("the end, found 'please'") != std::string::npos);
    CHECK(refused("solve x - for x in [0, 1], p in [0, 1]").find("parse: at column 11") != std::string::npos);
}

TEST_CASE("query: a limit that moves with the parameter, and special functions as integrals", "[rsc][query]") {
    // Si(x) = integral of sin(t)/t over [0, x]: the parameter is the upper limit. Its reference,
    // the library's integral, is held against the series sum (-1)^n x^(2n+1) / ((2n+1)(2n+1)!)
    const auto si = rsc::ask("integrate sin(t)/t for t in [0, x], x in [0.5, 6]");
    REQUIRE(si.has_value());
    CHECK((si->x_lo == 0 && si->x_hi_at.has_value() && !si->x_lo_at.has_value()));
    const auto series = [](double x) {
        double sum = 0, term = x;   // x^(2n+1) / (2n+1)!
        for (int n = 0; n < 40; ++n) {
            sum += (n % 2 ? -term : term) / (2 * n + 1);
            term *= x * x / ((2 * n + 2) * (2 * n + 3));
        }
        return sum;
    };
    const rsc::Query si_bound = rsc::bound(*si);
    for (double x : {0.5, 1.0, 2.0, 4.0, 6.0})
        CHECK_THAT(rsc::reference(si_bound, x), WithinAbs(series(x), 1e-13));
    auto asked = *si;
    asked.tolerance = 1e-8;
    const auto a = rsc::answer(asked);
    REQUIRE(a.answered);
    holds(a, 0.5, 6, series, 1e-8);

    // erf(x) = 2 / sqrt(pi) times the integral of exp(-t^2) over [0, x], held against std::erf
    const auto erf = rsc::ask("integrate 2 / sqrt(pi) * exp(-t^2) for t in [0, x], x in [0.25, 3]");
    REQUIRE(erf.has_value());
    const rsc::Query erf_bound = rsc::bound(*erf);
    for (double x : {0.25, 0.5, 1.0, 2.0, 3.0})
        CHECK_THAT(rsc::reference(erf_bound, x), WithinAbs(std::erf(x), 1e-15));

    // the Fresnel integral C(x), held against its series sum (-1)^n (pi/2)^(2n) x^(4n+1) / ((2n)! (4n+1))
    const auto fresnel = rsc::ask("integrate cos(pi * t^2 / 2) for t in [0, x], x in [0.25, 2]");
    REQUIRE(fresnel.has_value());
    const auto c_series = [](double x) {
        double sum = 0, term = x;   // (pi/2)^(2n) x^(4n+1) / (2n)!
        const double h = std::numbers::pi / 2;
        for (int n = 0; n < 30; ++n) {
            sum += (n % 2 ? -term : term) / (4 * n + 1);
            term *= h * h * x * x * x * x / ((2 * n + 1) * (2 * n + 2));
        }
        return sum;
    };
    const rsc::Query c_bound = rsc::bound(*fresnel);
    for (double x : {0.25, 1.0, 2.0})
        CHECK_THAT(rsc::reference(c_bound, x), WithinAbs(c_series(x), 1e-13));

    // the cases it must fail: a bound reading another name, a parameter's bound reading one, and an
    // interval that is empty somewhere the answer is wanted
    CHECK(rsc::ask("integrate t for t in [0, y], x in [0, 1]").error().message.find("I reads 'y', which no line above defines") !=
          std::string::npos);
    CHECK(rsc::ask("integrate t for t in [0, x], x in [0, y]").error().message.find("an end of the parameter x reads 'y'") !=
          std::string::npos);
    CHECK(rsc::ask("integrate t for t in [1, x], x in [0, 2]").error().message.find("is empty at x = 0") !=
          std::string::npos);
}
