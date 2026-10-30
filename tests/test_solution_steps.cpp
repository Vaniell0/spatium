// A quadratic solved with the steps it was solved by, and an equation typed as text solved so.
//
// The steps are the solver's own record, so the roots they end in are the roots `solve_quadratic`
// returns, to the bit; they are held against a second route, Vieta's formulas (x1 + x2 = -b/a,
// x1 x2 = c/a) and the residual of the equation at each root. An equation typed as text is held
// against the same equation built from its coefficients. Each check gets a case it must fail: a
// cubic and an equation in two unknowns are refused by name.

#include <catch2/catch_test_macros.hpp>
#include <rsc/spatium/closed_form.hpp>
#include <spatium/algebra/polynomial.hpp>
#include <spatium/io/steps_text.hpp>

#include <bit>
#include <cmath>
#include <cstdint>
#include <random>
#include <string>

namespace {

bool same(double x, double y) { return std::bit_cast<std::uint64_t>(x) == std::bit_cast<std::uint64_t>(y); }

} // namespace

TEST_CASE("solution steps: the roots the steps end in are the solver's, and Vieta's", "[algebra][steps]") {
    using spatium::QuadraticCase;
    std::mt19937 rng(11);
    std::uniform_int_distribution<int> k(-6, 6);
    int seen[6] = {};
    for (int i = 0; i < 20000; ++i) {
        const double a = k(rng), b = k(rng), c = k(rng);
        const auto s = spatium::solve_quadratic_steps(a, b, c);
        const auto r = spatium::solve_quadratic(a, b, c);
        REQUIRE(s.roots.size() == r.size());
        for (std::size_t j = 0; j < r.size(); ++j)
            CHECK((same(s.roots[j].re, r[j].re) && same(s.roots[j].im, r[j].im)));
        ++seen[static_cast<int>(s.kind)];
        INFO(a << " x^2 + " << b << " x + " << c);
        switch (s.kind) {
            case QuadraticCase::Identity: CHECK((a == 0 && b == 0 && c == 0)); break;
            case QuadraticCase::None: CHECK((a == 0 && b == 0 && c != 0 && s.roots.size() == 0)); break;
            case QuadraticCase::Linear: CHECK(std::abs(b * s.roots[0].re + c) <= 1e-12 * (std::abs(c) + 1)); break;
            default: {
                CHECK(s.discriminant == b * b - 4 * a * c);
                CHECK((s.kind == QuadraticCase::ComplexPair) == (s.discriminant < 0));
                CHECK((s.kind == QuadraticCase::Double) == (s.discriminant == 0));
                const auto& x1 = s.roots[0];
                const auto& x2 = s.roots[1];
                // Vieta, in complex arithmetic: the sum and the product of the roots
                CHECK(std::abs(x1.re + x2.re + b / a) <= 1e-12 * (std::abs(b / a) + 1));
                CHECK(std::abs(x1.im + x2.im) <= 1e-12);
                CHECK(std::abs(x1.re * x2.re - x1.im * x2.im - c / a) <= 1e-12 * (std::abs(c / a) + 1));
                CHECK(std::abs(x1.re * x2.im + x1.im * x2.re) <= 1e-12 * (std::abs(c / a) + 1));
            }
        }
    }
    for (int kind = 0; kind < 6; ++kind) CHECK(seen[kind] > 0);   // every case was reached
}

TEST_CASE("solution steps: Vieta's formulas check the roots, and the stable form keeps the small one", "[algebra][steps]") {
    std::mt19937 rng(5);
    std::uniform_int_distribution<int> k(-9, 9);
    for (int i = 0; i < 20000; ++i) {
        const double a = k(rng), b = k(rng), c = k(rng);
        if (a == 0) continue;
        const auto v = spatium::check_vieta(spatium::solve_quadratic_steps(a, b, c));
        REQUIRE(v.has_value());
        INFO(a << " x^2 + " << b << " x + " << c);
        CHECK(v->error <= 1e-12);
    }
    // x^2 + 1e8 x + 1, whose small root the textbook form lost to cancellation (-7.45e-9 against
    // -1e-8, the product of the roots 0.745): the solver's stable form keeps it
    const auto steep = spatium::solve_quadratic_steps(1.0, 1e8, 1.0);
    const auto kept = spatium::check_vieta(steep);
    REQUIRE(kept.has_value());
    CHECK(kept->error <= 1e-12);
    CHECK(std::abs(steep.roots[0].re + 1e-8) <= 1e-22);
    // the case it must fail: roots that are not the equation's, one off by the textbook's loss
    auto off = steep;
    off.roots = {spatium::Complex<double>{-7.450580596923828e-09}, steep.roots[1]};
    const auto lost = spatium::check_vieta(off);
    REQUIRE(lost.has_value());
    CHECK(lost->error > 0.1);
    const auto step = spatium::io::build::to_step(*lost);
    CHECK(step.label == "check: Vieta's formulas do not hold");
    CHECK(step.latex.find("\\ne") != std::string::npos);
    // a linear equation has one root and nothing to check this way
    CHECK_FALSE(spatium::check_vieta(spatium::solve_quadratic_steps(0.0, 2.0, 1.0)).has_value());
}

TEST_CASE("solution steps: a quadratic as a chain a page reads", "[io][steps]") {
    namespace bd = spatium::io::build;
    const auto two = bd::to_steps(spatium::solve_quadratic_steps(2.0, -3.0, -5.0));
    REQUIRE(two.size() == 5);
    CHECK(bd::to_step(*spatium::check_vieta(spatium::solve_quadratic_steps(2.0, -3.0, -5.0))).latex ==
          "x_{1} + x_{2} = 1.5 = -\\frac{b}{a} = 1.5,\\quad x_{1} x_{2} = -2.5 = \\frac{c}{a} = -2.5");
    CHECK(two[0].latex == "2x^{2} - 3x - 5 = 0");
    CHECK(two[1].latex == "a = 2,\\quad b = -3,\\quad c = -5");
    CHECK(two[2].latex == "D = b^{2} - 4ac = \\left(-3\\right)^{2} - 4 \\cdot 2 \\cdot \\left(-5\\right) = 49");
    CHECK(two[3].label == "D > 0: two real roots");
    CHECK(two[3].latex == "x_{1,2} = \\frac{-b \\pm \\sqrt{D}}{2a} = \\frac{3 \\pm \\sqrt{49}}{4}");
    CHECK(two[4].latex == "x_{1} = 2.5,\\quad x_{2} = -1");

    const auto pair = bd::to_steps(spatium::solve_quadratic_steps(1.0, 2.0, 5.0), "t");
    CHECK(pair[0].latex == "t^{2} + 2t + 5 = 0");
    CHECK(pair[3].label == "D < 0: two complex conjugate roots");
    CHECK(pair.back().latex == "t_{1} = -1 + 2\\,i,\\quad t_{2} = -1 - 2\\,i");

    const auto twice = bd::to_steps(spatium::solve_quadratic_steps(1.0, -2.0, 1.0));
    CHECK(twice[3].label == "D = 0: one root, counted twice");
    CHECK(twice.back().latex == "x = 1");
    const auto line = bd::to_steps(spatium::solve_quadratic_steps(0.0, 2.0, -1.0));
    CHECK(line.back().latex == "x = -\\frac{c}{b} = -\\frac{-1}{2} = 0.5");
}

TEST_CASE("solution steps: an equation typed as text is solved in closed form, or refused by name", "[rsc][steps]") {
    const auto typed = rsc::solve_closed_form("2x^2 - 3x - 5 = 0");
    REQUIRE(typed.has_value());
    const auto built = spatium::solve_quadratic_steps(2.0, -3.0, -5.0);
    CHECK((typed->solution.a == 2 && typed->solution.b == -3 && typed->solution.c == -5));
    CHECK(same(typed->solution.roots[0].re, built.roots[0].re));
    CHECK(same(typed->solution.roots[1].re, built.roots[1].re));
    CHECK(typed->unknown == "x");

    // an equation not written in the form: its coefficients come from its values, not from rewriting
    const auto product = rsc::solve_closed_form("(y - 1) * (y + 2) = 4 - y");
    REQUIRE(product.has_value());
    CHECK(product->unknown == "y");
    CHECK((product->solution.a == 1 && product->solution.b == 2 && product->solution.c == -6));

    // the cases it must fail
    const auto cubic = rsc::solve_closed_form("x^3 - 1 = 0");
    REQUIRE_FALSE(cubic.has_value());
    CHECK(cubic.error().message.find("not a polynomial of degree at most two in x") != std::string::npos);
    const auto two_names = rsc::solve_closed_form("x + y = 0");
    REQUIRE_FALSE(two_names.has_value());
    CHECK(two_names.error().message.find("exactly one unknown") != std::string::npos);
    CHECK_FALSE(rsc::solve_closed_form("x^2 - 1").has_value());
}
