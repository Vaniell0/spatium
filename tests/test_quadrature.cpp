#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/quadrature.hpp>
#include <spatium/algebra/dual.hpp>
#include <spatium/core/precision.hpp>
#include <cmath>
#include <limits>
#include <numbers>

using namespace spatium;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
constexpr double kPi = std::numbers::pi;
const double kNaN = std::numeric_limits<double>::quiet_NaN();
}

TEST_CASE("tanh-sinh on a finite interval agrees with the closed form", "[quadrature]") {
    auto r = quadrature([](double x) { return x * x; }, Finite<double>{0.0, 1.0});
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(1.0 / 3.0, 1e-12));

    r = quadrature([](double x) { return std::sin(x); }, Finite<double>{0.0, kPi});
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(2.0, 1e-12));

    // Reversed ends: the sign of the oriented integral.
    r = quadrature([](double x) { return x * x; }, Finite<double>{1.0, 0.0});
    CHECK_THAT(r.value, WithinAbs(-1.0 / 3.0, 1e-12));
    CHECK(quadrature([](double) { return 1.0; }, Finite<double>{2.0, 2.0}).value == 0.0);
}

TEST_CASE("tanh-sinh takes an integrable singularity at an end, which Simpson answered NaN", "[quadrature]") {
    // 1/sqrt(x) on [0,1] = 2, log x = -1: the end is never sampled.
    auto inv_sqrt = [](double x) { return 1.0 / std::sqrt(x); };
    auto r = quadrature(inv_sqrt, Finite<double>{0.0, 1.0});
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(2.0, 1e-9));
    CHECK(!std::isfinite(integrate_with_error(inv_sqrt, 0.0, 1.0).value));

    r = quadrature([](double x) { return std::log(x); }, Finite<double>{0.0, 1.0});
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(-1.0, 1e-10));

    // Both ends at once: 1/sqrt(x(1-x)) integrates to pi. (The near-1 end is
    // an ulp of 1 away at best: the honest limit stated in the header.)
    r = quadrature([](double x) { return 1.0 / std::sqrt(x * (1.0 - x)); }, Finite<double>{0.0, 1.0});
    CHECK_THAT(r.value, WithinAbs(kPi, 1e-6));
}

TEST_CASE("exp-sinh covers a half line, either way", "[quadrature]") {
    auto r = quadrature([](double x) { return std::exp(-x); }, HalfLine<double>{0.0});
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(1.0, 1e-11));

    r = quadrature([](double x) { return 1.0 / (1.0 + x * x); }, HalfLine<double>{0.0});
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(kPi / 2, 1e-10));

    // Toward minus infinity, from a point that is not 0.
    r = quadrature([](double x) { return std::exp(x - 2.0); }, HalfLine<double>{2.0, false});
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(1.0, 1e-11));

    // Singular at the finite end AND infinite: x^(-1/2) e^-x = sqrt(pi).
    r = quadrature([](double x) { return std::exp(-x) / std::sqrt(x); }, HalfLine<double>{0.0});
    CHECK_THAT(r.value, WithinAbs(std::sqrt(kPi), 1e-9));
}

TEST_CASE("sinh-sinh covers the whole line", "[quadrature]") {
    auto r = quadrature([](double x) { return std::exp(-x * x); }, WholeLine{});
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(std::sqrt(kPi), 1e-11));

    r = quadrature([](double x) { return 1.0 / (1.0 + x * x); }, WholeLine{});
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(kPi, 1e-10));
}

TEST_CASE("the error estimate covers the actual error on the smooth cases", "[quadrature]") {
    const auto r = quadrature([](double x) { return std::exp(x); }, Finite<double>{0.0, 1.0});
    const double actual = std::abs(r.value - (std::exp(1.0) - 1.0));
    CHECK(actual <= 10.0 * r.error_estimate + 1e-15);
    CHECK(r.evaluations > 0);
}

TEST_CASE("quadrature does not claim an answer it does not have", "[quadrature]") {
    // 1/x on (0,1] diverges: no rule may call that converged.
    const auto d = quadrature([](double x) { return 1.0 / x; }, Finite<double>{0.0, 1.0});
    CHECK(!d.trusted());

    // A bound that is NaN or infinite, and a function that is nowhere finite.
    CHECK(quadrature([](double x) { return x; }, Finite<double>{0.0, kNaN}).status == IntegralStatus::Failed);
    CHECK(quadrature([](double x) { return x; }, HalfLine<double>{kNaN}).status == IntegralStatus::Failed);
    const auto nowhere = quadrature([](double) { return kNaN; }, WholeLine{});
    CHECK(nowhere.status == IntegralStatus::Failed);
    CHECK(std::isnan(nowhere.value));

    // exp(x) over the whole line diverges at +inf.
    CHECK(!quadrature([](double x) { return std::exp(x); }, WholeLine{}).trusted());
}

TEST_CASE("the tolerance follows the scalar: float and double both converge", "[quadrature]") {
    const auto f = quadrature([](float x) { return std::exp(-x); }, HalfLine<float>{0.0f});
    CHECK(f.trusted());
    CHECK_THAT(static_cast<double>(f.value), WithinAbs(1.0, 1e-4));
}

TEST_CASE("derivative of an integral through the domain, on Dual", "[quadrature][dual]") {
    using D = Dual<double>;
    // I(p) = int_0^inf exp(-p x) dx = 1/p, I'(p) = -1/p^2.
    const D p = D::variable(2.0);
    const auto r = quadrature([p](D x) { return exp(D(-1.0) * p * x); }, HalfLine<D>{D(0.0)});
    CHECK_THAT(r.value.value, WithinAbs(0.5, 1e-10));
    CHECK_THAT(r.value.deriv, WithinAbs(-0.25, 1e-9));

    // And the finite case: d/dp int_0^1 x^p dx = d/dp 1/(p+1) = -1/(p+1)^2.
    const D q = D::variable(1.5);
    const auto s = quadrature([q](D x) { return exp(q * log(x)); }, Finite<D>{D(0.0), D(1.0)});
    CHECK_THAT(s.value.value, WithinAbs(1.0 / 2.5, 1e-9));
    CHECK_THAT(s.value.deriv, WithinAbs(-1.0 / (2.5 * 2.5), 1e-8));
}

#if SPATIUM_HAS_BOOST_MULTIPRECISION
TEST_CASE("fifty digits on Real50 where the node table of a fixed rule would end", "[quadrature][precision]") {
    // The nodes are formulas in the scalar: nothing to copy to fifty digits.
    const auto r = quadrature([](Real50 x) { return log(x); }, Finite<Real50>{Real50(0), Real50(1)});
    CHECK(r.trusted());
    CHECK(abs(r.value + Real50(1)) < Real50(1e-30));

    const auto h = quadrature([](Real50 x) { return exp(-x); }, HalfLine<Real50>{Real50(0)});
    CHECK(h.trusted());
    CHECK(abs(h.value - Real50(1)) < Real50(1e-30));
}
#endif

// ── the second witness ─────────────────────────────────────────

TEST_CASE("Gauss-Kronrod is exact to its degree and agrees with the closed form", "[quadrature]") {
    // K15 integrates polynomials to degree 22 exactly: x^14 over [-1,1] = 2/15.
    auto r = gauss_kronrod([](double x) { return std::pow(x, 14); }, -1.0, 1.0);
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(2.0 / 15.0, 1e-14));

    r = gauss_kronrod([](double x) { return std::sin(x); }, 0.0, kPi);
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(2.0, 1e-12));

    // A kink: adaptive subdivision, not a lucky global rule.
    r = gauss_kronrod([](double x) { return std::abs(x - 0.3); }, 0.0, 1.0);
    CHECK_THAT(r.value, WithinAbs(0.5 * 0.3 * 0.3 + 0.5 * 0.7 * 0.7, 1e-11));
    CHECK(r.evaluations > 15);
}

TEST_CASE("two witnesses that share no node agree within their estimates", "[quadrature][symmetry]") {
    const auto f = [](double x) { return std::exp(-x * x) * std::cos(3.0 * x); };
    const auto a = quadrature(f, Finite<double>{-2.0, 3.0});
    const auto b = gauss_kronrod(f, -2.0, 3.0);
    const auto c = integrate_with_error(f, -2.0, 3.0);
    REQUIRE(a.trusted());
    REQUIRE(b.trusted());
    CHECK(std::abs(a.value - b.value) <= a.error_estimate + b.error_estimate + 1e-12);
    CHECK_THAT(c.value, WithinAbs(a.value, 1e-8));
}

TEST_CASE("Gauss-Kronrod refuses a bad bound, a non-finite value, and a scalar it cannot serve", "[quadrature]") {
    CHECK(gauss_kronrod([](double x) { return x; }, 0.0, kNaN).status == IntegralStatus::Failed);
    CHECK(gauss_kronrod([](double) { return kNaN; }, 0.0, 1.0).status == IntegralStatus::Failed);
    CHECK(gauss_kronrod([](double x) { return x; }, 1.0, 1.0).value == 0.0);
    // Real50 is refused at compile time (static_assert), not answered to 19 digits.
}

TEST_CASE("Gauss-Kronrod carries a derivative", "[quadrature][dual]") {
    using D = Dual<double>;
    const D p = D::variable(1.5);
    const auto r = gauss_kronrod([p](D x) { return exp(p * log(x)); }, D(0.0), D(1.0));
    CHECK_THAT(r.value.deriv, WithinAbs(-1.0 / (2.5 * 2.5), 1e-8));
}
