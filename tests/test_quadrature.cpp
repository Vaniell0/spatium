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

// Gauss-Kronrod through the one door.
template<class F, class T>
auto gk(F&& f, T a, T b) {
    return integrate(f, Finite<T>{a, b}, IntegralOptions{.method = Method::GaussKronrod});
}
}

TEST_CASE("tanh-sinh on a finite interval agrees with the closed form", "[quadrature]") {
    auto r = integrate([](double x) { return x * x; }, Finite<double>{0.0, 1.0});
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(1.0 / 3.0, 1e-12));

    r = integrate([](double x) { return std::sin(x); }, Finite<double>{0.0, kPi});
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(2.0, 1e-12));

    // Reversed ends: the sign of the oriented integral.
    r = integrate([](double x) { return x * x; }, Finite<double>{1.0, 0.0});
    CHECK_THAT(r.value, WithinAbs(-1.0 / 3.0, 1e-12));
    CHECK(integrate([](double) { return 1.0; }, Finite<double>{2.0, 2.0}).value == 0.0);
}

TEST_CASE("tanh-sinh takes an integrable singularity at an end, which Simpson answered NaN", "[quadrature]") {
    // 1/sqrt(x) on [0,1] = 2, log x = -1: the end is never sampled.
    auto inv_sqrt = [](double x) { return 1.0 / std::sqrt(x); };
    auto r = integrate(inv_sqrt, Finite<double>{0.0, 1.0});
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(2.0, 1e-9));
    CHECK(!std::isfinite(integrate_with_error(inv_sqrt, 0.0, 1.0).value));

    r = integrate([](double x) { return std::log(x); }, Finite<double>{0.0, 1.0});
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(-1.0, 1e-10));

    // Both ends at once: 1/sqrt(x(1-x)) integrates to pi. (The near-1 end is
    // an ulp of 1 away at best: the honest limit stated in the header.)
    r = integrate([](double x) { return 1.0 / std::sqrt(x * (1.0 - x)); }, Finite<double>{0.0, 1.0});
    CHECK_THAT(r.value, WithinAbs(kPi, 1e-6));
}

TEST_CASE("exp-sinh covers a half line, either way", "[quadrature]") {
    auto r = integrate([](double x) { return std::exp(-x); }, HalfLine<double>{0.0});
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(1.0, 1e-11));

    r = integrate([](double x) { return 1.0 / (1.0 + x * x); }, HalfLine<double>{0.0});
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(kPi / 2, 1e-10));

    // Toward minus infinity, from a point that is not 0.
    r = integrate([](double x) { return std::exp(x - 2.0); }, HalfLine<double>{2.0, false});
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(1.0, 1e-11));

    // Singular at the finite end AND infinite: x^(-1/2) e^-x = sqrt(pi).
    r = integrate([](double x) { return std::exp(-x) / std::sqrt(x); }, HalfLine<double>{0.0});
    CHECK_THAT(r.value, WithinAbs(std::sqrt(kPi), 1e-9));
}

TEST_CASE("sinh-sinh covers the whole line", "[quadrature]") {
    auto r = integrate([](double x) { return std::exp(-x * x); }, WholeLine{});
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(std::sqrt(kPi), 1e-11));

    r = integrate([](double x) { return 1.0 / (1.0 + x * x); }, WholeLine{});
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(kPi, 1e-10));
}

TEST_CASE("the error estimate covers the actual error on the smooth cases", "[quadrature]") {
    const auto r = integrate([](double x) { return std::exp(x); }, Finite<double>{0.0, 1.0});
    const double actual = std::abs(r.value - (std::exp(1.0) - 1.0));
    CHECK(actual <= 10.0 * r.error_estimate + 1e-15);
    CHECK(r.evaluations > 0);
}

TEST_CASE("quadrature does not claim an answer it does not have", "[quadrature]") {
    // 1/x on (0,1] diverges: no rule may call that converged.
    const auto d = integrate([](double x) { return 1.0 / x; }, Finite<double>{0.0, 1.0});
    CHECK(!d.trusted());

    // A bound that is NaN or infinite, and a function that is nowhere finite.
    CHECK(integrate([](double x) { return x; }, Finite<double>{0.0, kNaN}).status == IntegralStatus::Failed);
    CHECK(integrate([](double x) { return x; }, HalfLine<double>{kNaN}).status == IntegralStatus::Failed);
    const auto nowhere = integrate([](double) { return kNaN; }, WholeLine{});
    CHECK(nowhere.status == IntegralStatus::Failed);
    CHECK(std::isnan(nowhere.value));

    // exp(x) over the whole line diverges at +inf.
    CHECK(!integrate([](double x) { return std::exp(x); }, WholeLine{}).trusted());
}

TEST_CASE("a value that is not finite where the integrand still matters fails the integral", "[quadrature]") {
    // NaN on half of [0, 1]: this used to be read as the end of that side and answered
    // 0.5008 with DepthCapped. It is not an end; the integral has no value.
    const auto half = integrate([](double x) { return x < 0.5 ? 1.0 : kNaN; }, Finite<double>{0.0, 1.0});
    CHECK(half.status == IntegralStatus::Failed);
    CHECK(std::isnan(half.value));

    // sqrt of a negative number is the same thing from the other side of the middle.
    const auto root = integrate([](double x) { return std::sqrt(x); }, Finite<double>{-1.0, 1.0});
    CHECK(root.status == IntegralStatus::Failed);

    // And on the half line, where the nodes far out are the ones that count.
    const auto tail = integrate([](double x) { return x < 3.0 ? std::exp(-x) : kNaN; }, HalfLine<double>{0.0});
    CHECK(tail.status == IntegralStatus::Failed);
}

TEST_CASE("a value that overflows where the integrand is already nothing does not fail it", "[quadrature]") {
    // x^2 exp(-x): far out, x^2 overflows and exp(-x) underflows, so the product is inf * 0 = NaN
    // at nodes that weigh nothing. The integral is 2.
    auto r = integrate([](double x) { return x * x * std::exp(-x); }, HalfLine<double>{0.0});
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(2.0, 1e-9));

    // sin(x)/x over [0, 1]: no node of the rule is at 0, so it is an ordinary integral.
    r = integrate([](double x) { return std::sin(x) / x; }, Finite<double>{0.0, 1.0});
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(0.946083070367183, 1e-10));
}

TEST_CASE("a middle node that is not finite fails the integral instead of leaving its weight out", "[quadrature]") {
    // Over [-1, 1] the middle node is x = 0, where sin(x)/x is 0/0. Dropping that node dropped its
    // weight from every level's sum: 1.8891 for 1.8922, labelled untrusted. It is Failed now; the
    // removable point is the caller's to remove (shift the interval, or return the limit).
    const auto r = integrate([](double x) { return std::sin(x) / x; }, Finite<double>{-1.0, 1.0});
    CHECK(r.status == IntegralStatus::Failed);
    const auto fixed = integrate([](double x) { return x == 0.0 ? 1.0 : std::sin(x) / x; }, Finite<double>{-1.0, 1.0});
    CHECK(fixed.trusted());
    CHECK_THAT(fixed.value, WithinAbs(1.8921661407343662, 1e-10));
}

TEST_CASE("the tolerance follows the scalar: float and double both converge", "[quadrature]") {
    const auto f = integrate([](float x) { return std::exp(-x); }, HalfLine<float>{0.0f});
    CHECK(f.trusted());
    CHECK_THAT(static_cast<double>(f.value), WithinAbs(1.0, 1e-4));
}

TEST_CASE("derivative of an integral through the domain, on Dual", "[quadrature][dual]") {
    using D = Dual<double>;
    // I(p) = int_0^inf exp(-p x) dx = 1/p, I'(p) = -1/p^2.
    const D p = D::variable(2.0);
    const auto r = integrate([p](D x) { return exp(D(-1.0) * p * x); }, HalfLine<D>{D(0.0)});
    CHECK_THAT(r.value.value, WithinAbs(0.5, 1e-10));
    CHECK_THAT(r.value.deriv, WithinAbs(-0.25, 1e-9));

    // And the finite case: d/dp int_0^1 x^p dx = d/dp 1/(p+1) = -1/(p+1)^2.
    const D q = D::variable(1.5);
    const auto s = integrate([q](D x) { return exp(q * log(x)); }, Finite<D>{D(0.0), D(1.0)});
    CHECK_THAT(s.value.value, WithinAbs(1.0 / 2.5, 1e-9));
    CHECK_THAT(s.value.deriv, WithinAbs(-1.0 / (2.5 * 2.5), 1e-8));
}

#if SPATIUM_HAS_BOOST_MULTIPRECISION
TEST_CASE("fifty digits on Real50 where the node table of a fixed rule would end", "[quadrature][precision]") {
    // The nodes are formulas in the scalar: nothing to copy to fifty digits.
    const auto r = integrate([](Real50 x) { return log(x); }, Finite<Real50>{Real50(0), Real50(1)});
    CHECK(r.trusted());
    CHECK(abs(r.value + Real50(1)) < Real50(1e-30));

    const auto h = integrate([](Real50 x) { return exp(-x); }, HalfLine<Real50>{Real50(0)});
    CHECK(h.trusted());
    CHECK(abs(h.value - Real50(1)) < Real50(1e-30));
}
#endif

// ── the second witness ─────────────────────────────────────────

TEST_CASE("Gauss-Kronrod is exact to its degree and agrees with the closed form", "[quadrature]") {
    // K15 integrates polynomials to degree 22 exactly: x^14 over [-1,1] = 2/15.
    auto r = gk([](double x) { return std::pow(x, 14); }, -1.0, 1.0);
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(2.0 / 15.0, 1e-14));

    r = gk([](double x) { return std::sin(x); }, 0.0, kPi);
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(2.0, 1e-12));

    // A kink: adaptive subdivision, not a lucky global rule.
    r = gk([](double x) { return std::abs(x - 0.3); }, 0.0, 1.0);
    CHECK_THAT(r.value, WithinAbs(0.5 * 0.3 * 0.3 + 0.5 * 0.7 * 0.7, 1e-11));
    CHECK(r.evaluations > 15);
}

TEST_CASE("two witnesses that share no node agree within their estimates", "[quadrature][symmetry]") {
    const auto f = [](double x) { return std::exp(-x * x) * std::cos(3.0 * x); };
    const auto a = integrate(f, Finite<double>{-2.0, 3.0});
    const auto b = gk(f, -2.0, 3.0);
    const auto c = integrate_with_error(f, -2.0, 3.0);
    REQUIRE(a.trusted());
    REQUIRE(b.trusted());
    CHECK(std::abs(a.value - b.value) <= a.error_estimate + b.error_estimate + 1e-12);
    CHECK_THAT(c.value, WithinAbs(a.value, 1e-8));
}

TEST_CASE("Gauss-Kronrod refuses a bad bound, a non-finite value, and a scalar it cannot serve", "[quadrature]") {
    CHECK(gk([](double x) { return x; }, 0.0, kNaN).status == IntegralStatus::Failed);
    CHECK(gk([](double) { return kNaN; }, 0.0, 1.0).status == IntegralStatus::Failed);
    CHECK(gk([](double x) { return x; }, 1.0, 1.0).value == 0.0);
    // Real50 is refused at compile time (static_assert), not answered to 19 digits.
}

TEST_CASE("Gauss-Kronrod carries a derivative", "[quadrature][dual]") {
    using D = Dual<double>;
    const D p = D::variable(1.5);
    const auto r = gk([p](D x) { return exp(p * log(x)); }, D(0.0), D(1.0));
    CHECK_THAT(r.value.deriv, WithinAbs(-1.0 / (2.5 * 2.5), 1e-8));
}

// ── the ends decide whether it exists ──────────────────────────

TEST_CASE("the p-test reads the order of f at the ends, exactly", "[quadrature][series]") {
    // 1/sqrt(x) on [0,1]: order -1/2 at 0, integrable; finite at 1.
    auto inv_root = [](auto x) { return 1.0 / sqrt(x); };
    auto e = improper_ends(inv_root, Finite<double>{0.0, 1.0});
    CHECK(e.integrable());
    CHECK(e.first.order == -0.5);

    // 1/x on [0,1]: order -1 is the boundary, and it diverges (logarithmically).
    auto inv = [](auto x) { return 1.0 / x; };
    e = improper_ends(inv, Finite<double>{0.0, 1.0});
    CHECK(e.first.verdict == EndBehaviour::NonIntegrable);
    CHECK(e.second.verdict == EndBehaviour::Integrable);

    // 1/x^2 on [1, inf): decays like x^-2, integrable; 1/x on [1, inf): x^-1, not.
    e = improper_ends([](auto x) { return 1.0 / (x * x); }, HalfLine<double>{1.0});
    CHECK(e.integrable());
    CHECK(e.second.order == 2.0);
    CHECK(improper_ends(inv, HalfLine<double>{1.0}).second.verdict == EndBehaviour::NonIntegrable);

    // x^-1.5 integrates at infinity (ramification 2), 1/sqrt(x) does not.
    CHECK(improper_ends([](auto x) { return pow(x, -1.5); }, HalfLine<double>{1.0}).second.verdict == EndBehaviour::Integrable);
    CHECK(improper_ends(inv_root, HalfLine<double>{1.0}).second.verdict == EndBehaviour::NonIntegrable);

    // a growing function and a constant do not decay at all
    CHECK(improper_ends([](auto x) { return x * 1.0; }, HalfLine<double>{0.0}).non_integrable());
    CHECK(improper_ends([](auto x) { return x * 0.0 + 1.0; }, HalfLine<double>{0.0}).non_integrable());

    // the whole line: 1/(1+x^2) integrable at both infinities, x/(1+x^2) not
    CHECK(improper_ends([](auto x) { return 1.0 / (x * x + 1.0); }, WholeLine{}).integrable());
    CHECK(improper_ends([](auto x) { return x / (x * x + 1.0); }, WholeLine{}).non_integrable());
}

TEST_CASE("a checked domain names a divergence instead of running out of levels", "[quadrature][series]") {
    auto inv = [](auto x) { return 1.0 / x; };
    const auto d = integrate(inv, checked(Finite<double>{0.0, 1.0}));
    CHECK(d.status == IntegralStatus::Divergent);
    CHECK(d.evaluations == 0);
    CHECK(std::isnan(d.value));
    CHECK(!d.trusted());

    CHECK(integrate(inv, checked(HalfLine<double>{1.0})).status == IntegralStatus::Divergent);
    CHECK(integrate([](auto x) { return x * 1.0; }, checked(WholeLine{})).status == IntegralStatus::Divergent);
}

TEST_CASE("where the series cannot say, the rule goes on alone", "[quadrature][series]") {
    // log x at 0 is not a Laurent series: undetermined, and the numerical answer stands.
    auto lg = [](auto x) { return log(x); };
    CHECK(improper_ends(lg, Finite<double>{0.0, 1.0}).first.verdict == EndBehaviour::Undetermined);
    const auto r = integrate(lg, checked(Finite<double>{0.0, 1.0}));
    CHECK(r.trusted());
    CHECK_THAT(r.value, WithinAbs(-1.0, 1e-10));

    // exp(-x) decays faster than any power: not a series at infinity either.
    auto ex = [](auto x) { return exp(x * -1.0); };
    CHECK(improper_ends(ex, HalfLine<double>{0.0}).second.verdict == EndBehaviour::Undetermined);
    CHECK_THAT(integrate(ex, checked(HalfLine<double>{0.0})).value, WithinAbs(1.0, 1e-11));

    // conditionally convergent: sin(x)/x is not a series at infinity; no verdict either way.
    CHECK(improper_ends([](auto x) { return sin(x) / x; }, HalfLine<double>{1.0}).second.verdict == EndBehaviour::Undetermined);

    // and where the ends are fine the number is the one `quadrature` gives
    const auto g = integrate([](auto x) { return 1.0 / (x * x + 1.0); }, checked(WholeLine{}));
    CHECK(g.trusted());
    CHECK_THAT(g.value, WithinAbs(kPi, 1e-10));
}

// ── the door's options ─────────────────────────────────────────

TEST_CASE("one door, and each method is the same answer through it", "[quadrature][door]") {
    const auto f = [](double x) { return std::exp(-x * x) * std::cos(3.0 * x); };
    const Finite<double> d{-2.0, 3.0};
    const auto a = integrate(f, d);
    const auto b = integrate(f, d, {.method = Method::GaussKronrod});
    const auto c = integrate(f, d, {.method = Method::Simpson});
    CHECK_THAT(a.value, WithinAbs(b.value, 1e-9));
    CHECK_THAT(a.value, WithinAbs(c.value, 1e-8));
    // the Simpson method is the rule integrate(f, a, b) of calculus.hpp, unchanged
    CHECK(c.value == integrate_with_error<double>(f, -2.0, 3.0).value);
    CHECK(integrate<double>(f, -2.0, 3.0) == c.value);
    // an infinite domain has one rule; asking for a finite-interval rule is a refusal, not a guess
    CHECK(integrate([](double x) { return std::exp(-x); }, HalfLine<double>{0.0}, {.method = Method::Simpson}).status
          == IntegralStatus::Failed);
}

TEST_CASE("a witness sees what one rule's estimate cannot", "[quadrature][door][symmetry]") {
    // x^0.7 by Simpson at a loose tolerance: Converged, with a true error several times its
    // estimate (test_calculus.cpp). A second rule that samples differently disagrees.
    const auto f = [](double x) { return std::pow(x, 0.7); };
    const Finite<double> d{0.0, 1.0};
    const auto alone = integrate(f, d, {.tolerance = 1e-3, .method = Method::Simpson});
    CHECK(alone.status == IntegralStatus::Converged);
    const auto witnessed = integrate(f, d, {.tolerance = 1e-3, .method = Method::Simpson, .witness = true});
    CHECK(witnessed.status == IntegralStatus::Suspicious);
    CHECK(witnessed.evaluations > alone.evaluations);
    CHECK(witnessed.error_estimate > alone.error_estimate);

    // and where the primary is right, the witness agrees and leaves its status alone
    const auto smooth = integrate([](double x) { return std::sin(x); }, Finite<double>{0.0, kPi}, {.witness = true});
    CHECK(smooth.trusted());
    CHECK_THAT(smooth.value, WithinAbs(2.0, 1e-12));
    // Gauss-Kronrod as the primary is witnessed by the doubly exponential rule
    CHECK(integrate([](double x) { return std::sin(x); }, Finite<double>{0.0, kPi},
                    {.method = Method::GaussKronrod, .witness = true}).trusted());
}

#if SPATIUM_HAS_BOOST_MULTIPRECISION
TEST_CASE("the domain carries Real50 like any other scalar", "[quadrature][door][precision]") {
    const auto w = integrate([](Real50 x) { return exp(-x * x); }, WholeLine<Real50>{});
    CHECK(w.trusted());
    CHECK(abs(w.value - sqrt(Real50(std::numbers::pi_v<long double>))) < Real50(1e-17));   // the constant above is long double's
}
#endif
