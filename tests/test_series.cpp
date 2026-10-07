#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/series.hpp>
#include <spatium/algebra/dual.hpp>
#include <spatium/core/precision.hpp>
#include <cmath>
#include <limits>
#include <numbers>

using namespace spatium;
using Catch::Matchers::WithinAbs;

namespace {
using S = Series<double, 8>;
const double kInf = std::numeric_limits<double>::infinity();
const S eps = S::monomial(1.0, 1);
}

static_assert(Scalar<Series<double, 8>>);
static_assert(Scalar<Series<float, 4>>);

TEST_CASE("series arithmetic is exact on what it keeps", "[series]") {
    // (1 + eps)(1 - eps) = 1 - eps^2
    const S p = (S(1.0) + eps) * (S(1.0) - eps);
    CHECK(p.valuation() == 0);
    CHECK(p.coefficient(0) == 1.0);
    CHECK(p.coefficient(1) == 0.0);
    CHECK(p.coefficient(2) == -1.0);

    // 1 / (1 - eps) = 1 + eps + eps^2 + ...
    const S g = S(1.0) / (S(1.0) - eps);
    for (int k = 0; k < 8; ++k) CHECK(g.coefficient(k) == 1.0);

    // a pole keeps its order: 1/eps
    const S pole = S(1.0) / eps;
    CHECK(pole.valuation() == -1);
    CHECK(pole.coefficient(-1) == 1.0);

    // a pole minus itself is the exact zero, not O(anything)
    CHECK((pole - pole).is_exact_zero());
}

TEST_CASE("a cancellation spends precision and says so", "[series]") {
    // sin(eps)/eps^3 - 1/eps^2 = -1/6 + ..., known to fewer terms than sin(eps) itself.
    const S s = sin(eps);
    const S d = s / (eps * eps * eps) - S(1.0) / (eps * eps);
    CHECK_THAT(d.coefficient(0), WithinAbs(-1.0 / 6.0, 1e-15));
    CHECK(d.valuation() == 0);
    CHECK(d.known_to() < S::kExact);
    CHECK(d.known_to() > 3);
}

TEST_CASE("limits of the forms Simpson, extrapolation and a division of zeros cannot take", "[series][limit]") {
    const auto sinc = limit([](auto x) { return sin(x) / x; }, 0.0);
    CHECK(sinc.kind == LimitKind::Finite);
    CHECK_THAT(sinc.value, WithinAbs(1.0, 1e-15));

    const auto one_minus_cos = limit([](auto x) { return (S(1.0) - cos(x)) / (x * x); }, 0.0);
    CHECK_THAT(one_minus_cos.value, WithinAbs(0.5, 1e-15));

    const auto expm1 = limit([](auto x) { return (exp(x) - S(1.0)) / x; }, 0.0);
    CHECK_THAT(expm1.value, WithinAbs(1.0, 1e-15));

    // Three terms cancelled: (sin x - x) / x^3 -> -1/6.
    const auto cubic = limit([](auto x) { return (sin(x) - x) / (x * x * x); }, 0.0);
    CHECK_THAT(cubic.value, WithinAbs(-1.0 / 6.0, 1e-14));

    // 1^infinity shaped: (1 + x)^(1/x) -> e, through exp(log(1+x)/x).
    const auto e = limit([](auto x) { return exp(log(S(1.0) + x) / x); }, 0.0);
    CHECK_THAT(e.value, WithinAbs(std::numbers::e, 1e-14));

    // At a point that is not 0.
    const auto at_two = limit([](auto x) { return (x * x - S(4.0)) / (x - S(2.0)); }, 2.0);
    CHECK_THAT(at_two.value, WithinAbs(4.0, 1e-14));
}

TEST_CASE("a pole comes back as a pole, with its order and its sign", "[series][limit]") {
    const auto inv = limit([](auto x) { return S(1.0) / x; }, 0.0);
    CHECK(inv.kind == LimitKind::Infinite);
    CHECK(inv.value == kInf);
    CHECK(inv.order() == -1.0);

    // the two sides disagree on 1/x and agree on 1/x^2
    CHECK(limit([](auto x) { return S(1.0) / x; }, 0.0, Approach::Left).value == -kInf);
    CHECK(!two_sided_limit([](auto x) { return S(1.0) / x; }, 0.0).determined());
    const auto sq = two_sided_limit([](auto x) { return S(1.0) / (x * x); }, 0.0);
    CHECK(sq.kind == LimitKind::Infinite);
    CHECK(sq.value == kInf);
    CHECK(sq.order() == -2.0);

    // cot(x) at 0+ is +inf of order -1
    const auto cot = limit([](auto x) { return cos(x) / sin(x); }, 0.0);
    CHECK(cot.kind == LimitKind::Infinite);
    CHECK(cot.order() == -1.0);

    // two poles that cancel: 1/x - 1/sin(x) -> 0, and of order 1.
    const auto d = limit([](auto x) { return S(1.0) / x - S(1.0) / sin(x); }, 0.0);
    CHECK(d.kind == LimitKind::Finite);
    CHECK_THAT(d.value, WithinAbs(0.0, 1e-15));
    CHECK(d.exponent == 1);
}

TEST_CASE("fractional powers by ramification", "[series][limit]") {
    const auto root = limit([](auto x) { return sqrt(x); }, 0.0);
    CHECK(root.kind == LimitKind::Finite);
    CHECK(root.value == 0.0);
    CHECK(root.ramification == 2);
    CHECK(root.order() == 0.5);

    // 1/sqrt(x): the integrable singularity of the quadrature tests, now with its order
    const auto inv_root = limit([](auto x) { return S(1.0) / sqrt(x); }, 0.0);
    CHECK(inv_root.kind == LimitKind::Infinite);
    CHECK(inv_root.order() == -0.5);

    const auto two_thirds = limit([](auto x) { return pow(x, 2.0 / 3.0); }, 0.0);
    CHECK(two_thirds.ramification == 3);
    CHECK(two_thirds.order() > 0.66);
    CHECK(two_thirds.order() < 0.67);
}

TEST_CASE("the point at infinity is a chart", "[series][limit]") {
    const auto ratio = limit_at_infinity([](auto x) { return (x * x + S(1.0)) / (x * x - S(1.0)); });
    CHECK(ratio.kind == LimitKind::Finite);
    CHECK_THAT(ratio.value, WithinAbs(1.0, 1e-15));

    // x sin(1/x) -> 1
    const auto xs = limit_at_infinity([](auto x) { return x * sin(S(1.0) / x); });
    CHECK_THAT(xs.value, WithinAbs(1.0, 1e-15));

    // (1 + 1/x)^x -> e
    const auto euler = limit_at_infinity([](auto x) { return exp(x * log(S(1.0) + S(1.0) / x)); });
    CHECK_THAT(euler.value, WithinAbs(std::numbers::e, 1e-14));

    // x/(x^2 + 1) decays like 1/x: the order that decides an improper integral
    const auto decay = limit_at_infinity([](auto x) { return x / (x * x + S(1.0)); });
    CHECK(decay.kind == LimitKind::Finite);
    CHECK(decay.value == 0.0);
    CHECK(decay.order() == 1.0);

    // sqrt(x^2 + 1) - x -> 0 like 1/(2x): a difference of two infinities
    const auto diff = limit_at_infinity([](auto x) { return sqrt(x * x + S(1.0)) - x; });
    CHECK(diff.value == 0.0);
    CHECK(diff.order() == 1.0);
    CHECK_THAT(diff.coefficient, WithinAbs(0.5, 1e-15));

    // x grows: a pole of the chart
    CHECK(limit_at_infinity([](auto x) { return x * x; }).kind == LimitKind::Infinite);
    CHECK(limit_at_infinity([](auto x) { return x * x * x; }, true).value == -kInf);
}

TEST_CASE("what is not a Laurent series is refused, not guessed", "[series][limit]") {
    CHECK(!limit([](auto x) { return exp(S(1.0) / x); }, 0.0).determined());     // essential singularity
    CHECK(!limit([](auto x) { return log(x); }, 0.0).determined());              // log(eps) is not a series
    CHECK(!limit([](auto x) { return sin(S(1.0) / x); }, 0.0).determined());     // oscillation
    CHECK(!limit_at_infinity([](auto x) { return sin(x); }).determined());
    // division by an exact zero
    CHECK(!limit([](auto x) { return S(1.0) / (x - x); }, 0.0).determined());
}

TEST_CASE("the series limit and the number at a small step agree", "[series][limit][symmetry]") {
    // Two paths to one answer that share nothing: the expansion, and the function
    // evaluated near the point. Where f is smooth the gap is the step's own error.
    const auto check = [](auto f_series, auto f_double, double a) {
        const auto r = limit(f_series, a);
        REQUIRE(r.kind == LimitKind::Finite);
        CHECK_THAT(f_double(a + 1e-5), WithinAbs(r.value, 1e-4));
    };
    check([](auto x) { return sin(x) / x; }, [](double x) { return std::sin(x) / x; }, 0.0);
    check([](auto x) { return (exp(x) - S(1.0)) / x; }, [](double x) { return std::expm1(x) / x; }, 0.0);
    check([](auto x) { return tan(x) / x; }, [](double x) { return std::tan(x) / x; }, 0.0);
    check([](auto x) { return (cosh(x) - S(1.0)) / (x * x); }, [](double x) { return (std::cosh(x) - 1.0) / (x * x); }, 0.0);
}

TEST_CASE("the order is that of eps to 0+, and equality is of the known terms", "[series]") {
    CHECK(eps > S(0.0));
    CHECK(eps < S(1.0));
    CHECK(S(1.0) / eps > S(1e9));       // a pole is larger than any number
    CHECK(-eps < S(0.0));
    CHECK(S(2.0) == S(2.0));
    CHECK(abs(-eps) == eps);
    CHECK(!(S::variable(0.0, 1) < S::variable(0.0, 1)));
}

TEST_CASE("the derivative of a limit, with Dual as the coefficient", "[series][dual]") {
    using D = Dual<double>;
    using SD = Series<D, 8>;
    // lim_{x->0} sin(p x)/x = p, so d/dp = 1.
    const D p = D::variable(2.0);
    const auto r = limit<8>([p](auto x) { return sin(SD(p) * x) / x; }, D(0.0));
    REQUIRE(r.kind == LimitKind::Finite);
    CHECK_THAT(r.value.value, WithinAbs(2.0, 1e-14));
    CHECK_THAT(r.value.deriv, WithinAbs(1.0, 1e-14));
}

#if SPATIUM_HAS_BOOST_MULTIPRECISION
TEST_CASE("fifty digits of a limit on Real50", "[series][precision]") {
    using R = Real50;
    using SR = Series<R, 8>;
    const auto r = limit<8>([](auto x) { return sin(x) / x; }, R(0));
    REQUIRE(r.kind == LimitKind::Finite);
    CHECK(abs(r.value - R(1)) < R(1e-45));
    const auto c = limit<8>([](auto x) { return (sin(x) - x) / (x * x * x); }, R(0));
    CHECK(abs(c.value + R(1) / R(6)) < R(1e-45));
    (void)sizeof(SR);
}
#endif
