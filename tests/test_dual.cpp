#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/dual.hpp>
#include <spatium/algebra/vector.hpp>

using namespace spatium;
using Catch::Matchers::WithinAbs;

TEST_CASE("Dual default is zero", "[dual]") {
    Dual64 d;
    CHECK(d.value == 0.0);
    CHECK(d.deriv == 0.0);
}

TEST_CASE("Dual variable seeds derivative to 1", "[dual]") {
    auto x = Dual64::variable(3.0);
    CHECK(x.value == 3.0);
    CHECK(x.deriv == 1.0);
}

TEST_CASE("Dual arithmetic follows product/quotient rule", "[dual]") {
    auto x = Dual64::variable(2.0);
    auto y = x * x; // d/dx x^2 at x=2 -> 4
    CHECK(y.value == 4.0);
    CHECK(y.deriv == 4.0);

    auto z = x / Dual64::constant(2.0); // d/dx (x/2) -> 0.5, everywhere
    CHECK(z.value == 1.0);
    CHECK(z.deriv == 0.5);
}

TEST_CASE("derivative() recovers known closed-form derivatives", "[dual]") {
    // f(x) = 3x + 2 -> f'(x) = 3, everywhere
    auto linear = [](Dual64 x) { return Dual64::constant(3.0) * x + Dual64::constant(2.0); };
    CHECK_THAT(derivative(linear, 5.0), WithinAbs(3.0, 1e-12));

    // f(x) = x^3 -> f'(x) = 3x^2
    auto cube = [](Dual64 x) { return x * x * x; };
    CHECK_THAT(derivative(cube, 2.0), WithinAbs(12.0, 1e-12));

    // f(x) = sqrt(x) -> f'(x) = 1/(2*sqrt(x))
    auto root = [](Dual64 x) { return sqrt(x); };
    CHECK_THAT(derivative(root, 4.0), WithinAbs(0.25, 1e-12));

    // f(x) = sin(x) -> f'(x) = cos(x)
    auto sine = [](Dual64 x) { return sin(x); };
    CHECK_THAT(derivative(sine, 0.0), WithinAbs(1.0, 1e-12));

    // f(x) = exp(x) -> f'(x) = exp(x)
    auto expf = [](Dual64 x) { return exp(x); };
    CHECK_THAT(derivative(expf, 1.0), WithinAbs(std::exp(1.0), 1e-12));

    // chain rule: f(x) = sin(x^2) -> f'(x) = 2x*cos(x^2)
    auto composed = [](Dual64 x) { return sin(x * x); };
    CHECK_THAT(derivative(composed, 1.5), WithinAbs(2.0 * 1.5 * std::cos(1.5 * 1.5), 1e-9));
}

TEST_CASE("Dual supports scalar-on-the-left arithmetic", "[dual]") {
    auto x = Dual64::variable(2.0);
    auto a = 10.0 * x;      // d/dx (10x) = 10
    auto b = 5.0 - x;       // d/dx (5-x) = -1
    auto c = 3.0 + x;       // d/dx (3+x) = 1
    CHECK(a.value == 20.0); CHECK(a.deriv == 10.0);
    CHECK(b.value == 3.0);  CHECK(b.deriv == -1.0);
    CHECK(c.value == 5.0);  CHECK(c.deriv == 1.0);
}

TEST_CASE("Dual satisfies Scalar and drops into existing generic code unchanged", "[dual]") {
    // Vec<T,N>::norm() is written as `using std::sqrt; return sqrt(...)`
    // (ADL-friendly, per convention) — substituting Dual64 for double here
    // differentiates the Euclidean norm with zero changes to vector.hpp.
    static_assert(Scalar<Dual64>);

    auto x = Dual64::variable(3.0); // differentiate w.r.t. this component
    Vec<Dual64, 3> v{x, Dual64::constant(4.0), Dual64::constant(0.0)};
    auto n = v.norm(); // sqrt(3^2+4^2) = 5

    CHECK_THAT(n.value, WithinAbs(5.0, 1e-12));
    // d/dx sqrt(x^2+16) at x=3 -> x/sqrt(x^2+16) = 3/5
    CHECK_THAT(n.deriv, WithinAbs(0.6, 1e-12));
}
