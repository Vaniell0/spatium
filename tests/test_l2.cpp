// L^2([a, b]) (spaces/l2.hpp): the first space with no finite basis.
//
// Parseval's identity gives the norm from the coefficients with no quadrature, and a quadrature
// of f^2 gives it independently: those two are held to each other. A series of an analytic
// function reproduces it; one of a function it does not suit (a kink) says so in its tail.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/core/verify.hpp>
#include <spatium/spaces/l2.hpp>
#include <cmath>
#include <vector>

using namespace spatium;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("a series reproduces an analytic function", "[l2]") {
    const L2Interval<double> space{-1.0, 2.0};
    const auto f = [](double x) { return std::exp(std::sin(x)) + 0.3 * x; };
    const auto s = space.from_function(f, 24);
    for (const double x : {-1.0, -0.4, 0.0, 0.9, 1.7, 2.0})
        CHECK_THAT(space.evaluate(s, x), WithinAbs(f(x), 1e-9));
    CHECK(s.tail < 1e-9);
}

TEST_CASE("Parseval against quadrature: the norm from the coefficients is the norm of the function", "[l2][symmetry]") {
    const L2Interval<double> space{-1.0, 2.0};
    struct Case { const char* name; double (*f)(double); };
    const Case cases[] = {
        {"exp", [](double x) { return std::exp(x); }},
        {"cos 3x", [](double x) { return std::cos(3 * x); }},
        {"1/(1+x^2)", [](double x) { return 1.0 / (1.0 + x * x); }},
    };
    for (const auto& c : cases) {
        const auto s = space.from_function(c.f, 40);
        const double by_coefficients = space.inner(s, s);
        const double by_quadrature = integrate([&](double x) { return c.f(x) * c.f(x); }, Finite<double>{-1.0, 2.0}).value;
        INFO(c.name);
        CHECK_THAT(by_coefficients, WithinRel(by_quadrature, 1e-9));
    }
}

TEST_CASE("the Legendre polynomials are an orthogonal basis of known norms", "[l2]") {
    const L2Interval<double> space{0.0, 3.0};
    for (std::size_t j = 0; j < 5; ++j)
        for (std::size_t k = 0; k < 5; ++k) {
            const double expected = j == k ? 3.0 / (2.0 * j + 1.0) : 0.0;     // (b - a)/(2 n + 1)
            CHECK_THAT(space.inner(space.basis(j), space.basis(k)), WithinAbs(expected, 1e-13));
        }
}

TEST_CASE("it is a Hilbert space: the metric axioms and Pythagoras", "[l2]") {
    const L2Interval<double> space{-1.0, 1.0};
    const std::vector<LegendreSeries<double>> pts{
        space.from_function([](double x) { return std::exp(x); }, 20),
        space.from_function([](double x) { return std::cos(2 * x); }, 20),
        space.from_function([](double x) { return x * x * x; }, 20),
        space.from_function([](double x) { return 1.0 / (2.0 + x); }, 20)};
    const auto r = verify_metric(space, pts, 1e-9);
    INFO(r.message());
    CHECK(r.passed);
    // orthogonal vectors: |u + v|^2 = |u|^2 + |v|^2
    const auto u = space.basis(2), v = space.basis(5);
    CHECK_THAT(space.inner(u + v, u + v), WithinAbs(space.inner(u, u) + space.inner(v, v), 1e-13));
    // Cauchy-Schwarz
    CHECK(std::abs(space.inner(pts[0], pts[1])) <= space.norm(pts[0]) * space.norm(pts[1]) + 1e-12);
}

TEST_CASE("a function the series does not suit says so in its tail", "[l2]") {
    const L2Interval<double> space{-1.0, 1.0};
    const auto smooth = space.from_function([](double x) { return std::exp(x); }, 24);
    const auto kink = space.from_function([](double x) { return std::abs(x); }, 24);
    CHECK(smooth.tail < 1e-10);
    CHECK(kink.tail > 1e-4);                      // a kink: the coefficients decay like 1/n^2, not geometrically
    CHECK(kink.tail > 1e4 * smooth.tail);
}
