// An ImplicitSurface whose function is a type (spaces/implicit.hpp,
// make_implicit): exact gradient and Hessian through nested Duals, the
// Gaussian curvature that follows from the two, and geodesics in the ambient
// space, x'' = -(x'^T H x') / |grad F|^2 grad F.
//
// Three independent routes meet here and must agree: the unit sphere as a
// level set against the closed forms of Sphere<2>; a torus as a level set
// (the ambient geodesic equation) against the same torus as a chart (the
// induced metric's Christoffel symbols, spaces/parametric.hpp); and the
// Gaussian curvature of the level set against the closed form -- and, with
// the chart's area element, against Gauss-Bonnet.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/calculus.hpp>
#include <spatium/core/access.hpp>
#include <spatium/spaces/implicit.hpp>
#include <spatium/spaces/parametric.hpp>
#include <spatium/spaces/sphere.hpp>
#include <cmath>
#include <numbers>

using namespace spatium;
using Catch::Matchers::WithinAbs;

namespace {

constexpr double R = 2.0, q = 0.7;
constexpr double two_pi = 2 * std::numbers::pi;

auto implicit_sphere() {
    return make_implicit<double>(
        [](auto x, auto y, auto z) { using std::sqrt; return sqrt(x * x + y * y + z * z) - 1.0; },
        {-2.0, 2.0, -2.0, 2.0, -2.0, 2.0});
}

auto implicit_torus() {
    return make_implicit<double>(
        [](auto x, auto y, auto z) {
            using std::sqrt;
            const auto d = sqrt(x * x + y * y) - R;
            return d * d + z * z - q * q;
        },
        {-4.0, 4.0, -4.0, 4.0, -2.0, 2.0});
}

auto chart_torus() {
    return make_parametric<double>(
        [](auto u, auto v) {
            using std::cos; using std::sin;
            return Vec<decltype(u), 3>{(R + q * cos(v)) * cos(u), (R + q * cos(v)) * sin(u), q * sin(v)};
        },
        {0.0, two_pi, 0.0, two_pi}, true, true);
}

}  // namespace

TEST_CASE("A generic function makes the gradient and the Hessian exact; a std::function keeps a difference quotient",
          "[implicit][typed]") {
    const auto typed = make_implicit<double>(
        [](auto x, auto y, auto z) { return x * x + y * y + z * z - 1.0; }, {-2.0, 2.0, -2.0, 2.0, -2.0, 2.0});
    const ImplicitSurface<double> erased(
        [](double x, double y, double z) { return x * x + y * y + z * z - 1.0; }, {-2.0, 2.0, -2.0, 2.0, -2.0, 2.0});
    static_assert(decltype(typed)::exact);
    static_assert(!ImplicitSurface<double>::exact);

    const Vec<double, 3> p{0.6, 0.0, 0.8};
    const auto g = typed.gradient(p);
    CHECK((g - Vec<double, 3>{1.2, 0.0, 1.6}).norm() < 1e-14);                 // 2p, exactly
    const auto H = typed.hessian(p);
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j) CHECK_THAT(H(i, j), WithinAbs(i == j ? 2.0 : 0.0, 1e-14));
    // The erased gradient has a difference quotient's accuracy at a step of
    // 128 eps * 1000, far from exact.
    CHECK((erased.gradient(p) - g).norm() > 1e-9);
}

TEST_CASE("The Gaussian curvature of a level set is exact: the unit sphere is 1, the torus is cos v / (q (R + q cos v))",
          "[implicit][typed]") {
    CHECK_THAT(implicit_sphere().gaussian_curvature(Vec<double, 3>{0.6, 0.0, 0.8}), WithinAbs(1.0, 1e-13));

    const auto implicit = implicit_torus();
    const auto chart = chart_torus();
    for (const double u : {0.3, 1.7, 4.0})
        for (const double v : {0.2, 1.1, 2.9, 4.5}) {
            const double closed = std::cos(v) / (q * (R + q * std::cos(v)));
            CHECK_THAT(implicit.gaussian_curvature(chart.evaluate(u, v)), WithinAbs(closed, 1e-12));
        }
}

TEST_CASE("Gauss-Bonnet through two descriptions: the level set's curvature, the chart's area", "[implicit][typed]") {
    // The integral of K dA is 2 pi chi: 4 pi on the sphere, 0 on the torus.
    // K from the implicit function, dA from the parametrization -- neither
    // knows the answer -- by the library's own quadrature, whose status is read.
    const auto sphere_chart = make_parametric<double>(
        [](auto u, auto v) {
            using std::cos; using std::sin;
            return Vec<decltype(u), 3>{sin(v) * cos(u), sin(v) * sin(u), cos(v)};
        },
        {0.0, two_pi, 0.0, std::numbers::pi}, true, false);
    const auto sphere = implicit_sphere();
    // K dA does not depend on u on a surface of revolution: the u integral is 2 pi.
    const auto over_sphere = integrate_with_error<double>(
        [&](double v) { return sphere.gaussian_curvature(sphere_chart.evaluate(0.7, v)) * sphere_chart.area_element(0.7, v); },
        1e-6, std::numbers::pi - 1e-6);
    CHECK(over_sphere.trusted());
    CHECK_THAT(two_pi * over_sphere.value, WithinAbs(4 * std::numbers::pi, 1e-4));

    const auto torus = implicit_torus();
    const auto torus_chart = chart_torus();
    const auto over_torus = integrate_with_error<double>(
        [&](double v) { return torus.gaussian_curvature(torus_chart.evaluate(0.7, v)) * torus_chart.area_element(0.7, v); },
        0.0, two_pi);
    CHECK(over_torus.trusted());
    CHECK_THAT(two_pi * over_torus.value, WithinAbs(0.0, 1e-9));
}

TEST_CASE("The unit sphere as a level set has the geodesics of Sphere<2>", "[implicit][typed]") {
    const auto level = implicit_sphere();
    const Sphere<2> closed;
    const Vec<double, 3> p{0.6, 0.0, 0.8};
    for (const Vec<double, 3>& w : {Vec<double, 3>{0.0, 1.0, 0.0}, Vec<double, 3>{0.8, 0.0, -0.6},
                                     Vec<double, 3>{0.5, 0.4, -0.375}}) {
        // exp: a point reached by the flow against the great circle
        const auto a = level.exp_map(p, w, 1.3);
        const auto b = closed.exp_map(p, closed.project_tangent(p, w), 1.3);
        CHECK((a - b).norm() < 1e-8);
        // log and distance
        const auto target = closed.exp_map(p, closed.project_tangent(p, w), 0.9);
        CHECK_THAT(spaces::distance(level, p, target), WithinAbs(closed.distance(p, target), 1e-8));
    }
    // exp after log is the identity
    const Vec<double, 3> a{0.6, 0.0, 0.8}, b{0.0, 0.8, 0.6};
    CHECK((level.exp_map(a, level.log_map(a, b), 1.0) - b).norm() < 1e-8);
}

TEST_CASE("The torus as a level set and as a chart have the same geodesics", "[implicit][typed]") {
    // The ambient geodesic equation and the induced metric's Christoffel
    // symbols share nothing but the surface.
    const auto level = implicit_torus();
    const auto chart = chart_torus();
    const double u0 = 0.4, v0 = 0.9;
    const auto p = chart.evaluate(u0, v0);
    const Vec<double, 3> fu{-(R + q * std::cos(v0)) * std::sin(u0), (R + q * std::cos(v0)) * std::cos(u0), 0};
    const Vec<double, 3> fv{-q * std::sin(v0) * std::cos(u0), -q * std::sin(v0) * std::sin(u0), q * std::cos(v0)};
    const Vec<double, 3> w{fu * (0.6 / (R + q * std::cos(v0))) + fv * (0.8 / q)};
    for (const double t : {0.5, 1.0, 2.5}) {
        CHECK((level.exp_map(p, w, t) - chart.exp_map(p, w, t)).norm() < 1e-8);
    }
    const auto b = chart.evaluate(0.9, 1.4);
    CHECK_THAT(spaces::distance(level, p, b), WithinAbs(spaces::distance(chart, p, b), 1e-8));
    CHECK((level.log_map(p, b) - chart.log_map(p, b)).norm() < 1e-7);
}

TEST_CASE("A typed implicit surface is erased at the boundary, and NaN stays NaN", "[implicit][typed][infinity]") {
    const auto level = implicit_sphere();
    const ImplicitSurface<double> erased = level.erased();
    CHECK_THAT(erased(0.6, 0.0, 0.8), WithinAbs(0.0, 1e-15));
    CHECK_THAT(erased.distance(Vec<double, 3>{1, 0, 0}, Vec<double, 3>{0, 1, 0}), WithinAbs(std::sqrt(2.0), 1e-15));   // the chord
    CHECK(spaces::distance(level, Vec<double, 3>{1, 0, 0}, Vec<double, 3>{0, 1, 0}) > std::sqrt(2.0));          // |log|
    const double nan = std::numeric_limits<double>::quiet_NaN();
    CHECK_FALSE(std::isfinite(level.exp_map(Vec<double, 3>{nan, 0, 1}, Vec<double, 3>{1, 0, 0}, 1.0)[0]));
}
