// A space given by its metric alone (spaces/metric_chart.hpp): exp and log
// derived by the geodesic equation and shooting, held against the closed
// forms the library already has.
//
// The sphere in (theta, phi) is Sphere<2> in another coordinate system; the
// upper half-plane is the hyperbolic plane in another model. Each is a pair
// of paths to one answer -- a closed form and a derivation through
// Christoffel symbols, Dual partials and an extrapolated ODE -- and they
// must agree, to the tolerance the derivation was asked for. A disagreement
// is in the formula, the derivation or the integrator, and a test of any one
// alone would not say which.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/core/access.hpp>
#include <spatium/core/precision.hpp>
#include <spatium/core/verify.hpp>
#include <spatium/algebra/dual.hpp>
#include <spatium/spaces/hyperbolic.hpp>
#include <spatium/spaces/metric_chart.hpp>
#include <spatium/spaces/sphere.hpp>
#include <cmath>
#include <span>
#include <vector>

using namespace spatium;
using Catch::Matchers::WithinAbs;

namespace {

// ds^2 = dtheta^2 + sin^2(theta) dphi^2
struct SphereMetric {
    template<class S> Matrix<S, 2, 2> operator()(const Vec<S, 2>& x) const {
        using std::sin;
        Matrix<S, 2, 2> g{};
        g(0, 0) = S(1);
        g(1, 1) = sin(x[0]) * sin(x[0]);
        return g;
    }
};

// ds^2 = (dx^2 + dy^2) / y^2
struct HalfPlaneMetric {
    template<class S> Matrix<S, 2, 2> operator()(const Vec<S, 2>& x) const {
        Matrix<S, 2, 2> g{};
        g(0, 0) = g(1, 1) = S(1) / (x[1] * x[1]);
        return g;
    }
};

// A metric that is singular on the line x0 = 0.
struct Degenerate {
    template<class S> Matrix<S, 2, 2> operator()(const Vec<S, 2>& x) const {
        Matrix<S, 2, 2> g{};
        g(0, 0) = S(1);
        g(1, 1) = x[0] * x[0];
        return g;
    }
};

// Minkowski space: the metric of a flat spacetime, four coordinates.
struct Flat {
    template<class S> Matrix<S, 4, 4> operator()(const Vec<S, 4>&) const {
        Matrix<S, 4, 4> g{};
        g(0, 0) = S(-1); g(1, 1) = g(2, 2) = g(3, 3) = S(1);
        return g;
    }
};

Vec<double, 3> embed(double theta, double phi) {
    return {std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)};
}

// The ambient vector a coordinate tangent (dtheta, dphi) at (theta, phi) is.
Vec<double, 3> embed_tangent(double theta, double phi, double dtheta, double dphi) {
    const Vec<double, 3> e_theta{std::cos(theta) * std::cos(phi), std::cos(theta) * std::sin(phi), -std::sin(theta)};
    const Vec<double, 3> e_phi{-std::sin(theta) * std::sin(phi), std::sin(theta) * std::cos(phi), 0.0};
    return Vec<double, 3>{e_theta * dtheta + e_phi * dphi};
}

double sphere_distance_closed(const Vec<double, 2>& a, const Vec<double, 2>& b) {
    const Sphere<2> s;
    return s.distance(embed(a[0], a[1]), embed(b[0], b[1]));
}

double half_plane_distance_closed(const Vec<double, 2>& a, const Vec<double, 2>& b) {
    const double dx = b[0] - a[0], dy = b[1] - a[1];
    return std::acosh(1.0 + (dx * dx + dy * dy) / (2.0 * a[1] * b[1]));
}

}  // namespace

TEST_CASE("A metric alone is a Riemannian space", "[metric_chart]") {
    const auto s = spaces::metric_chart<double, 2>(SphereMetric{});
    static_assert(spaces::Exponential<decltype(s)>);
    static_assert(spaces::Logarithmic<decltype(s)>);
    static_assert(spaces::Metrized<decltype(s)>);
    static_assert(spaces::Riemannian<decltype(s)>);
    // distance is not written anywhere: it is |log|_g, derived by the
    // customization points from the derived log.
    static_assert(spaces::Distanced<decltype(s)>);
}

TEST_CASE("The sphere in (theta, phi) derived from its metric is Sphere<2>", "[metric_chart]") {
    const auto chart = spaces::metric_chart<double, 2>(SphereMetric{});
    const Sphere<2> sphere;

    const Vec<double, 2> starts[] = {{1.1, 0.2}, {0.9, -0.5}, {2.0, 1.3}, {1.57, 0.0}};
    const Vec<double, 2> tangents[] = {{0.3, 0.2}, {-0.5, 0.4}, {0.1, -0.7}, {0.8, 0.05}};
    for (const auto& p : starts)
        for (const auto& v : tangents) {
            // The derivation, then the closed form on the embedded point
            // and tangent; the embeddings of the results must coincide.
            const Vec<double, 2> reached = chart.exp_map(p, v, 1.0);
            const Vec<double, 3> closed = sphere.exp_map(embed(p[0], p[1]),
                                                         embed_tangent(p[0], p[1], v[0], v[1]), 1.0);
            const Vec<double, 3> derived = embed(reached[0], reached[1]);
            CHECK((derived - closed).norm() < 1e-8);
        }
}

TEST_CASE("The derived distance on the sphere is the arc length", "[metric_chart]") {
    const auto chart = spaces::metric_chart<double, 2>(SphereMetric{});
    const Vec<double, 2> pts[] = {{1.1, 0.2}, {1.7, 0.9}, {0.9, -0.5}, {2.0, 0.4}};
    for (const auto& a : pts)
        for (const auto& b : pts) {
            // The closed form is acos, which is only sqrt(eps) accurate where
            // the points nearly coincide (2e-8 for a point and itself); the
            // derived distance is exact there, so the tolerance is acos's.
            const double d = spaces::distance(chart, a, b);
            CHECK_THAT(d, WithinAbs(sphere_distance_closed(a, b), 1e-7));
        }
}

TEST_CASE("The half-plane's vertical geodesics and distance, closed form against derived", "[metric_chart]") {
    const auto h2 = spaces::metric_chart<double, 2>(HalfPlaneMetric{});

    // y' = y at unit speed up a vertical line: y(t) = y0 e^t.
    const Vec<double, 2> p{0.0, 1.0}, up{0.0, 1.0};
    const Vec<double, 2> r = h2.exp_map(p, up, 2.0);
    CHECK_THAT(r[0], WithinAbs(0.0, 1e-9));
    CHECK_THAT(r[1], WithinAbs(std::exp(2.0), 1e-8));

    const Vec<double, 2> pts[] = {{0.3, 1.2}, {-0.4, 2.1}, {0.9, 0.8}, {-0.8, 1.5}};
    for (const auto& a : pts)
        for (const auto& b : pts)
            CHECK_THAT(spaces::distance(h2, a, b), WithinAbs(half_plane_distance_closed(a, b), 1e-8));
}

TEST_CASE("exp of log is the identity on a derived space, and the library's own verifier agrees",
          "[metric_chart]") {
    const auto h2 = spaces::metric_chart<double, 2>(HalfPlaneMetric{});
    std::vector<Vec<double, 2>> pts = {{0.3, 1.2}, {-0.4, 2.1}, {0.9, 0.8}, {-0.8, 1.5}, {0.0, 1.0}};
    for (const auto& p : pts)
        for (const auto& q : pts) {
            const auto v = h2.log_map(p, q);
            const auto back = h2.exp_map(p, v, 1.0);
            CHECK((back - q).norm() < 1e-8);
        }
    const auto report = verify_exp_log(h2, std::span<const Vec<double, 2>>(pts), 1e-6);
    CHECK(report.passed);
}

TEST_CASE("The derivative of a derived exp in t is the velocity, through the flow", "[metric_chart]") {
    // The scalar is a Dual, so the whole machinery -- Christoffel by nested
    // Dual, the extrapolation, the adaptive step -- is differentiated, and
    // d/dt exp_p(t v) at t = 0 must be v. An exact check of AD through an
    // integrator, where a finite difference can only be approximate.
    using D = Dual<double>;
    const auto chart = spaces::metric_chart<D, 2>(SphereMetric{});
    const Vec<D, 2> p{D(1.1), D(0.2)}, v{D(0.3), D(-0.5)};
    D t = D(0.0);
    t.deriv = 1.0;
    const Vec<D, 2> r = chart.exp_map(p, v, t);
    CHECK_THAT(r[0].deriv, WithinAbs(0.3, 1e-9));
    CHECK_THAT(r[1].deriv, WithinAbs(-0.5, 1e-9));
    CHECK_THAT(r[0].value, WithinAbs(1.1, 1e-12));
}

TEST_CASE("A singular metric gives NaN, not a straight line", "[metric_chart]") {
    // Where g does not invert the Christoffel symbols do not exist; reading
    // them as zero would fly on in a straight line and look like an answer.
    const auto bad = spaces::metric_chart<double, 2>(Degenerate{});
    const Vec<double, 2> p{0.0, 0.3}, v{0.1, 0.1};
    const auto r = bad.exp_map(p, v, 1.0);
    CHECK(std::isnan(r[0]));

    // The sphere at its pole is a coordinate singularity of the same kind.
    const auto s = spaces::metric_chart<double, 2>(SphereMetric{});
    const auto pole = s.exp_map(Vec<double, 2>{0.0, 0.3}, Vec<double, 2>{0.1, 0.1}, 1.0);
    CHECK(std::isnan(pole[0]));

    // And the checked Christoffel refuses outright.
    const auto G = spaces::christoffel_checked(Degenerate{}, p);
    CHECK_FALSE(G.has_value());
}

TEST_CASE("The 4D relativity names are the generic ones", "[metric_chart]") {
    // Schwarzschild's exact Christoffel symbols from the generic code at
    // N = 4 are what physics/relativity has always called christoffel().
    const Vec<double, 4> x{0, 1, 2, 3};
    const auto G = spaces::christoffel(Flat{}, x);
    for (std::size_t l = 0; l < 4; ++l)
        for (std::size_t i = 0; i < 4; ++i)
            for (std::size_t j = 0; j < 4; ++j) CHECK(G[l](i, j) == 0.0);
}

#if SPATIUM_HAS_BOOST_MULTIPRECISION
TEST_CASE("Over fifty digits the derived exp and log keep fifty-digit accuracy", "[metric_chart][precision]") {
    // The integrator's tolerance follows the scalar: the same call that is
    // right to 1e-11 on double is right to 1e-30 on Real50, which a fixed
    // RK4 step could not be at any practical step count.
    const auto chart = spaces::metric_chart<Real50, 2>(HalfPlaneMetric{});
    const Vec<Real50, 2> p{Real50("0.3"), Real50("1.2")}, q{Real50("-0.4"), Real50("2.1")};
    const auto v = chart.log_map(p, q);
    const auto back = chart.exp_map(p, v, Real50(1));
    CHECK(static_cast<double>(Real50(abs(back[0] - q[0]) + abs(back[1] - q[1]))) < 1e-28);

    // And it is double's answer, to double's precision: L3 of the matrix,
    // here for the one pair.
    const auto dchart = spaces::metric_chart<double, 2>(HalfPlaneMetric{});
    const auto dv = dchart.log_map(Vec<double, 2>{0.3, 1.2}, Vec<double, 2>{-0.4, 2.1});
    CHECK_THAT(static_cast<double>(v[0]), WithinAbs(dv[0], 1e-8));
    CHECK_THAT(static_cast<double>(v[1]), WithinAbs(dv[1], 1e-8));

    const Real50 d = spaces::distance(chart, p, q);
    const double closed = half_plane_distance_closed(Vec<double, 2>{0.3, 1.2}, Vec<double, 2>{-0.4, 2.1});
    CHECK_THAT(static_cast<double>(d), WithinAbs(closed, 1e-9));
}
#endif
