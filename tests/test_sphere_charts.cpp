// The two stereographic charts of S^N: embedding <-> chart maps for both, the
// inversion between them, velocities through the conformal factor, geodesics
// that are the same curve in the embedding and in either chart, and a flight
// that crosses the equator from the north chart into the south one.
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <spatium/physics/mechanics/manifold_body.hpp>
#include <spatium/spaces/sphere_charts.hpp>
#include <cmath>
#include <random>

using namespace spatium;
using namespace spatium::physics::mechanics;
using Catch::Approx;

namespace {

template<std::size_t N>
Vec<double, N + 1> random_point(std::mt19937& g, double R, int hemisphere) {   // +1 north, -1 south
    std::normal_distribution<double> n;
    Vec<double, N + 1> y;
    for (std::size_t i = 0; i <= N; ++i) y[i] = n(g);
    y = Vec<double, N + 1>{y * (R / y.norm())};
    if (y[N] * hemisphere < 0) y[N] = -y[N];
    return y;
}

template<std::size_t N>
Vec<double, N + 1> random_tangent(std::mt19937& g, const Vec<double, N + 1>& y, double speed) {
    std::normal_distribution<double> n;
    Vec<double, N + 1> w;
    for (std::size_t i = 0; i <= N; ++i) w[i] = n(g);
    w = Vec<double, N + 1>{w - y * (y.dot(w) / y.dot(y))};
    return Vec<double, N + 1>{w * (speed / w.norm())};
}

template<std::size_t N>
void check_charts(double kappa) {
    std::mt19937 g(31);
    SphereCharts<N> C{kappa};
    const double R = C.radius();
    for (int hemi : {1, -1}) {
        for (int trial = 0; trial < 50; ++trial) {
            const auto y = random_point<N>(g, R, hemi);
            if (std::abs(y[N]) > 0.97 * R && false) continue;
            const auto xn = C.to_north(y), xs = C.to_south(y);
            // round trip, both charts
            REQUIRE(Vec<double, N + 1>{C.from_north(xn) - y}.norm() < 1e-12 * R);
            REQUIRE(Vec<double, N + 1>{C.from_south(xs) - y}.norm() < 1e-12 * R);
            // the transition is the inversion and an involution
            REQUIRE(Vec<double, N>{C.transition(xn) - xs}.norm() < 1e-12 * R * (1 + xs.norm() / R));
            REQUIRE(Vec<double, N>{C.transition(C.transition(xn)) - xn}.norm() < 1e-12 * R * (1 + xn.norm() / R));
            // the chart of the point's own hemisphere is inside the ball |x| < R
            const auto own = hemi == 1 ? xn : xs;
            REQUIRE(own.norm() <= R * (1 + 1e-14));
            // the chart's distance is the sphere's: from the pole of the chart to the point
            const auto sp = C.sphere();
            const auto pole = Vec<double, N + 1>{[&] { Vec<double, N + 1> p{}; p[N] = hemi * R; return p; }()};
            REQUIRE(C.space().distance(Vec<double, N>{}, own) == Approx(sp.distance(pole, y)).epsilon(1e-12).margin(1e-12));
            // velocity: push then pull is the identity, and the length is the conformal factor's
            const auto w = random_tangent<N>(g, y, 1.3);
            const auto un = C.push_north(y, w), us = C.push_south(y, w);
            REQUIRE(Vec<double, N + 1>{C.pull_north(xn, un) - w}.norm() < 1e-10);
            REQUIRE(Vec<double, N + 1>{C.pull_south(xs, us) - w}.norm() < 1e-10);
            REQUIRE(C.space().conformal_factor(xn) * un.norm() == Approx(w.norm()).epsilon(1e-10));
            REQUIRE(C.space().conformal_factor(xs) * us.norm() == Approx(w.norm()).epsilon(1e-10));
            // the transition's derivative carries one chart vector to the other
            REQUIRE(Vec<double, N>{C.transition_push(xn, un) - us}.norm() < 1e-9 * (1 + us.norm()));
        }
    }
}

} // namespace

TEST_CASE("sphere charts: round trips, inversion, velocities", "[charts][sphere]") {
    check_charts<2>(1.0);
    check_charts<3>(1.0);
    check_charts<3>(2.25);
    check_charts<2>(0.3);
}

TEST_CASE("sphere charts: a geodesic is the geodesic of either chart", "[charts][sphere]") {
    std::mt19937 g(8);
    for (double kappa : {1.0, 0.5}) {
        SphereCharts<3> C{kappa};
        const double R = C.radius();
        const auto S = C.sphere();
        const auto M = C.space();
        // start in the north hemisphere near the pole, fly a quarter circle (well inside the chart)
        for (int trial = 0; trial < 20; ++trial) {
            auto y = random_point<3>(g, R, +1);
            const auto w = random_tangent<3>(g, y, 0.7);
            const auto xn = C.to_north(y);
            const auto un = C.push_north(y, w);
            for (double t : {0.1, 0.4, 0.9}) {
                auto body = PointOnManifold<Sphere<3>>{1.0, y, w};
                geodesic_step(body, S, t);
                if (body.position[3] <= 0) continue;                 // left the hemisphere
                const auto in_chart = M.exp_map(xn, un, t);
                REQUIRE(Vec<double, 3>{C.to_north(body.position) - in_chart}.norm() < 1e-10 * R);
                // and the velocity at the end is the pushed forward chart velocity (derivative of exp)
                const double h = 1e-6;
                const auto fwd = M.exp_map(xn, un, t + h), bwd = M.exp_map(xn, un, t - h);
                const Vec<double, 3> du{(fwd - bwd) / (2 * h)};
                REQUIRE(Vec<double, 3>{C.push_north(body.position, body.velocity) - du}.norm() < 1e-6);
            }
        }
    }
}

TEST_CASE("sphere charts: a flight across the equator reappears at the same coordinates", "[charts][sphere]") {
    // A body starts in the north hemisphere and crosses y_N = 0 at the point p. In the north chart
    // it reaches |x| = R at p; in the south chart it enters at the SAME coordinates p, and from
    // there the south chart's own geodesic continues the flight.
    for (double kappa : {1.0, 4.0}) {
        SphereCharts<3> C{kappa};
        const double R = C.radius();
        const auto S = C.sphere();
        const auto M = C.space();
        Vec<double, 4> y{0.3 * R, 0.0, 0.1 * R, 0.0};
        y[3] = std::sqrt(R * R - y[0] * y[0] - y[2] * y[2]);
        const Vec<double, 4> w{0.2, 0.5, -0.6, -0.9};
        auto v = Vec<double, 4>{w - y * (y.dot(w) / y.dot(y))};
        PointOnManifold<Sphere<3>> body{1.0, y, v};
        // the crossing time: y_N(t) = y_N cos(a) + v_N R sin(a)/|v| = 0, a = |v| t / R
        const double speed = v.norm();
        const double a = std::atan2(-y[3] * speed, v[3] * R);
        const double tc = (a > 0 ? a : a + std::numbers::pi) * R / speed;
        auto cross = body;
        geodesic_step(cross, S, tc);
        REQUIRE(std::abs(cross.position[3]) < 1e-12 * R);
        const auto pn = C.to_north(cross.position), ps = C.to_south(cross.position);
        REQUIRE(Vec<double, 3>{pn - ps}.norm() < 1e-12 * R);                 // the same coordinates
        REQUIRE(pn.norm() == Approx(R).epsilon(1e-12));                       // on the rim of both discs
        // Just before: in the north chart's ball; just after: in the south chart's ball.
        auto before = body, after = body;
        geodesic_step(before, S, tc * 0.999);
        geodesic_step(after, S, tc * 1.001);
        REQUIRE(before.position[3] > 0);
        REQUIRE(after.position[3] < 0);
        REQUIRE(C.to_north(before.position).norm() < R);
        REQUIRE(C.to_south(after.position).norm() < R);
        // Continue the flight inside the south chart alone, from the crossing point and velocity.
        const auto us = C.push_south(cross.position, cross.velocity);
        const double dt = 0.2 * R;
        auto later = cross;
        geodesic_step(later, S, dt);
        const auto in_chart = M.exp_map(ps, us, dt);
        REQUIRE(Vec<double, 3>{C.to_south(later.position) - in_chart}.norm() < 1e-10 * R);
    }
}
