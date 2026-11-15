#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/spaces/hyperbolic.hpp>
#include <spatium/point.hpp>
#include <spatium/algebra/calculus.hpp>
#include <spatium/core/access.hpp>
#include <cmath>
#include <random>
#include <vector>

using namespace spatium;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Hyperbolic concept satisfaction", "[hyperbolic]") {
    static_assert(RiemannianManifold<H2>);
    static_assert(Surface<H2>);
    static_assert(Complete<H2>);
    SUCCEED();
}

TEST_CASE("Hyperbolic origin", "[hyperbolic]") {
    H2 space;
    auto o = H2::origin();
    CHECK(o[0] == 1.0);
    CHECK(o[1] == 0.0);
    CHECK(o[2] == 0.0);
    CHECK(space.contains(o));
}

TEST_CASE("Hyperbolic contains", "[hyperbolic]") {
    H2 space;
    // (cosh(1), sinh(1), 0) is on the hyperboloid
    Vec3 p{std::cosh(1.0), std::sinh(1.0), 0.0};
    CHECK(space.contains(p));

    // Not on hyperboloid
    CHECK_FALSE(space.contains(Vec3{1.0, 1.0, 0.0}));
}

TEST_CASE("Hyperbolic distance: same point", "[hyperbolic]") {
    H2 space;
    auto o = H2::origin();
    CHECK_THAT(space.distance(o, o), WithinAbs(0.0, 1e-12));
}

TEST_CASE("Hyperbolic distance: a point is at distance 0 from itself far from the origin", "[hyperbolic][distance_self]") {
    // -<p,p> is 1 plus rounding once p is a few units out, and acosh of that
    // was 1e-8 (d(p,p) = d(exp(p,w), exp(p,w)) failed on H^n).
    H2 space;
    auto p = H2::origin();
    for (double s : {0.5, 3.0, 8.0, 12.0}) {
        p = space.exp_map(p, Vec3{0.0, 0.8, 0.6}, s);
        CHECK(space.distance(p, p) == 0.0);
    }
    Hyperbolic<4> h4;
    auto q4 = h4.exp_map(Hyperbolic<4>::origin(), Vec<double, 5>{0.0, 0.3, -0.5, 0.7, 0.1}, 9.0);
    CHECK(h4.distance(q4, q4) == 0.0);
}

TEST_CASE("Hyperbolic distance: far from the origin agrees with the closed form", "[hyperbolic][distance_self]") {
    // Two points on one geodesic through the origin are |s - t| apart.
    H2 space;
    for (auto [s, t] : {std::pair{7.0, 7.5}, std::pair{10.0, -4.0}, std::pair{1e-4, 2e-4}}) {
        Vec3 a{std::cosh(s), std::sinh(s), 0.0}, b{std::cosh(t), std::sinh(t), 0.0};
        CHECK_THAT(space.distance(a, b), WithinRel(std::abs(s - t), 1e-9));
    }
}

TEST_CASE("Hyperbolic distance: known values", "[hyperbolic]") {
    H2 space;
    auto o = H2::origin(); // (1, 0, 0)
    Vec3 p{std::cosh(2.0), std::sinh(2.0), 0.0}; // distance 2 from origin
    CHECK_THAT(space.distance(o, p), WithinAbs(2.0, 1e-10));
}

TEST_CASE("Hyperbolic exp/log roundtrip", "[hyperbolic]") {
    H2 space;
    auto o = H2::origin();
    Vec3 p{std::cosh(1.5), std::sinh(1.5), 0.0};

    auto v = space.log_map(o, p);
    auto p2 = space.exp_map(o, v, 1.0);

    CHECK_THAT(space.distance(p, p2), WithinAbs(0.0, 1e-6));
}

TEST_CASE("Hyperbolic exp_map partial t", "[hyperbolic]") {
    H2 space;
    auto o = H2::origin();
    Vec3 p{std::cosh(2.0), std::sinh(2.0), 0.0};

    auto v = space.log_map(o, p);
    auto mid = space.exp_map(o, v, 0.5);

    CHECK_THAT(space.distance(o, mid), WithinAbs(1.0, 1e-8));
    CHECK_THAT(space.distance(mid, p), WithinAbs(1.0, 1e-8));
}

TEST_CASE("Hyperbolic project", "[hyperbolic]") {
    H2 space;
    // Project a spatial point (0, 3, 4) onto hyperboloid
    auto proj = space.project(Vec3{0.0, 3.0, 4.0});
    CHECK(space.contains(proj));
    CHECK_THAT(proj[1], WithinAbs(3.0, 1e-12));
    CHECK_THAT(proj[2], WithinAbs(4.0, 1e-12));
    CHECK_THAT(proj[0], WithinAbs(std::sqrt(1.0 + 9.0 + 16.0), 1e-12));
}

TEST_CASE("Hyperbolic typed Point", "[hyperbolic]") {
    H2 space;
    auto p = pt<H2>(H2::origin());
    Vec3 q_raw{std::cosh(1.0), std::sinh(1.0), 0.0};
    auto q = pt<H2>(q_raw);
    CHECK_THAT(p.distance_to(q, space), WithinAbs(1.0, 1e-10));
}

TEST_CASE("Far from the origin exp after log keeps every digit the model can hold", "[hyperbolic]") {
    // Found by enumerating laws over the space's operations: exp_q(log_q p) = p
    // held near the origin and failed once q was 8 away. A tangent at q has
    // coordinates of size d e^d and length d, and both the length exp read as
    // -v0^2 + |v_s|^2 and the norm log read off its projection cancelled: the
    // error grew like e^{3.5 d} -- 2.6e-4 at d = 8, 1.4e3 at d = 12. What the
    // model cannot avoid is p rebuilt from q cosh d, about e^{2d} of rounding.
    const Hyperbolic<3> h;
    const auto p = h.exp_map(Hyperbolic<3>::origin(), Vec<double, 4>{0, 0.3, -0.2, 0.1}, 1.0);
    Vec<double, 4> dir{0, 0.6, 0.8, 0.0};
    dir[0] = (dir[1] * p[1] + dir[2] * p[2] + dir[3] * p[3]) / p[0];   // tangent at p
    dir = dir * (1.0 / std::sqrt(h.metric_at(p, dir, dir)));
    for (const double d : {2.0, 4.0, 8.0, 12.0}) {
        const auto q = h.exp_map(p, dir, d);
        INFO("d = " << d);
        CHECK_THAT(h.distance(p, q), WithinAbs(d, 1e-12 * d));
        // log at the far point has the length of the way, to the last digits
        const auto v = h.log_map(q, p);
        CHECK_THAT(std::sqrt(Hyperbolic<3>::tangent_norm_sq(q, v)), WithinAbs(d, 1e-12 * d));
        const auto back = h.exp_map(q, v, 1.0);
        const double floor = 64 * std::numeric_limits<double>::epsilon() * std::cosh(d) * std::cosh(d);
        for (std::size_t i = 0; i < 4; ++i) CHECK_THAT(back[i], WithinAbs(p[i], floor));
    }
}

TEST_CASE("The Frechet mean converges on widely spread hyperbolic points", "[hyperbolic][frechet]") {
    // Found by the connectivity matrix's dimension axis, at H^8, where random
    // tangents are longer: the Karcher step of 1 overshoots when d coth d
    // exceeds 2 (a spread past about 1.5), so the iteration diverged to NaN.
    // The mean of a point set symmetric about the origin is the origin.
    std::mt19937_64 g(7);
    std::uniform_real_distribution<double> u(-1, 1);
    for (const double spread : {1.0, 2.0, 3.0}) {
        const Hyperbolic<8> h;
        const auto c = Hyperbolic<8>::origin();
        std::vector<Vec<double, 9>> pts;
        for (int k = 0; k < 3; ++k) {
            Vec<double, 9> w;
            w[0] = 0;
            for (std::size_t i = 1; i < 9; ++i) w[i] = u(g) * spread;
            pts.push_back(h.exp_map(c, w, 1.0));
            pts.push_back(h.exp_map(c, w, -1.0));
        }
        const auto m = frechet_mean(h, pts, pts.front(), 1e-12, 500);
        INFO("spread " << spread);
        CHECK(std::isfinite(h.distance(m, c)));
        CHECK(h.distance(m, c) < 1e-6);
    }
}
