#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/spaces/sphere.hpp>
#include <spatium/point.hpp>
#include <numbers>

using namespace spatium;
using Catch::Matchers::WithinAbs;

TEST_CASE("Sphere concept satisfaction", "[sphere]") {
    static_assert(RiemannianManifold<S2>);
    static_assert(Surface<S2>);
    static_assert(Complete<S2>);
    SUCCEED();
}

TEST_CASE("Sphere contains", "[sphere]") {
    S2 sphere;
    CHECK(sphere.contains(Vec3{1.0, 0.0, 0.0}));
    CHECK(sphere.contains(Vec3{0.0, 1.0, 0.0}));
    CHECK_FALSE(sphere.contains(Vec3{1.0, 1.0, 0.0})); // not on unit sphere
}

TEST_CASE("Sphere distance: poles", "[sphere]") {
    S2 sphere;
    Vec3 north{0.0, 0.0, 1.0};
    Vec3 south{0.0, 0.0, -1.0};
    CHECK_THAT(sphere.distance(north, south), WithinAbs(std::numbers::pi, 1e-10));
}

TEST_CASE("Sphere distance: same point", "[sphere]") {
    S2 sphere;
    Vec3 p{1.0, 0.0, 0.0};
    CHECK_THAT(sphere.distance(p, p), WithinAbs(0.0, 1e-12));
}

TEST_CASE("Sphere distance: quarter turn", "[sphere]") {
    S2 sphere;
    Vec3 a{1.0, 0.0, 0.0};
    Vec3 b{0.0, 1.0, 0.0};
    CHECK_THAT(sphere.distance(a, b), WithinAbs(std::numbers::pi / 2.0, 1e-10));
}

TEST_CASE("Sphere exp/log roundtrip", "[sphere]") {
    S2 sphere;
    Vec3 p{1.0, 0.0, 0.0};
    Vec3 q{0.0, 1.0, 0.0};

    auto v = sphere.log_map(p, q);
    auto q2 = sphere.exp_map(p, v, 1.0);

    CHECK_THAT(sphere.distance(q, q2), WithinAbs(0.0, 1e-8));
}

TEST_CASE("Sphere exp_map partial t", "[sphere]") {
    S2 sphere;
    Vec3 p{1.0, 0.0, 0.0};
    Vec3 q{0.0, 1.0, 0.0};

    auto v = sphere.log_map(p, q);
    auto mid = sphere.exp_map(p, v, 0.5);

    // Midpoint should be at pi/4 from both
    CHECK_THAT(sphere.distance(p, mid), WithinAbs(std::numbers::pi / 4.0, 1e-8));
    CHECK_THAT(sphere.distance(mid, q), WithinAbs(std::numbers::pi / 4.0, 1e-8));
}

TEST_CASE("Sphere project", "[sphere]") {
    S2 sphere;
    auto proj = sphere.project(Vec3{3.0, 0.0, 0.0});
    CHECK_THAT(proj[0], WithinAbs(1.0, 1e-12));
    CHECK_THAT(proj[1], WithinAbs(0.0, 1e-12));
    CHECK_THAT(proj[2], WithinAbs(0.0, 1e-12));
}

TEST_CASE("Sphere with radius", "[sphere]") {
    Sphere<2> sphere{.radius = 5.0};
    Vec3 a{5.0, 0.0, 0.0};
    Vec3 b{0.0, 5.0, 0.0};
    CHECK_THAT(sphere.distance(a, b), WithinAbs(5.0 * std::numbers::pi / 2.0, 1e-8));
}

TEST_CASE("Sphere typed Point", "[sphere]") {
    S2 sphere;
    auto p = pt<S2>(Vec3{1.0, 0.0, 0.0});
    auto q = pt<S2>(Vec3{0.0, 1.0, 0.0});
    CHECK_THAT(p.distance_to(q, sphere), WithinAbs(std::numbers::pi / 2.0, 1e-10));
}

// log_map returned a tangent of length theta, the angle, where exp_map
// reads its argument as a length along the sphere and divides by the
// radius: the two were inverse only on the unit sphere, which is the only
// one the tests used. On a sphere of radius 2, exp(log(q)) missed q by
// 0.99; the log's length must be the distance, radius * theta.
TEST_CASE("Sphere exp_map inverts log_map at any radius", "[sphere]") {
    for (double R : {0.5, 1.0, 2.0, 1000.0}) {
        const Sphere<2, double> s{R};
        const Vec<double, 3> p{0.0, 0.0, R};
        for (double angle : {0.3, 1.0, 2.5, std::numbers::pi}) {
            const Vec<double, 3> q{R * std::sin(angle), 0.0, R * std::cos(angle)};
            const auto v = s.log_map(p, q);
            INFO("R " << R << " angle " << angle);
            CHECK_THAT(v.norm(), WithinAbs(s.distance(p, q), 1e-12 * R));
            const double miss = Vec<double, 3>{s.exp_map(p, v, 1.0) - q}.norm();
            CHECK_THAT(miss, WithinAbs(0.0, 1e-9 * R));
        }
    }
}
