#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/spaces/hyperbolic.hpp>
#include <spatium/point.hpp>

using namespace spatium;
using Catch::Matchers::WithinAbs;

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
