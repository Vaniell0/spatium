#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/spaces/euclidean.hpp>

using namespace spatium;
using Catch::Matchers::WithinAbs;

TEST_CASE("Euclidean concept satisfaction", "[euclidean]") {
    // All verified via static_assert in euclidean.hpp
    // This test exists to confirm compilation
    static_assert(EuclideanSpace<E3>);
    static_assert(RiemannianManifold<E3>);
    SUCCEED();
}

TEST_CASE("Euclidean distance", "[euclidean]") {
    E3 space;
    Vec3 a{0.0, 0.0, 0.0};
    Vec3 b{3.0, 4.0, 0.0};
    CHECK_THAT(space.distance(a, b), WithinAbs(5.0, 1e-12));
}

TEST_CASE("Euclidean inner product", "[euclidean]") {
    E3 space;
    Vec3 u{1.0, 0.0, 0.0};
    Vec3 v{0.0, 1.0, 0.0};
    CHECK(space.inner(u, v) == 0.0);  // orthogonal
    CHECK(space.inner(u, u) == 1.0);  // unit vector
}

TEST_CASE("Euclidean norm", "[euclidean]") {
    E3 space;
    Vec3 v{3.0, 4.0, 0.0};
    CHECK_THAT(space.norm(v), WithinAbs(5.0, 1e-12));
}

TEST_CASE("Euclidean exp_map/log_map roundtrip", "[euclidean]") {
    E3 space;
    Vec3 p{1.0, 2.0, 3.0};
    Vec3 q{4.0, 6.0, 8.0};

    // log_map gives the displacement vector
    auto v = space.log_map(p, q);
    CHECK(v[0] == 3.0);
    CHECK(v[1] == 4.0);
    CHECK(v[2] == 5.0);

    // exp_map with t=1 should recover q
    auto q2 = space.exp_map(p, v, 1.0);
    CHECK_THAT(q2[0], WithinAbs(q[0], 1e-12));
    CHECK_THAT(q2[1], WithinAbs(q[1], 1e-12));
    CHECK_THAT(q2[2], WithinAbs(q[2], 1e-12));
}

TEST_CASE("Euclidean exp_map partial t", "[euclidean]") {
    E3 space;
    Vec3 p{0.0, 0.0, 0.0};
    Vec3 v{10.0, 0.0, 0.0};

    auto mid = space.exp_map(p, v, 0.5);
    CHECK_THAT(mid[0], WithinAbs(5.0, 1e-12));
    CHECK(mid[1] == 0.0);
    CHECK(mid[2] == 0.0);
}

TEST_CASE("Euclidean metric_at is constant", "[euclidean]") {
    E3 space;
    Vec3 p1{0.0, 0.0, 0.0};
    Vec3 p2{100.0, 200.0, 300.0};
    Vec3 u{1.0, 0.0, 0.0};
    Vec3 v{0.0, 1.0, 0.0};

    // Flat metric: same at every point
    CHECK(space.metric_at(p1, u, v) == space.metric_at(p2, u, v));
    CHECK(space.metric_at(p1, u, u) == 1.0);
}

TEST_CASE("Euclidean contains everything", "[euclidean]") {
    E3 space;
    CHECK(space.contains(Vec3{0.0, 0.0, 0.0}));
    CHECK(space.contains(Vec3{1e100, -1e100, 0.0}));
}

TEST_CASE("Euclidean works with different dimensions", "[euclidean]") {
    Euclidean<1> e1;
    Vec<double, 1> a{3.0};
    Vec<double, 1> b{7.0};
    CHECK_THAT(e1.distance(a, b), WithinAbs(4.0, 1e-12));

    static_assert(EuclideanSpace<Euclidean<5>>);
    SUCCEED();
}
