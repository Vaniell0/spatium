#include <catch2/catch_test_macros.hpp>

#if SPATIUM_HAS_EIGEN

#include <spatium/algebra/eigen_interop.hpp>
#include <Eigen/Geometry>

using namespace spatium;

TEST_CASE("Vec3 round-trip through Eigen", "[eigen]") {
    Vec3 v{1.0, 2.0, 3.0};
    auto ev = to_eigen(v);
    auto v2 = from_eigen<double, 3>(ev);
    CHECK(v2[0] == 1.0);
    CHECK(v2[1] == 2.0);
    CHECK(v2[2] == 3.0);
}

TEST_CASE("Vec3f round-trip through Eigen", "[eigen]") {
    Vec<float, 3> v{1.0f, 2.5f, 3.5f};
    auto ev = to_eigen(v);
    auto v2 = from_eigen<float, 3>(ev);
    CHECK(v2[0] == 1.0f);
    CHECK(v2[1] == 2.5f);
    CHECK(v2[2] == 3.5f);
}

TEST_CASE("Matrix3x3 round-trip through Eigen", "[eigen]") {
    Matrix<double, 3, 3> m;
    m(0, 0) = 1.0; m(0, 1) = 2.0; m(0, 2) = 3.0;
    m(1, 0) = 4.0; m(1, 1) = 5.0; m(1, 2) = 6.0;
    m(2, 0) = 7.0; m(2, 1) = 8.0; m(2, 2) = 9.0;

    auto em = to_eigen(m);
    auto m2 = from_eigen<double, 3, 3>(em);

    for (std::size_t r = 0; r < 3; ++r)
        for (std::size_t c = 0; c < 3; ++c)
            CHECK(m2(r, c) == m(r, c));
}

TEST_CASE("eigen_view Vec: zero-copy read", "[eigen]") {
    Vec3 v{10.0, 20.0, 30.0};
    auto view = eigen_view(v);
    CHECK(view(0) == 10.0);
    CHECK(view(1) == 20.0);
    CHECK(view(2) == 30.0);
}

TEST_CASE("eigen_view Vec: mutable write-through", "[eigen]") {
    Vec3 v{1.0, 2.0, 3.0};
    auto view = eigen_view(v);
    view(1) = 99.0;
    CHECK(v[1] == 99.0);
}

TEST_CASE("eigen_view Matrix: zero-copy read", "[eigen]") {
    Matrix<double, 2, 2> m;
    m(0, 0) = 1.0; m(0, 1) = 2.0;
    m(1, 0) = 3.0; m(1, 1) = 4.0;

    auto view = eigen_view(m);
    CHECK(view(0, 0) == 1.0);
    CHECK(view(0, 1) == 2.0);
    CHECK(view(1, 0) == 3.0);
    CHECK(view(1, 1) == 4.0);
}

TEST_CASE("eigen_view Matrix: mutable write-through", "[eigen]") {
    Matrix<double, 2, 2> m;
    m(0, 0) = 1.0; m(0, 1) = 2.0;
    m(1, 0) = 3.0; m(1, 1) = 4.0;

    auto view = eigen_view(m);
    view(1, 0) = 42.0;
    CHECK(m(1, 0) == 42.0);
}

TEST_CASE("to_eigen preserves Eigen operations", "[eigen]") {
    Vec3 a{1.0, 0.0, 0.0};
    Vec3 b{0.0, 1.0, 0.0};
    auto ea = to_eigen(a);
    auto eb = to_eigen(b);

    // Eigen dot product
    CHECK(ea.dot(eb) == 0.0);

    // Eigen cross product
    auto ec = ea.cross(eb);
    CHECK(ec(2) == 1.0);
}

TEST_CASE("Non-square matrix round-trip", "[eigen]") {
    Matrix<double, 2, 3> m;
    m(0, 0) = 1.0; m(0, 1) = 2.0; m(0, 2) = 3.0;
    m(1, 0) = 4.0; m(1, 1) = 5.0; m(1, 2) = 6.0;

    auto em = to_eigen(m);
    auto m2 = from_eigen<double, 2, 3>(em);

    for (std::size_t r = 0; r < 2; ++r)
        for (std::size_t c = 0; c < 3; ++c)
            CHECK(m2(r, c) == m(r, c));
}

#endif // SPATIUM_HAS_EIGEN
