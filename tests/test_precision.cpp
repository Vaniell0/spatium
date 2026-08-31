#include <catch2/catch_test_macros.hpp>
#include <spatium/core/precision.hpp>
#include <spatium/algebra/vector.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <spatium/geometry/triangle.hpp>

using namespace spatium;
using namespace spatium::geometry;

TEST_CASE("Real50 satisfies Scalar", "[precision]") {
    static_assert(Scalar<Real50>);
    static_assert(Scalar<Real100>);
    SUCCEED();
}

TEST_CASE("Vec with Real50", "[precision]") {
    Vec<Real50, 3> a{Real50{1}, Real50{2}, Real50{3}};
    Vec<Real50, 3> b{Real50{4}, Real50{5}, Real50{6}};
    auto c = a + b;
    CHECK(c[0] == Real50{5});
    CHECK(c[1] == Real50{7});
    CHECK(c[2] == Real50{9});
}

TEST_CASE("Euclidean with Real50", "[precision]") {
    Euclidean<3, Real50> space;
    Vec<Real50, 3> a{Real50{0}, Real50{0}, Real50{0}};
    Vec<Real50, 3> b{Real50{3}, Real50{4}, Real50{0}};
    auto d = space.distance(a, b);
    CHECK(d == Real50{5});
}

TEST_CASE("Triangle with Real50 — high precision area", "[precision]") {
    // Right triangle with known area
    Triangle<3, Real50> tri(
        Vec<Real50, 3>{Real50{0}, Real50{0}, Real50{0}},
        Vec<Real50, 3>{Real50{1}, Real50{0}, Real50{0}},
        Vec<Real50, 3>{Real50{0}, Real50{1}, Real50{0}}
    );
    auto area = tri.area();
    // 0.5 with 50 digits of precision
    CHECK(area == Real50{"0.5"});
}

TEST_CASE("High precision pi approximation via polygon", "[precision]") {
    // Inscribed regular N-gon in unit circle: area = N/2 * sin(2π/N)
    // As N→∞, area → π
    // With Real100, we can compute this very precisely

    using T = Real100;
    using boost::multiprecision::asin;

    // π via asin(1)*2
    T pi = asin(T{1}) * T{2};

    // 10000-gon inscribed in unit circle
    int n = 10000;
    T angle = T{2} * pi / T{n};
    T polygon_area = T{n} / T{2} * sin(angle);

    // Should be close to π
    T error = abs(polygon_area - pi);
    CHECK(error < T{"0.00001"});
}
