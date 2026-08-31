#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/spaces/product.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <spatium/spaces/sphere.hpp>
#include <spatium/point.hpp>
#include <numbers>

using namespace spatium;
using Catch::Matchers::WithinAbs;

TEST_CASE("ProductSpace Euclidean x Euclidean", "[product]") {
    ProductSpace<Euclidean<1>, Euclidean<1>> r2;
    static_assert(MetricSpace<decltype(r2)>);

    Vec<double, 2> a{0, 0};
    Vec<double, 2> b{3, 4};
    CHECK_THAT(r2.distance(a, b), WithinAbs(5.0, 1e-10));
}

TEST_CASE("ProductSpace contains", "[product]") {
    ProductSpace<Euclidean<1>, Euclidean<1>> r2;
    Vec<double, 2> p{1, 2};
    CHECK(r2.contains(p));
}

TEST_CASE("ProductSpace Sphere x Sphere (Torus)", "[product]") {
    ProductSpace<Sphere<1>, Sphere<1>> torus;
    static_assert(MetricSpace<ProductSpace<Sphere<1>, Sphere<1>>>);

    // S1 points are Vec<double,2> on unit circle
    // Torus point = (cos a, sin a, cos b, sin b) in Vec<double,4>
    Vec<double, 4> a{1, 0, 1, 0};  // (0, 0) on torus
    Vec<double, 4> b{0, 1, 1, 0};  // (pi/2, 0) on torus
    auto d = torus.distance(a, b);
    CHECK_THAT(d, WithinAbs(std::numbers::pi / 2, 1e-6));
}

TEST_CASE("ProductSpace split/join", "[product]") {
    ProductSpace<Euclidean<1>, Euclidean<1>> r2;
    Vec<double, 2> p{3, 7};
    auto f = r2.first(p);
    auto s = r2.second(p);
    CHECK(f[0] == 3.0);
    CHECK(s[0] == 7.0);
    auto joined = r2.join(f, s);
    CHECK(joined == p);
}

TEST_CASE("ProductSpace Cylinder", "[product]") {
    ProductSpace<Euclidean<1>, Sphere<1>> cyl;
    static_assert(MetricSpace<ProductSpace<Euclidean<1>, Sphere<1>>>);
    // E1 point = Vec<double,1>, S1 point = Vec<double,2>
    // Cylinder point = (height, cos theta, sin theta) in Vec<double,3>
    Vec<double, 3> a{0, 1, 0};
    Vec<double, 3> b{1, 1, 0};
    auto d = cyl.distance(a, b);
    CHECK_THAT(d, WithinAbs(1.0, 1e-10));  // only height differs
}
