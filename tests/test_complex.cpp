#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/complex.hpp>
#include <numbers>

using namespace spatium;
using Catch::Matchers::WithinAbs;

TEST_CASE("Complex default is zero", "[complex]") {
    Complex64 z;
    CHECK(z.re == 0.0);
    CHECK(z.im == 0.0);
}

TEST_CASE("Complex from real", "[complex]") {
    Complex64 z{3.0};
    CHECK(z.re == 3.0);
    CHECK(z.im == 0.0);
    CHECK(z.is_real());
}

TEST_CASE("Complex addition", "[complex]") {
    Complex64 a{1, 2}, b{3, 4};
    auto c = a + b;
    CHECK(c.re == 4.0);
    CHECK(c.im == 6.0);
}

TEST_CASE("Complex subtraction", "[complex]") {
    Complex64 a{5, 3}, b{2, 1};
    auto c = a - b;
    CHECK(c.re == 3.0);
    CHECK(c.im == 2.0);
}

TEST_CASE("Complex multiplication", "[complex]") {
    Complex64 a{1, 2}, b{3, 4};
    auto c = a * b;
    // (1+2i)(3+4i) = 3+4i+6i+8i² = 3+10i-8 = -5+10i
    CHECK_THAT(c.re, WithinAbs(-5.0, 1e-12));
    CHECK_THAT(c.im, WithinAbs(10.0, 1e-12));
}

TEST_CASE("Complex division", "[complex]") {
    Complex64 a{-5, 10}, b{3, 4};
    auto c = a / b;
    // Should recover (1, 2)
    CHECK_THAT(c.re, WithinAbs(1.0, 1e-12));
    CHECK_THAT(c.im, WithinAbs(2.0, 1e-12));
}

TEST_CASE("Complex conjugate", "[complex]") {
    Complex64 z{3, 4};
    auto c = z.conjugate();
    CHECK(c.re == 3.0);
    CHECK(c.im == -4.0);
}

TEST_CASE("Complex magnitude", "[complex]") {
    Complex64 z{3, 4};
    CHECK_THAT(z.magnitude(), WithinAbs(5.0, 1e-12));
    CHECK_THAT(z.magnitude_sq(), WithinAbs(25.0, 1e-12));
}

TEST_CASE("Complex phase", "[complex]") {
    Complex64 z{1, 1};
    CHECK_THAT(z.phase(), WithinAbs(std::numbers::pi / 4.0, 1e-12));
}

TEST_CASE("Complex from_polar", "[complex]") {
    auto z = Complex64::from_polar(5.0, std::numbers::pi / 2.0);
    CHECK_THAT(z.re, WithinAbs(0.0, 1e-12));
    CHECK_THAT(z.im, WithinAbs(5.0, 1e-12));
}

TEST_CASE("Complex sqrt", "[complex]") {
    // sqrt(-1) = i
    auto z = sqrt(Complex64{-1.0, 0.0});
    CHECK_THAT(z.re, WithinAbs(0.0, 1e-12));
    CHECK_THAT(z.im, WithinAbs(1.0, 1e-12));
}

TEST_CASE("Complex sqrt positive real", "[complex]") {
    auto z = sqrt(Complex64{4.0, 0.0});
    CHECK_THAT(z.re, WithinAbs(2.0, 1e-12));
    CHECK_THAT(z.im, WithinAbs(0.0, 1e-12));
}

TEST_CASE("Complex is_real", "[complex]") {
    CHECK(Complex64{3.0, 0.0}.is_real());
    CHECK_FALSE(Complex64{3.0, 1.0}.is_real());
    CHECK(Complex64{3.0, 1e-16}.is_real());
}

TEST_CASE("Complex scalar multiplication", "[complex]") {
    Complex64 z{2, 3};
    auto a = z * 2.0;
    auto b = 2.0 * z;
    CHECK(a == b);
    CHECK(a.re == 4.0);
    CHECK(a.im == 6.0);
}

TEST_CASE("Complex constexpr arithmetic", "[complex]") {
    constexpr Complex64 a{1, 2};
    constexpr Complex64 b{3, 4};
    constexpr auto c = a + b;
    static_assert(c.re == 4.0);
    static_assert(c.im == 6.0);
    SUCCEED();
}
