#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/polynomial.hpp>
#include <spatium/core/precision.hpp>
#include <algorithm>

using namespace spatium;
using Catch::Matchers::WithinAbs;

// ── Quadratic ────────────────────────────────────────────────

TEST_CASE("Quadratic: two real roots", "[polynomial]") {
    // x² - 5x + 6 = 0 → x = 2, 3
    auto roots = solve_quadratic(1.0, -5.0, 6.0);
    CHECK_THAT(roots[0].re, WithinAbs(3.0, 1e-10));
    CHECK_THAT(roots[1].re, WithinAbs(2.0, 1e-10));
    CHECK(roots[0].is_real());
    CHECK(roots[1].is_real());
}

TEST_CASE("Quadratic: repeated root", "[polynomial]") {
    // x² - 4x + 4 = 0 → x = 2, 2
    auto roots = solve_quadratic(1.0, -4.0, 4.0);
    CHECK_THAT(roots[0].re, WithinAbs(2.0, 1e-10));
    CHECK_THAT(roots[1].re, WithinAbs(2.0, 1e-10));
}

TEST_CASE("Quadratic: complex roots", "[polynomial]") {
    // x² + 1 = 0 → x = ±i
    auto roots = solve_quadratic(1.0, 0.0, 1.0);
    CHECK_THAT(roots[0].re, WithinAbs(0.0, 1e-10));
    CHECK_THAT(roots[0].im, WithinAbs(1.0, 1e-10));
    CHECK_THAT(roots[1].re, WithinAbs(0.0, 1e-10));
    CHECK_THAT(roots[1].im, WithinAbs(-1.0, 1e-10));
    CHECK_FALSE(roots[0].is_real());
}

TEST_CASE("Quadratic: real_roots_quadratic filters correctly", "[polynomial]") {
    auto reals = real_roots_quadratic(1.0, -5.0, 6.0);
    CHECK(reals.size() == 2);

    auto complex = real_roots_quadratic(1.0, 0.0, 1.0);
    CHECK(complex.empty());
}

// ── Cubic ────────────────────────────────────────────────────

TEST_CASE("Cubic: three real roots", "[polynomial]") {
    // (x-1)(x-2)(x-3) = x³ - 6x² + 11x - 6
    auto roots = solve_cubic(1.0, -6.0, 11.0, -6.0);
    std::vector<double> reals;
    for (auto& r : roots) {
        CHECK(r.is_real(1e-8));
        reals.push_back(r.re);
    }
    std::sort(reals.begin(), reals.end());
    CHECK_THAT(reals[0], WithinAbs(1.0, 1e-8));
    CHECK_THAT(reals[1], WithinAbs(2.0, 1e-8));
    CHECK_THAT(reals[2], WithinAbs(3.0, 1e-8));
}

TEST_CASE("Cubic: one real + two complex", "[polynomial]") {
    // x³ + 1 = 0 → x = -1, (1±i√3)/2
    auto roots = solve_cubic(1.0, 0.0, 0.0, 1.0);
    int real_count = 0;
    for (auto& r : roots) {
        if (r.is_real(1e-8)) {
            CHECK_THAT(r.re, WithinAbs(-1.0, 1e-8));
            real_count++;
        }
    }
    CHECK(real_count == 1);
}

TEST_CASE("Cubic: repeated root", "[polynomial]") {
    // (x-1)²(x+2) = x³ - 3x + 2 → wait: x³ + 0x² - 3x + 2
    // Actually: (x-1)²(x+2) = x³ - 0x² - 3x + 2
    auto roots = solve_cubic(1.0, 0.0, -3.0, 2.0);
    std::vector<double> reals;
    for (auto& r : roots) {
        CHECK(r.is_real(1e-6));
        reals.push_back(r.re);
    }
    std::sort(reals.begin(), reals.end());
    CHECK_THAT(reals[0], WithinAbs(-2.0, 1e-6));
    CHECK_THAT(reals[1], WithinAbs(1.0, 1e-6));
    CHECK_THAT(reals[2], WithinAbs(1.0, 1e-6));
}

TEST_CASE("real_roots_cubic filters correctly", "[polynomial]") {
    auto reals = real_roots_cubic(1.0, -6.0, 11.0, -6.0);
    CHECK(reals.size() == 3);
}

// ── Quartic ──────────────────────────────────────────────────

TEST_CASE("Quartic: biquadratic with real roots", "[polynomial]") {
    // x⁴ - 5x² + 4 = 0 → (x²-1)(x²-4) → x = ±1, ±2
    auto roots = solve_quartic(1.0, 0.0, -5.0, 0.0, 4.0);
    std::vector<double> reals;
    for (auto& r : roots) {
        if (r.is_real(1e-6)) reals.push_back(r.re);
    }
    std::sort(reals.begin(), reals.end());
    REQUIRE(reals.size() == 4);
    CHECK_THAT(reals[0], WithinAbs(-2.0, 1e-6));
    CHECK_THAT(reals[1], WithinAbs(-1.0, 1e-6));
    CHECK_THAT(reals[2], WithinAbs(1.0, 1e-6));
    CHECK_THAT(reals[3], WithinAbs(2.0, 1e-6));
}

TEST_CASE("Quartic: all complex", "[polynomial]") {
    // x⁴ + 1 = 0 → no real roots
    auto reals = real_roots_quartic(1.0, 0.0, 0.0, 0.0, 1.0);
    CHECK(reals.empty());
}

TEST_CASE("Quartic: (x-1)(x-2)(x-3)(x-4)", "[polynomial]") {
    // x⁴ - 10x³ + 35x² - 50x + 24
    auto roots = solve_quartic(1.0, -10.0, 35.0, -50.0, 24.0);
    std::vector<double> reals;
    for (auto& r : roots) {
        if (r.is_real(1e-4)) reals.push_back(r.re);
    }
    std::sort(reals.begin(), reals.end());
    REQUIRE(reals.size() == 4);
    CHECK_THAT(reals[0], WithinAbs(1.0, 1e-4));
    CHECK_THAT(reals[1], WithinAbs(2.0, 1e-4));
    CHECK_THAT(reals[2], WithinAbs(3.0, 1e-4));
    CHECK_THAT(reals[3], WithinAbs(4.0, 1e-4));
}

// ── Real50 (arbitrary precision) ───────────────────────────────
// Every solver here is a plain `template<Scalar T>`, implying it should work
// for any Scalar, Real50 included. It didn't: polynomial.hpp and complex.hpp
// called std::sqrt/cbrt/cos/acos/atan2/abs *qualified*, which only overloads
// for the built-in floating-point types -- Real50 (Boost.Multiprecision)
// provides its own via ADL, found only by unqualified calls. solve_cubic's
// casus-irreducibilis branch also used std::numbers::pi_v<T>, which the
// standard restricts to std::floating_point types outright; replaced with
// the portable 4*atan(1). Found while building an RSC dispatch domain meant
// to compare Spatium's own double vs. Real50 precision -- not a hypothetical
// bug, solve_cubic<Real50> failed to compile at all before this fix, on
// every branch.

TEST_CASE("Quadratic with Real50", "[polynomial]") {
    auto roots = solve_quadratic(Real50{1}, Real50{-5}, Real50{6});
    CHECK(abs(roots[0].re - Real50{3}) < Real50{1e-30});
    CHECK(abs(roots[1].re - Real50{2}) < Real50{1e-30});
}

TEST_CASE("Cubic with Real50: one real root, two complex", "[polynomial]") {
    // x³ - 1 = 0 -> disc > 0 branch
    auto roots = solve_cubic(Real50{1}, Real50{0}, Real50{0}, Real50{-1});
    CHECK(roots[0].is_real(Real50{1e-25}));
    CHECK(abs(roots[0].re - Real50{1}) < Real50{1e-25});
}

TEST_CASE("Cubic with Real50: repeated root", "[polynomial]") {
    // (x-2)^3 = x³ - 6x² + 12x - 8 -> |disc| ~ 0 branch
    auto roots = solve_cubic(Real50{1}, Real50{-6}, Real50{12}, Real50{-8});
    for (auto& r : roots) CHECK(abs(r.re - Real50{2}) < Real50{1e-25});
}

TEST_CASE("Cubic with Real50: casus irreducibilis (the branch that was broken)",
          "[polynomial]") {
    // (x-1)(x-2)(x-3) -> three distinct real roots, the trig-formula branch
    auto roots = solve_cubic(Real50{1}, Real50{-6}, Real50{11}, Real50{-6});
    std::vector<Real50> reals;
    for (auto& r : roots)
        if (r.is_real(Real50{1e-25})) reals.push_back(r.re);
    REQUIRE(reals.size() == 3);
    std::sort(reals.begin(), reals.end());
    // acos/cos are transcendental -- exact equality isn't reasonable to
    // expect even at 50 decimal digits, tolerance is.
    CHECK(abs(reals[0] - Real50{1}) < Real50{1e-25});
    CHECK(abs(reals[1] - Real50{2}) < Real50{1e-25});
    CHECK(abs(reals[2] - Real50{3}) < Real50{1e-25});
}

TEST_CASE("Quartic with Real50", "[polynomial]") {
    // (x-1)(x-2)(x-3)(x-4)
    auto roots = solve_quartic(Real50{1}, Real50{-10}, Real50{35}, Real50{-50}, Real50{24});
    std::vector<Real50> reals;
    for (auto& r : roots)
        if (r.is_real(Real50{1e-20})) reals.push_back(r.re);
    REQUIRE(reals.size() == 4);
    std::sort(reals.begin(), reals.end());
    CHECK(abs(reals[0] - Real50{1}) < Real50{1e-20});
    CHECK(abs(reals[1] - Real50{2}) < Real50{1e-20});
    CHECK(abs(reals[2] - Real50{3}) < Real50{1e-20});
    CHECK(abs(reals[3] - Real50{4}) < Real50{1e-20});
}
