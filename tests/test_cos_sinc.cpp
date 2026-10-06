// spaces/cos_sinc.hpp: cos(sqrt(s)) and sinc as functions of s, smooth through
// zero, and what they fixed -- the derivative of a geodesic step in its
// tangent at a zero tangent, which a `speed == 0` shortcut dropped for Dual.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/dual.hpp>
#include <spatium/core/precision.hpp>
#include <spatium/spaces/cos_sinc.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <spatium/spaces/hyperbolic.hpp>
#include <spatium/spaces/sphere.hpp>
#include <cmath>
#include <limits>
#include <numbers>

using namespace spatium;
using Catch::Matchers::WithinAbs;

TEST_CASE("cos_sinc agrees with the direct formula on both sides of its threshold", "[cos_sinc]") {
    // The series below the threshold and the closed form above it meet to
    // rounding: values straddling it must differ by no more than eps does.
    const double threshold = 0.1 * std::pow(machine_epsilon<double>(), 1.0 / 6.0);
    for (const double s : {threshold * 0.5, threshold * 0.999, threshold * 1.001, threshold * 2, 0.01, 0.5, 4.0}) {
        const double r = std::sqrt(s);
        const auto circ = cos_sinc_of_square(s);
        CHECK_THAT(circ.cos, WithinAbs(std::cos(r), 1e-15));
        CHECK_THAT(circ.sinc, WithinAbs(std::sin(r) / r, 1e-15));
        const auto hyp = cos_sinc_of_square(s, true);
        CHECK_THAT(hyp.cos, WithinAbs(std::cosh(r), 1e-15 * std::cosh(r)));
        CHECK_THAT(hyp.sinc, WithinAbs(std::sinh(r) / r, 1e-15 * std::sinh(r) / r));
    }
    // And exactly at zero: 1 and 1, no 0/0.
    CHECK(cos_sinc_of_square(0.0).cos == 1.0);
    CHECK(cos_sinc_of_square(0.0).sinc == 1.0);
}

TEST_CASE("The derivative of a geodesic step in its tangent at zero is t, on every space", "[cos_sinc][dual]") {
    // d/d(eps) exp_p(eps w, 1) at eps = 0 is w. The tangent is a Dual of
    // value zero and derivative w; a `== 0` on it compares the value alone,
    // and the old shortcut returned p with derivative 0 on Sphere and
    // Hyperbolic.
    using D = Dual<double>;
    const auto seeded = [](auto w) {
        Vec<D, w.data.size()> v;
        for (std::size_t i = 0; i < v.data.size(); ++i) { v[i] = D(0.0); v[i].deriv = w[i]; }
        return v;
    };
    const D one{1.0};
    {
        Euclidean<3, D> e;
        const Vec<D, 3> p{D(0.1), D(0.2), D(0.3)};
        const auto r = e.exp_map(p, seeded(Vec<double, 3>{0.3, -0.2, 0.5}), one);
        CHECK_THAT(r[0].deriv, WithinAbs(0.3, 1e-14));
        CHECK_THAT(r[2].deriv, WithinAbs(0.5, 1e-14));
    }
    {
        Sphere<2, D> s;
        const Vec<D, 3> p{D(0), D(0), D(1)};
        const auto r = s.exp_map(p, seeded(Vec<double, 3>{0.3, -0.2, 0.0}), one);
        CHECK_THAT(r[0].deriv, WithinAbs(0.3, 1e-14));
        CHECK_THAT(r[1].deriv, WithinAbs(-0.2, 1e-14));
        CHECK_THAT(r[2].deriv, WithinAbs(0.0, 1e-14));
        CHECK_THAT(r[2].value, WithinAbs(1.0, 1e-15));
    }
    {
        Hyperbolic<2, D> h;
        const Vec<D, 3> p{D(1), D(0), D(0)};
        const auto r = h.exp_map(p, seeded(Vec<double, 3>{0.0, 0.3, -0.2}), one);
        CHECK_THAT(r[1].deriv, WithinAbs(0.3, 1e-14));
        CHECK_THAT(r[2].deriv, WithinAbs(-0.2, 1e-14));
        CHECK_THAT(r[0].value, WithinAbs(1.0, 1e-15));
    }
}

TEST_CASE("A zero tangent leaves the point, and the closed forms are unchanged", "[cos_sinc]") {
    const Sphere<2> s;
    const Vec<double, 3> p{0.0, 0.6, 0.8}, zero{0.0, 0.0, 0.0}, v{1.0, 0.0, 0.0};   // v is tangent at p and a unit vector
    CHECK((s.exp_map(p, zero, 1.0) - p).norm() < 1e-15);
    // A quarter turn along a unit tangent: p cos + v sin.
    const auto q = s.exp_map(p, v, std::numbers::pi / 2);
    CHECK_THAT(q[0], WithinAbs(1.0, 1e-14));
    CHECK_THAT(q[1], WithinAbs(0.0, 1e-14));

    const Hyperbolic<2> h;
    const auto o = Hyperbolic<2>::origin();
    CHECK((h.exp_map(o, Vec<double, 3>{0.0, 0.0, 0.0}, 3.0) - o).norm() < 1e-15);
}

#if SPATIUM_HAS_BOOST_MULTIPRECISION
TEST_CASE("cos_sinc keeps fifty digits across its threshold", "[cos_sinc][precision]") {
    // The series is right to Real50's epsilon where it takes over: the same
    // function, with the threshold following the scalar.
    const double threshold = 0.1 * std::pow(machine_epsilon<Real50>(), 1.0 / 6.0);
    for (const double f : {0.5, 0.999, 1.001, 2.0}) {
        const Real50 s = Real50(threshold) * Real50(f);
        const Real50 r = sqrt(s);
        const auto circ = cos_sinc_of_square(s);
        CHECK(static_cast<double>(abs(circ.cos - cos(r)) / Real50("1e-48")) < 1.0);
        CHECK(static_cast<double>(abs(circ.sinc - sin(r) / r) / Real50("1e-48")) < 1.0);
    }
}
#endif

TEST_CASE("A NaN point has no distance: Hyperbolic does not clamp it to zero", "[hyperbolic][infinity]") {
    // max(NaN, 1) is NaN for double and 1 for a Boost number, so over Real50
    // a NaN point sat at distance 0 from the origin. The comparison form
    // keeps the NaN on every scalar.
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const Hyperbolic<2> h;
    const auto o = Hyperbolic<2>::origin();
    CHECK(std::isnan(h.distance(o, Vec<double, 3>{nan, nan, nan})));
#if SPATIUM_HAS_BOOST_MULTIPRECISION
    const Hyperbolic<2, Real50> hr;
    const Real50 rn = Real50(nan);
    const Vec<Real50, 3> bad{rn, rn, rn};
    CHECK(isnan(hr.distance(Hyperbolic<2, Real50>::origin(), bad)));
#endif
}

TEST_CASE("numeric_limits of a Dual are its value type's", "[dual][infinity]") {
    using D = Dual<double>;
    static_assert(std::numeric_limits<D>::is_specialized);
    CHECK(std::isnan(std::numeric_limits<D>::quiet_NaN().value));
    CHECK(std::isinf(std::numeric_limits<D>::infinity().value));
    CHECK(std::numeric_limits<D>::max().value == std::numeric_limits<double>::max());
    CHECK(std::numeric_limits<D>::epsilon().value == std::numeric_limits<double>::epsilon());
    // At any depth, and a derivative of a limit is zero.
    using D2 = Dual<D>;
    CHECK(std::numeric_limits<D2>::max().value.value == std::numeric_limits<double>::max());
    CHECK(std::numeric_limits<D>::max().deriv == 0.0);
}
