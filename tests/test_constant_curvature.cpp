// The spaces of constant curvature as one family (spaces/constant_curvature.hpp).
//
// Each check is a pair of paths to one answer that share nothing: the family at
// kappa = 0, 1, -1 against the closed forms the library already has (Euclidean,
// Sphere, Hyperbolic) through their stereographic maps; the family at kappa
// against itself at 1 under the scaling law; and the family against itself
// across kappa = 0, where every formula is 0/0.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/core/precision.hpp>
#include <spatium/core/verify.hpp>
#include <spatium/algebra/dual.hpp>
#include <spatium/spaces/constant_curvature.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <spatium/spaces/hyperbolic.hpp>
#include <spatium/spaces/sphere.hpp>
#include <cmath>
#include <numbers>
#include <vector>

using namespace spatium;
using Catch::Matchers::WithinAbs;

namespace {

using V2 = Vec<double, 2>;
const V2 a{0.12, -0.31}, b{0.40, 0.22}, c{-0.25, 0.18};

// Unit sphere S^2 from the stereographic coordinates of kappa = 1.
Vec<double, 3> to_sphere(const V2& x) {
    const double r2 = x.norm_squared();
    return {2 * x[0] / (1 + r2), 2 * x[1] / (1 + r2), (1 - r2) / (1 + r2)};
}
// The hyperboloid (time first) from the Poincare ball of kappa = -1.
Vec<double, 3> to_hyperboloid(const V2& x) {
    const double r2 = x.norm_squared();
    return {(1 + r2) / (1 - r2), 2 * x[0] / (1 - r2), 2 * x[1] / (1 - r2)};
}

}  // namespace

TEST_CASE("kappa = 0 is the Euclidean plane, in the same code", "[constant_curvature]") {
    const ConstantCurvature<2> flat{0.0};
    const Euclidean<2> e;
    // the conformal factor is 2, so the metric is 4 <u, v>: distance is twice Euclid's
    CHECK_THAT(flat.distance(a, b), WithinAbs(2 * e.distance(a, b), 1e-14));
    // exp and log: x + t v * (lambda/2) with lambda = 2, so a unit-coordinate step goes one unit
    const V2 v{0.3, -0.2};
    const V2 x = flat.exp_map(a, v, 1.0);
    CHECK_THAT(x[0], WithinAbs(a[0] + v[0], 1e-14));
    CHECK_THAT(x[1], WithinAbs(a[1] + v[1], 1e-14));
    const V2 back = flat.log_map(a, b);
    CHECK_THAT(back[0], WithinAbs(b[0] - a[0], 1e-14));
    CHECK_THAT(back[1], WithinAbs(b[1] - a[1], 1e-14));
}

TEST_CASE("kappa = 1 is Sphere<2>, through the stereographic map", "[constant_curvature]") {
    const ConstantCurvature<2> family{1.0};
    const Sphere<2> sphere;
    for (const auto& [p, q] : {std::pair{a, b}, std::pair{b, c}, std::pair{a, c}})
        CHECK_THAT(family.distance(p, q), WithinAbs(sphere.distance(to_sphere(p), to_sphere(q)), 1e-13));
}

TEST_CASE("kappa = -1 is Hyperbolic<2>, through the Poincare ball", "[constant_curvature]") {
    const ConstantCurvature<2> family{-1.0};
    const Hyperbolic<2> hyper;
    for (const auto& [p, q] : {std::pair{a, b}, std::pair{b, c}, std::pair{a, c}})
        CHECK_THAT(family.distance(p, q), WithinAbs(hyper.distance(to_hyperboloid(p), to_hyperboloid(q)), 1e-12));
}

TEST_CASE("scaling: d_kappa(x, y) = d_1(sqrt(kappa) x, sqrt(kappa) y) / sqrt(kappa)", "[constant_curvature][symmetry]") {
    const ConstantCurvature<2> unit{1.0};
    for (const double kappa : {0.25, 2.0, 9.0}) {
        const ConstantCurvature<2> family{kappa};
        const double s = std::sqrt(kappa);
        CHECK_THAT(family.distance(a, b), WithinAbs(unit.distance(a * s, b * s) / s, 1e-13));
    }
    // and the same for the negative side, against the unit hyperbolic ball
    const ConstantCurvature<2> hyp{-1.0};
    for (const double kappa : {-0.25, -2.0}) {
        const ConstantCurvature<2> family{kappa};
        const double s = std::sqrt(-kappa);
        CHECK_THAT(family.distance(a * 0.4, b * 0.4), WithinAbs(hyp.distance(a * 0.4 * s, b * 0.4 * s) / s, 1e-13));
    }
}

TEST_CASE("exp and log invert each other and the step has the length it was asked", "[constant_curvature]") {
    for (const double kappa : {-1.0, -0.3, 0.0, 0.3, 1.0}) {
        const ConstantCurvature<2> family{kappa};
        const V2 v{0.35, -0.2};
        const V2 y = family.exp_map(a, v, 1.0);
        REQUIRE(family.contains(y));
        const V2 back = family.log_map(a, y);
        CHECK_THAT(back[0], WithinAbs(v[0], 1e-12));
        CHECK_THAT(back[1], WithinAbs(v[1], 1e-12));
        // the geodesic's length is the metric norm of v: lambda |v|
        const double len = std::sqrt(family.metric_at(a, v, v));
        CHECK_THAT(family.distance(a, y), WithinAbs(len, 1e-12));
    }
}

TEST_CASE("the axioms hold at kappa of either sign and zero", "[constant_curvature]") {
    const std::vector<V2> pts{a, b, c, V2{0.0, 0.0}, V2{-0.1, -0.4}};
    for (const double kappa : {-1.0, 0.0, 1.0}) {
        const ConstantCurvature<2> family{kappa};
        const auto r = verify_metric(family, pts, 1e-9);
        INFO("kappa " << kappa << ": " << r.message());
        CHECK(r.passed);
        const auto e = verify_exp_log(family, pts, 1e-8);
        INFO("kappa " << kappa << ": " << e.message());
        CHECK(e.passed);
    }
}

TEST_CASE("continuity through kappa = 0, where every formula is 0/0", "[constant_curvature][limit]") {
    // The probe the matrix lacks: the function of kappa has no jump at 0, on either side,
    // at every scale of kappa a double can tell apart from 0.
    const double d0 = ConstantCurvature<2>{0.0}.distance(a, b);
    for (const double k : {1e-3, 1e-6, 1e-9, 1e-12, 1e-15, 1e-30, 1e-100}) {
        const double dp = ConstantCurvature<2>{k}.distance(a, b);
        const double dm = ConstantCurvature<2>{-k}.distance(a, b);
        INFO("kappa " << k);
        CHECK(std::isfinite(dp));
        CHECK(std::isfinite(dm));
        CHECK(std::abs(dp - d0) <= 4 * k + 1e-14);      // d is smooth in kappa: the gap is O(kappa)
        CHECK(std::abs(dm - d0) <= 4 * k + 1e-14);
    }
    // exp and log too
    const V2 v{0.2, 0.1};
    for (const double k : {1e-9, -1e-9, 1e-20}) {
        const ConstantCurvature<2> f{k}, flat{0.0};
        const V2 x = f.exp_map(a, v, 1.0), x0 = flat.exp_map(a, v, 1.0);
        CHECK_THAT(x[0], WithinAbs(x0[0], 1e-8));
        CHECK_THAT(x[1], WithinAbs(x0[1], 1e-8));
        const V2 l = f.log_map(a, b), l0 = flat.log_map(a, b);
        CHECK_THAT(l[0], WithinAbs(l0[0], 1e-8));
    }
}

TEST_CASE("the derivative in kappa: Dual against a central difference, through zero", "[constant_curvature][dual]") {
    using D = Dual<double>;
    const auto dist_at = [](double k) { return ConstantCurvature<2>{k}.distance(a, b); };
    const Vec<D, 2> da{D(a[0]), D(a[1])}, db{D(b[0]), D(b[1])};
    for (const double k0 : {-0.7, 0.0, 0.7}) {
        const ConstantCurvature<2, D> family{D::variable(k0)};
        const double derivative = family.distance(da, db).deriv;
        const double h = 1e-5;
        const double fd = (dist_at(k0 + h) - dist_at(k0 - h)) / (2 * h);
        INFO("kappa " << k0);
        CHECK_THAT(derivative, WithinAbs(fd, 1e-7));
    }
}

#if SPATIUM_HAS_BOOST_MULTIPRECISION
TEST_CASE("fifty digits across kappa = 0", "[constant_curvature][precision]") {
    using R = Real50;
    const Vec<R, 2> x{R("0.12"), R("-0.31")}, y{R("0.4"), R("0.22")};
    const R d0 = ConstantCurvature<2, R>{R(0)}.distance(x, y);
    for (const char* k : {"1e-20", "1e-40", "-1e-40"}) {
        const R dk = ConstantCurvature<2, R>{R(k)}.distance(x, y);
        CHECK(abs(dk - d0) < R(k[0] == '-' ? "1e-39" : k) * R(8));
    }
}
#endif
