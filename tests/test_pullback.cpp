// The pullback of a metric through a map (spaces/pullback.hpp).
//
// Each check sets the pulled-back metric against something that does not come from
// it: the torus's first fundamental form in closed form; the kappa family's conformal
// factor, for the sphere through the stereographic map and for the hyperbolic plane
// through the Poincare map (the second from a Lorentzian ambient metric); and the
// distance of the space it defines, derived by geodesic shooting, against the closed
// form of the family.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/core/access.hpp>
#include <spatium/spaces/constant_curvature.hpp>
#include <spatium/spaces/metric_chart.hpp>
#include <spatium/spaces/pullback.hpp>
#include <cmath>

using namespace spatium;
using Catch::Matchers::WithinAbs;

namespace {

// (u, v) -> a torus of radius R, tube r
struct Torus {
    double R, r;
    template<class S>
    Vec<S, 3> operator()(const Vec<S, 2>& x) const {
        using std::cos; using std::sin;
        const S ring = S(R) + S(r) * cos(x[1]);
        return {ring * cos(x[0]), ring * sin(x[0]), S(r) * sin(x[1])};
    }
};

// x in R^2 -> the unit sphere in R^3, stereographic from the south pole
struct Stereographic {
    template<class S>
    Vec<S, 3> operator()(const Vec<S, 2>& x) const {
        const S r2 = x[0] * x[0] + x[1] * x[1];
        return {S(2) * x[0] / (S(1) + r2), S(2) * x[1] / (S(1) + r2), (S(1) - r2) / (S(1) + r2)};
    }
};

// x in the unit disk -> the hyperboloid (time first)
struct PoincareToHyperboloid {
    template<class S>
    Vec<S, 3> operator()(const Vec<S, 2>& x) const {
        const S r2 = x[0] * x[0] + x[1] * x[1];
        return {(S(1) + r2) / (S(1) - r2), S(2) * x[0] / (S(1) - r2), S(2) * x[1] / (S(1) - r2)};
    }
};

}  // namespace

TEST_CASE("the pullback of R^3 through a torus is its first fundamental form", "[pullback]") {
    const double R = 2.0, r = 0.7;
    const auto g = spaces::pullback_metric<2, 3>(Torus{R, r});
    for (const auto& x : {Vec<double, 2>{0.3, 0.2}, Vec<double, 2>{1.7, -2.1}, Vec<double, 2>{-0.9, 3.0}}) {
        const auto m = g(x);
        const double ring = R + r * std::cos(x[1]);
        CHECK_THAT(m(0, 0), WithinAbs(ring * ring, 1e-13));
        CHECK_THAT(m(1, 1), WithinAbs(r * r, 1e-13));
        CHECK_THAT(m(0, 1), WithinAbs(0.0, 1e-13));
        CHECK_THAT(m(1, 0), WithinAbs(0.0, 1e-13));
    }
}

TEST_CASE("the sphere through the stereographic map is the kappa = 1 family's metric", "[pullback][symmetry]") {
    const auto g = spaces::pullback_metric<2, 3>(Stereographic{});
    const ConstantCurvature<2> family{1.0};
    for (const auto& x : {Vec<double, 2>{0.3, 0.2}, Vec<double, 2>{-1.2, 0.8}, Vec<double, 2>{2.5, -1.0}}) {
        const auto m = g(x);
        const Vec<double, 2> e0{1, 0}, e1{0, 1};
        CHECK_THAT(m(0, 0), WithinAbs(family.metric_at(x, e0, e0), 1e-12));
        CHECK_THAT(m(1, 1), WithinAbs(family.metric_at(x, e1, e1), 1e-12));
        CHECK_THAT(m(0, 1), WithinAbs(0.0, 1e-12));
    }
}

TEST_CASE("the hyperbolic plane through the Poincare map, from a Lorentzian ambient metric", "[pullback][symmetry]") {
    const auto g = spaces::pullback_metric<2, 3>(PoincareToHyperboloid{}, spaces::MinkowskiMetric{});
    const ConstantCurvature<2> family{-1.0};
    for (const auto& x : {Vec<double, 2>{0.3, 0.2}, Vec<double, 2>{-0.5, 0.6}, Vec<double, 2>{0.1, -0.8}}) {
        const auto m = g(x);
        const Vec<double, 2> e0{1, 0}, e1{0, 1};
        CHECK_THAT(m(0, 0), WithinAbs(family.metric_at(x, e0, e0), 1e-12));
        CHECK_THAT(m(1, 1), WithinAbs(family.metric_at(x, e1, e1), 1e-12));
        CHECK_THAT(m(0, 1), WithinAbs(0.0, 1e-12));
    }
}

TEST_CASE("the space a pullback defines has geodesics: shooting against the closed form", "[pullback][symmetry]") {
    // The distance of the pulled-back sphere, found by the geodesic equation, Christoffel
    // symbols by Dual and an extrapolated flow, against the kappa family's closed form.
    const auto chart = spaces::metric_chart<double, 2>(spaces::pullback_metric<2, 3>(Stereographic{}));
    const ConstantCurvature<2> family{1.0};
    const Vec<double, 2> a{0.10, 0.20}, b{0.45, -0.30};
    CHECK_THAT(spaces::distance(chart, a, b), WithinAbs(family.distance(a, b), 1e-7));
}
