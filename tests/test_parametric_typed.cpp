// A ParametricSurface whose function is a type (spaces/parametric.hpp,
// make_parametric): exact partials through a Dual, and geodesics.
//
// The torus has what a geodesic needs checked against something that is not
// the code under test: a meridian is a geodesic (v advances by t|w|/q), so is
// the outer equator (u advances by t|w|/(R+q)), and along any geodesic of a
// surface of revolution Clairaut's quantity rho(v) (T . u_hat) is conserved.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/dual.hpp>
#include <spatium/core/access.hpp>
#include <spatium/spaces/chart.hpp>
#include <spatium/spaces/parametric.hpp>
#include <cmath>
#include <numbers>

using namespace spatium;
using Catch::Matchers::WithinAbs;

namespace {

constexpr double R = 2.0, q = 0.7;
constexpr double two_pi = 2 * std::numbers::pi;

auto typed_torus() {
    return make_parametric<double>(
        [](auto u, auto v) {
            using std::cos; using std::sin;
            return Vec<decltype(u), 3>{(R + q * cos(v)) * cos(u), (R + q * cos(v)) * sin(u), q * sin(v)};
        },
        {0.0, two_pi, 0.0, two_pi}, true, true);
}

ParametricSurface<double> erased_torus() {
    return ParametricSurface<double>(
        [](double u, double v) {
            return Vec<double, 3>{(R + q * std::cos(v)) * std::cos(u), (R + q * std::cos(v)) * std::sin(u), q * std::sin(v)};
        },
        {0.0, two_pi, 0.0, two_pi}, true, true);
}

}  // namespace

TEST_CASE("A generic function makes the partials exact; a std::function keeps finite differences", "[parametric][typed]") {
    const auto typed = typed_torus();
    const auto erased = erased_torus();
    static_assert(decltype(typed)::exact);
    static_assert(!ParametricSurface<double>::exact);

    const double u = 0.6, v = 1.1;
    const double true_area = q * (R + q * std::cos(v));     // |f_u x f_v| = q (R + q cos v)
    CHECK_THAT(typed.area_element(u, v), WithinAbs(true_area, 1e-13));
    // The pair: erased and typed agree to the erased form's difference-quotient accuracy.
    CHECK_THAT(erased.area_element(u, v), WithinAbs(true_area, 1e-6));
    CHECK((typed.evaluate(u, v) - erased.evaluate(u, v)).norm() < 1e-15);
    CHECK((typed.normal_at(u, v) - erased.normal_at(u, v)).norm() < 1e-6);
}

TEST_CASE("A meridian and the outer equator of the torus are geodesics, at their own speeds", "[parametric][typed]") {
    const auto torus = typed_torus();
    const double u0 = 0.6, v0 = 1.1;
    const auto p = torus.evaluate(u0, v0);
    const Vec<double, 3> fv{-q * std::sin(v0) * std::cos(u0), -q * std::sin(v0) * std::sin(u0), q * std::cos(v0)};
    // along the meridian at unit speed 0.9: v advances by 0.9 / q
    const auto along = torus.exp_map(p, Vec<double, 3>{fv * (0.9 / q)}, 1.0);
    CHECK((along - torus.evaluate(u0, v0 + 0.9 / q)).norm() < 1e-12);

    const auto pe = torus.evaluate(0.3, 0.0);
    const Vec<double, 3> eq{-std::sin(0.3), std::cos(0.3), 0.0};
    const auto around = torus.exp_map(pe, Vec<double, 3>{eq * 1.5}, 1.0);
    CHECK((around - torus.evaluate(0.3 + 1.5 / (R + q), 0.0)).norm() < 1e-12);

    // and the intrinsic distance along the meridian is its arc length q * dv
    CHECK_THAT(spaces::distance(torus, p, torus.evaluate(u0, v0 + 0.5)), WithinAbs(q * 0.5, 1e-9));
}

TEST_CASE("exp after log is the identity on the typed torus, and the distance is not the chord", "[parametric][typed]") {
    const auto torus = typed_torus();
    const auto a = torus.evaluate(0.2, 0.3), b = torus.evaluate(0.9, 1.0);
    const auto back = torus.exp_map(a, torus.log_map(a, b), 1.0);
    CHECK((back - b).norm() < 1e-9);
    const double geodesic = spaces::distance(torus, a, b);
    CHECK(geodesic > torus.chord_distance(a, b));              // the chord is a lower bound
    CHECK_THAT(geodesic, WithinAbs(spaces::norm_at(torus, a, torus.log_map(a, b)), 1e-12));
}

TEST_CASE("Clairaut's quantity is conserved along a geodesic of the torus", "[parametric][typed]") {
    const auto torus = typed_torus();
    const double ru = 0.4, rv = 0.9;
    const auto start = torus.evaluate(ru, rv);
    const Vec<double, 3> fu{-(R + q * std::cos(rv)) * std::sin(ru), (R + q * std::cos(rv)) * std::cos(ru), 0};
    const Vec<double, 3> fv{-q * std::sin(rv) * std::cos(ru), -q * std::sin(rv) * std::sin(ru), q * std::cos(rv)};
    const Vec<double, 3> t0{fu * (0.6 / (R + q * std::cos(rv))) + fv * (0.8 / q)};     // unit speed, oblique
    const auto clairaut = [](const Vec<double, 3>& pt, const Vec<double, 3>& tv) {
        const double rho = std::hypot(pt[0], pt[1]);
        return rho * tv.dot(Vec<double, 3>{-pt[1] / rho, pt[0] / rho, 0.0});
    };
    const double c0 = clairaut(start, t0);
    for (const double t : {0.5, 1.0, 2.0, 3.5}) {
        // the velocity at t by a central difference of two flows, a step wide enough to keep the solver's noise small
        const double h = 1e-4;
        const auto ahead = torus.exp_map(start, t0, t + h), behind = torus.exp_map(start, t0, t - h);
        const Vec<double, 3> tv{(ahead - behind) * (0.5 / h)};
        CHECK_THAT(clairaut(torus.exp_map(start, t0, t), tv), WithinAbs(c0, 1e-5));
    }
}

TEST_CASE("On a periodic chart log takes the nearest representative of the other point", "[parametric][typed]") {
    // u = 0.1 and u = 6.2 on a 2 pi chart are 0.18 apart, not 6.1.
    const auto torus = typed_torus();
    const auto a = torus.evaluate(0.1, 0.0), b = torus.evaluate(6.2, 0.0);
    CHECK(spaces::distance(torus, a, b) < 0.6);
    CHECK_THAT(spaces::distance(torus, a, b), WithinAbs((two_pi - 6.2 + 0.1) * (R + q), 1e-8));
}

TEST_CASE("The derivative of the flow in t reaches the surface's exp through a Dual scalar", "[parametric][typed][dual]") {
    // d/dt exp_p(t w) at t = 0 is w: AD through the induced metric, its
    // Christoffel symbols, the extrapolated flow and f itself.
    using D = Dual<double>;
    const auto torus = make_parametric<D>(
        [](auto u, auto v) {
            using std::cos; using std::sin;
            return Vec<decltype(u), 3>{(D(R) + D(q) * cos(v)) * cos(u), (D(R) + D(q) * cos(v)) * sin(u), D(q) * sin(v)};
        },
        {D(0.0), D(two_pi), D(0.0), D(two_pi)}, true, true);
    const auto p = torus.evaluate(D(0.6), D(1.1));
    const Vec<double, 3> fv{-q * std::sin(1.1) * std::cos(0.6), -q * std::sin(1.1) * std::sin(0.6), q * std::cos(1.1)};
    const Vec<D, 3> w{D(fv[0]), D(fv[1]), D(fv[2])};
    D t = D(0.0);
    t.deriv = 1.0;
    const auto e = torus.exp_map(p, w, t);
    CHECK_THAT(e[0].deriv, WithinAbs(fv[0], 1e-7));
    CHECK_THAT(e[1].deriv, WithinAbs(fv[1], 1e-7));
    CHECK_THAT(e[2].deriv, WithinAbs(fv[2], 1e-7));
}

TEST_CASE("A typed surface enters the DSL through its erased form, and the erased form is a retraction", "[parametric][typed]") {
    const auto typed = typed_torus();
    const ParametricSurface<double> erased = chart_of(typed);
    CHECK((erased.evaluate(0.6, 1.1) - typed.evaluate(0.6, 1.1)).norm() < 1e-15);
    static_assert(Chart<decltype(typed)>);
    static_assert(Chart<ParametricSurface<double>>);
    // the erased distance is the chord; the typed has no member distance, so it is |log|
    const auto a = typed.evaluate(0.2, 0.3), b = typed.evaluate(0.9, 1.0);
    CHECK_THAT(erased.distance(a, b), WithinAbs(typed.chord_distance(a, b), 1e-15));
    CHECK(spaces::distance(typed, a, b) > erased.distance(a, b));
}

TEST_CASE("Where the chart degenerates the answer is NaN, not a wrong number", "[parametric][typed][infinity]") {
    // A sphere in latitude-longitude: at the pole f_u vanishes and the induced
    // metric does not invert.
    const auto sphere = make_parametric<double>(
        [](auto u, auto v) {
            using std::cos; using std::sin;
            return Vec<decltype(u), 3>{sin(v) * cos(u), sin(v) * sin(u), cos(v)};
        },
        {0.0, two_pi, 0.0, std::numbers::pi}, true, false);
    const Vec<double, 3> pole{0.0, 0.0, 1.0};
    const auto r = sphere.exp_map(pole, Vec<double, 3>{0.3, 0.0, 0.0}, 1.0);
    CHECK_FALSE(std::isfinite(r[0]));
    // Away from the pole it is the great circle: an eighth of a turn south along a meridian from the equator.
    const auto eq = sphere.evaluate(0.0, std::numbers::pi / 2);
    const auto eighth = sphere.exp_map(eq, Vec<double, 3>{0.0, 0.0, -1.0}, std::numbers::pi / 4);
    CHECK((eighth - sphere.evaluate(0.0, 3 * std::numbers::pi / 4)).norm() < 1e-9);   // south: v grows from the equator
}
