// The integral over a space given by a metric (spaces/integrate_volume.hpp), against
// closed forms that come from nowhere near the code: the torus's area, the whole sphere's
// area and a moment of it through the stereographic chart (two infinite domains at once),
// the area of a hyperbolic disk, the length of a circle, and what the status reports where
// the integral does not exist.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/spaces/integrate_volume.hpp>
#include <cmath>
#include <numbers>

using namespace spatium;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
constexpr double kPi = std::numbers::pi;

// The torus in (u, v): g = diag((R + r cos v)^2, r^2)
struct TorusMetric {
    double R, r;
    template<class S>
    Matrix<S, 2, 2> operator()(const Vec<S, 2>& x) const {
        using std::cos;
        Matrix<S, 2, 2> g{};
        const S ring = S(R) + S(r) * cos(x[1]);
        g(0, 0) = ring * ring;
        g(1, 1) = S(r * r);
        return g;
    }
};

// The unit sphere through the stereographic chart: g = 4 / (1 + |x|^2)^2 I
struct SphereChart {
    template<class S>
    Matrix<S, 2, 2> operator()(const Vec<S, 2>& x) const {
        const S r2 = x[0] * x[0] + x[1] * x[1];
        const S lam = S(2) / (S(1) + r2);
        Matrix<S, 2, 2> g{};
        g(0, 0) = lam * lam;
        g(1, 1) = lam * lam;
        return g;
    }
};

// The hyperbolic plane in geodesic polar coordinates (rho, theta): g = diag(1, sinh^2 rho)
struct HyperbolicPolar {
    template<class S>
    Matrix<S, 2, 2> operator()(const Vec<S, 2>& x) const {
        using std::sinh;
        Matrix<S, 2, 2> g{};
        g(0, 0) = S(1);
        g(1, 1) = sinh(x[0]) * sinh(x[0]);
        return g;
    }
};

// A circle of radius a, parametrised by angle: g = a^2
struct CircleMetric {
    double a;
    template<class S>
    Matrix<S, 1, 1> operator()(const Vec<S, 1>&) const {
        Matrix<S, 1, 1> g{};
        g(0, 0) = S(a * a);
        return g;
    }
};

// A metric under which the plane has infinite area: g = I / (1 + |x|^2)
struct Slow {
    template<class S>
    Matrix<S, 2, 2> operator()(const Vec<S, 2>& x) const {
        const S f = S(1) / (S(1) + x[0] * x[0] + x[1] * x[1]);
        Matrix<S, 2, 2> g{};
        g(0, 0) = f;
        g(1, 1) = f;
        return g;
    }
};
}  // namespace

TEST_CASE("the area of a torus is 4 pi^2 R r", "[integrate_volume]") {
    const double R = 2.0, r = 0.7;
    const auto a = spaces::volume_of(TorusMetric{R, r}, Finite<double>{0.0, 2 * kPi}, Finite<double>{0.0, 2 * kPi});
    CHECK(a.trusted());
    CHECK_THAT(a.value, WithinRel(4 * kPi * kPi * R * r, 1e-10));
}

TEST_CASE("the whole sphere through the stereographic chart: two infinite domains, 4 pi", "[integrate_volume]") {
    // Nothing is sampled at the pole the chart leaves out; the plane is the whole space.
    const auto a = spaces::volume_of(SphereChart{}, WholeLine<double>{}, WholeLine<double>{});
    CHECK(a.trusted());
    CHECK_THAT(a.value, WithinRel(4 * kPi, 1e-8));
    // A moment: the integral of z^2 over the sphere is 4 pi / 3, with z = (1 - r^2) / (1 + r^2)
    const auto z2 = spaces::integrate_volume(SphereChart{},
        [](const Vec<double, 2>& x) {
            const double r2 = x[0] * x[0] + x[1] * x[1];
            const double z = 1 - 2 / (1 + r2);          // not (1 - r2)/(1 + r2): inf/inf at the far nodes of the infinite rule
            return z * z;
        }, WholeLine<double>{}, WholeLine<double>{});
    CHECK(z2.trusted());
    CHECK_THAT(z2.value, WithinRel(4 * kPi / 3, 1e-8));
}

TEST_CASE("the area of a hyperbolic disk is 2 pi (cosh rho - 1)", "[integrate_volume]") {
    for (const double rho : {0.5, 1.0, 2.5}) {
        const auto a = spaces::volume_of(HyperbolicPolar{}, Finite<double>{0.0, rho}, Finite<double>{0.0, 2 * kPi});
        CHECK_THAT(a.value, WithinRel(2 * kPi * (std::cosh(rho) - 1), 1e-10));
    }
}

TEST_CASE("one dimension: the length of a circle", "[integrate_volume]") {
    const auto l = spaces::volume_of(CircleMetric{3.0}, Finite<double>{0.0, 2 * kPi});
    CHECK_THAT(l.value, WithinRel(2 * kPi * 3.0, 1e-12));
}

TEST_CASE("an integral that does not exist is not called converged", "[integrate_volume]") {
    // The plane's area under the metric 1/(1 + r^2) I diverges (logarithmically, in 2-D).
    const auto a = spaces::volume_of(Slow{}, WholeLine<double>{}, WholeLine<double>{});
    CHECK(!a.trusted());
}
