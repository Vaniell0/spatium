// Continuous collision held to the one rule that makes it a guarantee:
// a reported time of impact may be early, never late, and a miss is
// reported only when it is proved.
//
// The spike and the tube are cases the code before them got wrong, run
// against that code and seen to fail before the fix was written -- a test
// that passes on the broken implementation is not evidence of anything.
// The graze and the units passed on the old code too and stay as guards on
// the edges the fix had to touch. The rest hold the new path to the same
// rule: a degenerate chart point, a ball instead of a point, and a fuzz
// against the closed forms.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <spatium/geometry/ray_surface.hpp>
#include <spatium/physics/mechanics/rigid_contact.hpp>
#include <spatium/spaces/parametric.hpp>

#include <algorithm>
#include <cmath>
#include <format>
#include <random>

using namespace spatium;
using namespace spatium::physics::mechanics;
using V3 = Vec<double, 3>;

namespace {

// A plane with a tall, thin Gaussian spike: height h, width sigma, centred
// at (c, c) over [-1, 1]^2 -- between the nodes of the 17 x 17 grid a
// closest-point search seeds from, so no seed sees it.
ParametricSurface<double> spike(double h, double sigma, double c) {
    return ParametricSurface<double>(
        [h, sigma, c](double u, double v) {
            const double du = u - c, dv = v - c;
            return V3{u, v, h * std::exp(-(du * du + dv * dv) / (sigma * sigma))};
        },
        {-1.0, 1.0, -1.0, 1.0});
}

// A sphere-of-radius-r chart in angles with its per-cell reach, which a
// single constant overstates near the poles: |f_u| = r sin v there.
LipschitzChart<double> sphere_chart_bound(const ParametricSurface<double>& chart, double r) {
    return {chart, r, [r](double u0, double u1, double v0, double v1) {
                const double pi2 = std::numbers::pi / 2;
                const double s = (v0 <= pi2 && pi2 <= v1) ? 1.0
                                                          : std::max(std::sin(v0), std::sin(v1));
                return r * s * (u1 - u0) / 2 + r * (v1 - v0) / 2;
            }};
}

ParametricSurface<double> sphere_chart(double r) {
    return ParametricSurface<double>(
        [r](double u, double v) {
            return V3{r * std::sin(v) * std::cos(u), r * std::sin(v) * std::sin(u), r * std::cos(v)};
        },
        {0.0, 2.0 * std::numbers::pi, 0.0, std::numbers::pi}, true, false);
}

}  // namespace

// Flying sideways through the spike's shaft at mid-height. Every seed of
// the closest-point search sits on the plane a whole unit below, so the
// search reports the plane's distance (1.0) where the shaft is 0.28 away,
// and a step sized by that overshoots the shaft entirely: the code before
// this test reported no contact at all. The shaft's radius at a height is
// closed form, so the true first contact is known.
TEST_CASE("A sweep does not fly through a spike its closest-point search cannot see",
          "[physics][ccd][certified]") {
    const double h = 2.0, sigma = 0.005, c = 0.06, z = 1.0;
    const auto surf = spike(h, sigma, c);
    const double r_shaft = sigma * std::sqrt(std::log(h / z));
    const V3 start{c + 0.3, c, z}, disp{-0.6, 0.0, 0.0};
    const double toi_true = (0.3 - r_shaft) / 0.6;

    // The spike's steepest slope is h * sqrt(2/e) / sigma; a chart is
    // Lipschitz in its parameters with at most sqrt(1 + slope^2).
    const double L = std::sqrt(1.0 + std::pow(h * std::sqrt(2.0 / std::numbers::e) / sigma, 2));
    const auto s = sweep_point_surface<double>(start, disp, LipschitzChart<double>{surf, L});
    INFO(std::format("true toi {:.6f}, reported hit={} toi={:.6f}", toi_true, s.hit, s.toi));
    REQUIRE(s.hit);
    CHECK(s.toi <= toi_true + 1e-9);
    CHECK(s.toi > toi_true - 0.05);   // early is safe; this early is still useful
}

// A thin tube, the point three units away. The code before this test
// reported a contact at time zero, because its inside test trusts the
// chart's normal orientation and this chart's normal points inward. Safe
// by the one-sided rule and useless in practice: it would stop every step.
TEST_CASE("A sweep does not report contact with a surface it is far from",
          "[physics][ccd][certified]") {
    const double r = 0.01;
    const ParametricSurface<double> tube(
        [r](double u, double v) { return V3{u, r * std::cos(v), r * std::sin(v)}; },
        {-5.0, 5.0, 0.0, 2.0 * std::numbers::pi}, false, true);
    const auto s = sweep_point_surface<double>(V3{0.37, 0.0, 3.0}, V3{0.0, 0.0, -6.0},
                                               LipschitzChart<double>{tube, 1.0});
    const double toi_true = (3.0 - r) / 6.0;
    INFO(std::format("true toi {:.6f}, reported hit={} toi={:.6f}", toi_true, s.hit, s.toi));
    REQUIRE(s.hit);
    CHECK(s.toi <= toi_true + 1e-9);
    CHECK(s.toi > toi_true - 0.01);
}

// A path that enters a unit sphere by 1e-12 -- the quadratic's
// discriminant is tiny and positive, and rounding may push it either way.
// A graze that touches must not come back as a miss.
TEST_CASE("A swept quadric counts a graze that enters as a contact",
          "[physics][ccd][certified]") {
    const auto q = geometry::Quadric<double>::sphere(1.0);
    const double x = 1.0 - 1e-12;
    const auto s = sweep_point_quadric<double>(V3{x, 0.0, 3.0}, V3{0.0, 0.0, -6.0}, q);
    INFO(std::format("hit={} toi={}", s.hit, s.toi));
    REQUIRE(s.hit);
    CHECK(s.toi <= 0.5 + 1e-6);   // the true first touch is at z = sqrt(2e-12) ~ 0, t ~ 0.5
}

// The same scene at the bottom and the top of the range a scene is built
// in: an absolute tolerance is either never reached or always met.
TEST_CASE("A sweep's answer does not depend on the units", "[physics][ccd][certified]") {
    for (double scale : {1e-6, 1.0, 1e6}) {
        const auto q = geometry::Quadric<double>::sphere(scale);
        const V3 start{0.0, 0.0, 3.0 * scale}, disp{0.0, 0.0, -6.0 * scale};
        const auto exact = sweep_point_quadric<double>(start, disp, q);
        INFO(std::format("scale {} exact hit={} toi={}", scale, exact.hit, exact.toi));
        REQUIRE(exact.hit);
        CHECK_THAT(exact.toi, Catch::Matchers::WithinAbs(1.0 / 3.0, 1e-9));

        const auto chart = sphere_chart(scale);
        const auto swept = sweep_point_surface<double>(start, disp, sphere_chart_bound(chart, scale));
        INFO(std::format("scale {} swept hit={} toi={}", scale, swept.hit, swept.toi));
        REQUIRE(swept.hit);
        CHECK(swept.toi <= 1.0 / 3.0 + 1e-9);
        CHECK(swept.toi > 1.0 / 3.0 - 1e-3);
    }
}

// Head-on into a sphere chart's pole with only the single constant: every
// cell touching the pole spans a range of u whose image is tiny, so the
// floor stays loose and the search runs into its budget. That must cost
// evaluations and an early answer, never a late one or a miss.
TEST_CASE("A chart's degenerate point costs a sweep its budget, not its answer",
          "[physics][ccd][certified]") {
    const auto chart = sphere_chart(1.0);
    const std::size_t budget = std::size_t{1} << 18;
    const auto s = sweep_sphere_surface<double>(V3{0.0, 0.0, 3.0}, 0.0, V3{0.0, 0.0, -6.0},
                                                LipschitzChart<double>{chart, 1.0}, 64, budget);
    INFO(std::format("hit={} toi={} evaluations={}", s.hit, s.toi, s.evaluations));
    REQUIRE(s.hit);
    CHECK(s.toi <= 1.0 / 3.0 + 1e-12);
    CHECK(s.evaluations <= budget + 4);
}

namespace {

// A path that grazes: through a point `depth` below (depth > 0) or above
// the surface along its normal at a random surface point, tangent to it
// there. Half the fuzz is these, since uniform random segments almost
// never come near the case that is hard.
struct Path { V3 start, disp; bool graze; };

template<typename Rng>
Path random_path(Rng& rng, double reach, const V3& foot, const V3& normal) {
    std::uniform_real_distribution<double> uni(-1.0, 1.0);
    if (std::uniform_int_distribution<int>(0, 1)(rng) == 0) {
        V3 a{uni(rng), uni(rng), uni(rng)}, b{uni(rng), uni(rng), uni(rng)};
        return {V3{a * reach}, V3{(b - a) * reach}, false};
    }
    V3 t{uni(rng), uni(rng), uni(rng)};
    t = V3{t - normal * t.dot(normal)};
    t = V3{t * (1.0 / t.norm())};
    const double depth = std::pow(10.0, -std::uniform_real_distribution<double>(2.0, 12.0)(rng)) *
                         (uni(rng) < 0 ? -1.0 : 1.0);
    const V3 through{foot - normal * depth};
    return {V3{through - t * reach}, V3{t * (2.0 * reach)}, true};
}

}  // namespace

// The chart search against the closed forms, which are an independent
// route to the same answer: a quadric's first root for the sphere and the
// exact torus distance, advanced by, for the torus. The one-sided rule is
// the whole correctness check -- a late time or a miss where the oracle
// hits fails. Tightness is asked of the straight-through paths only: at a
// graze, advancement crawls (the distance falls as x^2 along the tangent)
// and the iteration cap ends it with an early contact, which is the
// honest answer and is only counted here. Measured on the sphere at 64
// iterations: every graze early by up to 0.0096 of the step; at 1024, none
// beyond 6e-4 but 235k evaluations a query. Closing that is the search in
// time, not more iterations.
TEST_CASE("A chart sweep never answers later than the closed form",
          "[physics][ccd][certified][fuzz]") {
    std::mt19937_64 rng(20260925);
    std::uniform_real_distribution<double> ang(0.0, 2.0 * std::numbers::pi);
    constexpr int kQueries = 1500;

    SECTION("sphere") {
        const double r = 1.0;
        const auto chart = sphere_chart(r);
        const auto bound = sphere_chart_bound(chart, r);
        const auto q = geometry::Quadric<double>::sphere(r);
        int oracle_hits = 0, early = 0, grazes = 0, grazes_early = 0;
        for (int i = 0; i < kQueries; ++i) {
            const double u = ang(rng), v = std::acos(std::uniform_real_distribution<double>(-1, 1)(rng));
            const V3 n{std::sin(v) * std::cos(u), std::sin(v) * std::sin(u), std::cos(v)};
            const auto path = random_path(rng, 3.0, V3{n * r}, n);
            if (path.start.norm() <= r * (1 + 1e-6)) continue;   // a chart has no inside
            const auto exact = sweep_point_quadric<double>(path.start, path.disp, q);
            const auto swept = sweep_point_surface<double>(path.start, path.disp, bound);
            INFO(std::format("query {} exact hit={} toi={:.12f}; swept hit={} toi={:.12f}", i,
                             exact.hit, exact.toi, swept.hit, swept.toi));
            if (exact.hit) {
                REQUIRE(swept.hit);
                REQUIRE(swept.toi <= exact.toi + 1e-12);
                const bool is_early = swept.toi < exact.toi - 1e-3;
                if (path.graze) { ++grazes; grazes_early += is_early; }
                else { ++oracle_hits; early += is_early; }
            }
        }
        INFO(std::format("straight: {} hits, {} more than 1e-3 early; grazes: {} hits, {} early",
                         oracle_hits, early, grazes, grazes_early));
        CHECK(oracle_hits > kQueries / 20);
        CHECK(grazes > kQueries / 4);
        CHECK(early <= oracle_hits / 20);   // measured 2-3%: random paths that happen to graze
    }

    SECTION("torus") {
        const geometry::Torus<double> torus{V3{}, V3{0.0, 0.0, 1.0}, 1.0, 0.3};
        const double R = torus.major_radius, r = torus.minor_radius;
        const ParametricSurface<double> chart(
            [R, r](double u, double v) {
                return V3{(R + r * std::cos(v)) * std::cos(u), (R + r * std::cos(v)) * std::sin(u),
                          r * std::sin(v)};
            },
            {0.0, 2.0 * std::numbers::pi, 0.0, 2.0 * std::numbers::pi}, true, true);
        const LipschitzChart<double> bound{chart, R + r, [R, r](double u0, double u1, double v0, double v1) {
                                               return (R + r) * (u1 - u0) / 2 + r * (v1 - v0) / 2;
                                           }};
        int oracle_hits = 0, early = 0, grazes = 0, grazes_early = 0;
        for (int i = 0; i < kQueries; ++i) {
            const double u = ang(rng), v = ang(rng);
            const V3 n{std::cos(v) * std::cos(u), std::cos(v) * std::sin(u), std::sin(v)};
            const V3 foot{V3{std::cos(u), std::sin(u), 0.0} * R + n * r};
            const auto path = random_path(rng, 2.0, foot, n);
            if (point_to(path.start, torus).inside || point_to(path.start, torus).distance < 1e-6) continue;
            const auto exact = sweep_sphere_surface<double>(path.start, 0.0, path.disp, torus);
            const auto swept = sweep_point_surface<double>(path.start, path.disp, bound);
            // Both stop within sqrt(eps) * scale of the surface, so either
            // may be that much earlier than the other.
            const double slack = 2.0 * std::sqrt(std::numeric_limits<double>::epsilon()) *
                                 std::max(path.start.norm(), V3{path.start + path.disp}.norm()) /
                                 path.disp.norm();
            INFO(std::format("query {} exact hit={} toi={:.12f}; swept hit={} toi={:.12f}", i,
                             exact.hit, exact.toi, swept.hit, swept.toi));
            if (exact.hit) {
                REQUIRE(swept.hit);
                REQUIRE(swept.toi <= exact.toi + slack);
                const bool is_early = swept.toi < exact.toi - 1e-3;
                if (path.graze) { ++grazes; grazes_early += is_early; }
                else { ++oracle_hits; early += is_early; }
            }
        }
        INFO(std::format("straight: {} hits, {} more than 1e-3 early; grazes: {} hits, {} early",
                         oracle_hits, early, grazes, grazes_early));
        CHECK(oracle_hits > kQueries / 20);
        CHECK(grazes > kQueries / 4);
        CHECK(early <= oracle_hits / 20);   // measured 2-3%: random paths that happen to graze
    }
}

// A ball rather than a point: first contact when the centre is a radius
// away, the same answer from the closed-form sphere and from its chart.
TEST_CASE("A swept ball touches when its centre is a radius away", "[physics][ccd][certified]") {
    const double radius = 0.25;
    const V3 start{0.0, 0.0, 3.0}, disp{0.0, 0.0, -6.0};
    const double toi_true = (3.0 - 1.0 - radius) / 6.0;

    const auto exact = sweep_sphere_surface<double>(start, radius, disp, Sphere<2, double>{1.0});
    REQUIRE(exact.hit);
    CHECK(exact.toi <= toi_true + 1e-12);
    CHECK(exact.toi > toi_true - 1e-6);

    // Through the chart a ball costs more than a point. A point's distance
    // to the surface has a sharp minimum, of curvature 1/gap, and a
    // first-order bound needs a fixed number of cells to certify it; a
    // ball's is as flat as the two surfaces are alike, and the cells needed
    // grow as 1/gap. So the budget, not the tolerance, ends the sweep, and
    // how early it answers falls as 1/budget: measured 1.0e-4 at 2^14
    // evaluations, 2.4e-5 at 2^16, 9.2e-8 at 2^24. A bound of second order
    // is what removes that; here the budget is chosen and the answer held
    // to it.
    const auto chart = sphere_chart(1.0);
    const auto swept = sweep_sphere_surface<double>(start, radius, disp, sphere_chart_bound(chart, 1.0),
                                                    64, std::size_t{1} << 16);
    INFO(std::format("chart hit={} toi={} evaluations={}", swept.hit, swept.toi, swept.evaluations));
    REQUIRE(swept.hit);
    CHECK(swept.toi <= toi_true + 1e-12);
    CHECK(swept.toi > toi_true - 1e-4);

    // Passing a radius and a half beside it is a miss, and a proved one.
    const auto beside = sweep_sphere_surface<double>(V3{1.0 + 1.5 * radius, 0.0, 3.0}, radius, disp,
                                                     Sphere<2, double>{1.0});
    CHECK_FALSE(beside.hit);
    CHECK(beside.certified);
}
