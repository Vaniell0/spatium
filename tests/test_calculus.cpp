#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/calculus.hpp>
#include <spatium/spaces/constant_curvature.hpp>
#include <spatium/spaces/hyperbolic.hpp>
#include <spatium/spaces/sphere.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numbers>

using namespace spatium;
using Catch::Matchers::WithinAbs;

TEST_CASE("gradient of a paraboloid matches the closed form", "[calculus]") {
    // f(p) = p.x^2 + p.y^2 -> grad = (2x, 2y)
    auto f = [](auto p) { return p[0] * p[0] + p[1] * p[1]; };
    auto g = gradient(f, Vec<double, 2>{1.0, 2.0});
    CHECK_THAT(g[0], WithinAbs(2.0, 1e-9));
    CHECK_THAT(g[1], WithinAbs(4.0, 1e-9));
}

TEST_CASE("gradient of an implicit sphere gives the outward surface normal", "[calculus]") {
    // F(p) = |p|^2 - r^2 = 0 defines a sphere; grad(F) at a point on the
    // surface points radially outward. Point (3,4,0), r=5 -> normal (0.6,0.8,0).
    auto sphere = [](auto p) { return p[0] * p[0] + p[1] * p[1] + p[2] * p[2] - 25.0; };
    auto g = gradient(sphere, Vec<double, 3>{3.0, 4.0, 0.0});
    auto n = g.normalized();
    CHECK_THAT(n[0], WithinAbs(0.6, 1e-9));
    CHECK_THAT(n[1], WithinAbs(0.8, 1e-9));
    CHECK_THAT(n[2], WithinAbs(0.0, 1e-9));
}

TEST_CASE("integrate recovers known closed-form definite integrals", "[calculus]") {
    // integral of x^2 dx from 0 to 1 = 1/3
    CHECK_THAT((integrate<double>([](double x) { return x * x; }, 0.0, 1.0)),
               WithinAbs(1.0 / 3.0, 1e-9));

    // integral of sin(x) dx from 0 to pi = 2
    CHECK_THAT((integrate<double>([](double x) { return std::sin(x); }, 0.0, std::numbers::pi)),
               WithinAbs(2.0, 1e-9));

    // integral of exp(x) dx from 0 to 1 = e - 1
    CHECK_THAT((integrate<double>([](double x) { return std::exp(x); }, 0.0, 1.0)),
               WithinAbs(std::exp(1.0) - 1.0, 1e-9));
}

TEST_CASE("Function concept accepts plain lambdas with no wrapper", "[calculus]") {
    auto f = [](double x) { return x; };
    static_assert(Function<decltype(f), double, double>);
}

TEST_CASE("minimize finds the minimum of a simple bowl", "[calculus]") {
    // f(x,y) = (x-3)^2 + (y+2)^2 -> minimum at (3,-2), value 0
    auto f = [](auto p) { return (p[0] - 3.0) * (p[0] - 3.0) + (p[1] + 2.0) * (p[1] + 2.0); };
    auto theta = minimize(f, Vec<double, 2>{0.0, 0.0});
    CHECK_THAT(theta[0], WithinAbs(3.0, 1e-4));
    CHECK_THAT(theta[1], WithinAbs(-2.0, 1e-4));
}

TEST_CASE("minimize handles an anisotropic (ill-conditioned) bowl", "[calculus]") {
    // f(x,y) = (x-1)^2 + 10*(y-2)^2 -> minimum at (1,2)
    // Steeper in y than x -- exercises the Armijo line search, not just a
    // fixed step size.
    auto f = [](auto p) {
        return (p[0] - 1.0) * (p[0] - 1.0) + 10.0 * (p[1] - 2.0) * (p[1] - 2.0);
    };
    auto theta = minimize(f, Vec<double, 2>{-5.0, 8.0});
    CHECK_THAT(theta[0], WithinAbs(1.0, 1e-3));
    CHECK_THAT(theta[1], WithinAbs(2.0, 1e-3));
}

TEST_CASE("minimize calibrates a toy contact-style coefficient", "[calculus]") {
    // Stand-in for the real XPBD case: a "penetration" loss shaped like a
    // 1-parameter compliance search -- min of (compliance - target)^2.
    auto loss = [](auto p) {
        auto compliance = p[0];
        auto target = 0.42;
        return (compliance - target) * (compliance - target);
    };
    auto theta = minimize(loss, Vec<double, 1>{0.0});
    CHECK_THAT(theta[0], WithinAbs(0.42, 1e-4));
}

TEST_CASE("project_tangent removes the normal component under the space's own metric", "[calculus]") {
    Sphere<2> sph;
    Vec<double, 3> p{1.0, 0.0, 0.0};
    Vec<double, 3> v{0.5, 0.5, 0.5};
    auto t = project_tangent(sph, p, v);
    auto n = sph.normal(p);
    CHECK_THAT(sph.metric_at(p, n, t), WithinAbs(0.0, 1e-9));

    Hyperbolic<2> hyp;
    auto hp = Hyperbolic<2>::origin();
    Vec<double, 3> hv{0.3, 0.7, -0.2};
    auto ht = project_tangent(hyp, hp, hv);
    auto hn = hyp.normal(hp);
    CHECK_THAT(hyp.metric_at(hp, hn, ht), WithinAbs(0.0, 1e-9));
}

TEST_CASE("raise_gradient is identity on a Euclidean ambient metric, sign-flip on Minkowski", "[calculus]") {
    Sphere<2> sph;
    Vec<double, 3> p{1.0, 0.0, 0.0};
    Vec<double, 3> ambient_grad{0.4, -0.2, 0.9};
    auto raised = raise_gradient(sph, p, ambient_grad);
    CHECK_THAT(raised[0], WithinAbs(ambient_grad[0], 1e-9));
    CHECK_THAT(raised[1], WithinAbs(ambient_grad[1], 1e-9));
    CHECK_THAT(raised[2], WithinAbs(ambient_grad[2], 1e-9));

    Hyperbolic<2> hyp;
    auto hp = Hyperbolic<2>::origin();
    Vec<double, 3> hgrad{1.0, 0.0, 0.0}; // d/dp of f(p) = p[0]
    auto hraised = raise_gradient(hyp, hp, hgrad);
    CHECK_THAT(hraised[0], WithinAbs(-1.0, 1e-9)); // Minkowski flips the time component
    CHECK_THAT(hraised[1], WithinAbs(0.0, 1e-9));
    CHECK_THAT(hraised[2], WithinAbs(0.0, 1e-9));
}

TEST_CASE("riemannian_minimize finds the closest point on a sphere to a target", "[calculus]") {
    // f(p) = -dot(p, target) is minimized on the sphere exactly at p = target
    // (the antipode is the max, not a competing minimum).
    Sphere<2> sph;
    Vec<double, 3> target{0.0, 0.0, 1.0};
    auto f = [target](auto p) {
        return -(p[0] * target[0] + p[1] * target[1] + p[2] * target[2]);
    };
    Vec<double, 3> start{1.0, 0.0, 0.0};
    auto result = riemannian_minimize(sph, f, start);

    CHECK_THAT(result[0], WithinAbs(0.0, 1e-3));
    CHECK_THAT(result[1], WithinAbs(0.0, 1e-3));
    CHECK_THAT(result[2], WithinAbs(1.0, 1e-3));
    CHECK_THAT(result.norm(), WithinAbs(1.0, 1e-6)); // stays on the manifold
}

TEST_CASE("frechet_mean of two sphere points is the normalized geodesic midpoint", "[calculus]") {
    // Independent ground truth, not self-consistency: for two points on a
    // unit sphere the minor-arc midpoint is exactly the normalized vector
    // sum (p+q)/|p+q| -- a separate closed form from log_map/exp_map.
    Sphere<2> sph;
    Vec<double, 3> p{1.0, 0.0, 0.0};
    Vec<double, 3> q{0.0, 1.0, 0.0};
    std::vector<Vec<double, 3>> points{p, q};

    auto mean = frechet_mean(sph, points);
    auto expected = (p + q).normalized();

    CHECK_THAT(mean[0], WithinAbs(expected[0], 1e-6));
    CHECK_THAT(mean[1], WithinAbs(expected[1], 1e-6));
    CHECK_THAT(mean[2], WithinAbs(expected[2], 1e-6));
    CHECK_THAT(mean.norm(), WithinAbs(1.0, 1e-9)); // stays on the manifold
}

TEST_CASE("frechet_mean of two hyperboloid points is the geodesic midpoint", "[calculus]") {
    // Same "two-point mean = midpoint of the connecting geodesic" theorem
    // holds on any Riemannian manifold -- cross-checked here against the
    // hyperboloid's own exp_map/log_map directly, proving frechet_mean()
    // is genuinely generic (Part of this PR's point: one implementation,
    // not one per space).
    Hyperbolic<2> hyp;
    auto p = Hyperbolic<2>::origin();
    Vec<double, 3> q{std::cosh(1.0), std::sinh(1.0), 0.0};
    std::vector<Vec<double, 3>> points{p, q};

    auto mean = frechet_mean(hyp, points);
    auto expected_midpoint = hyp.exp_map(p, hyp.log_map(p, q), 0.5);

    CHECK_THAT(mean[0], WithinAbs(expected_midpoint[0], 1e-6));
    CHECK_THAT(mean[1], WithinAbs(expected_midpoint[1], 1e-6));
    CHECK_THAT(mean[2], WithinAbs(expected_midpoint[2], 1e-6));
    CHECK_THAT(hyp.distance(mean, p), WithinAbs(hyp.distance(mean, q), 1e-6));
}

TEST_CASE("frechet_mean minimizes total squared geodesic distance", "[calculus]") {
    // Generic local-optimality sanity check, valid on any RiemannianManifold:
    // the mean's total loss must not exceed the loss centered at any single
    // sample (each sample has zero distance to itself but nonzero distance
    // to the other three, unless the samples happen to coincide).
    Sphere<2> sph;
    std::vector<Vec<double, 3>> points{
        Vec<double, 3>{1.0, 0.0, 0.0}.normalized(),
        Vec<double, 3>{0.3, 0.9, 0.2}.normalized(),
        Vec<double, 3>{-0.5, 0.1, 0.8}.normalized(),
        Vec<double, 3>{0.1, -0.6, 0.7}.normalized(),
    };

    auto total_loss = [&](const Vec<double, 3>& center) {
        double sum = 0.0;
        for (const auto& p : points) {
            auto d = sph.distance(center, p);
            sum += d * d;
        }
        return sum;
    };

    auto mean = frechet_mean(sph, points);
    double mean_loss = total_loss(mean);
    for (const auto& p : points)
        CHECK(mean_loss <= total_loss(p) + 1e-9);
}

TEST_CASE("riemannian_minimize finds the closest point on a hyperboloid to a target", "[calculus]") {
    // -minkowski(p, target) == cosh(distance(p, target)) on the hyperboloid,
    // minimized exactly at p = target.
    Hyperbolic<2> hyp;
    Vec<double, 3> target = Hyperbolic<2>::origin();
    auto f = [target](auto p) {
        auto mink = -p[0] * target[0] + p[1] * target[1] + p[2] * target[2];
        return -mink;
    };
    using std::cosh; using std::sinh;
    Vec<double, 3> start{cosh(1.0), sinh(1.0), 0.0};
    auto result = riemannian_minimize(hyp, f, start);

    CHECK_THAT(result[0], WithinAbs(1.0, 1e-3));
    CHECK_THAT(result[1], WithinAbs(0.0, 1e-3));
    CHECK_THAT(result[2], WithinAbs(0.0, 1e-3));
}

// ── The answer carries what the rule knows about it ──────────

TEST_CASE("integrate_with_error returns an estimate that covers the true error, and its cost", "[calculus]") {
    struct Case { std::function<double(double)> f; double a, b, truth; };
    const Case cases[] = {
        {[](double x) { return x * x; }, 0.0, 1.0, 1.0 / 3.0},
        {[](double x) { return std::sin(x); }, 0.0, std::numbers::pi, 2.0},
        {[](double x) { return std::exp(x); }, 0.0, 1.0, std::numbers::e - 1.0},
        {[](double x) { return 1.0 / (1.0 + x * x); }, 0.0, 1.0, std::numbers::pi / 4.0},
    };
    for (const auto& c : cases) {
        long calls = 0;
        const auto f = [&](double x) { ++calls; return c.f(x); };
        const auto r = integrate_with_error<double>(f, c.a, c.b);
        CHECK(r.trusted());
        CHECK(r.evaluations == calls);                       // the cost it reports is the cost it paid
        CHECK(std::abs(r.value - c.truth) <= std::max(r.error_estimate, 1e-14));
        // and the value alone is the same number
        CHECK(integrate<double>(c.f, c.a, c.b) == r.value);
    }
}

TEST_CASE("A bound or a value that is not finite stops the rule at once, with no answer", "[calculus][infinity]") {
    // Before, a NaN walked the recursion to its depth limit of 20 -- up to 2^20
    // evaluations, minutes over Real50 -- and came back NaN anyway.
    const double inf = std::numeric_limits<double>::infinity();
    {
        const auto r = integrate_with_error<double>([](double x) { return std::exp(-x); }, 0.0, inf);
        CHECK(r.status == IntegralStatus::Failed);
        CHECK(std::isnan(r.value));
        CHECK(r.evaluations == 0);                           // the integrand was never asked
    }
    {   // an integrable singularity at an endpoint: f(0) = infinity
        const auto r = integrate_with_error<double>([](double x) { return 1.0 / std::sqrt(x); }, 0.0, 1.0);
        CHECK(r.status == IntegralStatus::Failed);
        CHECK(std::isnan(r.value));
        CHECK(r.evaluations == 3);
    }
    {   // a NaN on the first refinement point, 0.25: found there, after 5 evaluations
        const auto f = [](double x) { return (x > 0.2 && x < 0.3) ? std::numeric_limits<double>::quiet_NaN() : 1.0; };
        const auto r = integrate_with_error<double>(f, 0.0, 1.0);
        CHECK(r.status == IntegralStatus::Failed);
        CHECK(std::isnan(r.value));
        CHECK(r.evaluations == 5);
    }
    {   // the limit, stated: an adaptive rule sees only where it samples. A NaN in
        // (0.4, 0.45) lies between the nodes of a constant integrand, which
        // converges at once -- the answer is 1 and the status is Converged.
        const auto f = [](double x) { return (x > 0.4 && x < 0.45) ? std::numeric_limits<double>::quiet_NaN() : 1.0; };
        const auto r = integrate_with_error<double>(f, 0.0, 1.0);
        CHECK(r.trusted());
        CHECK_THAT(r.value, WithinAbs(1.0, 1e-14));
    }
}

TEST_CASE("A depth cap that decided the answer is reported, not hidden", "[calculus]") {
    // sin(1/(x + 1e-6)) oscillates faster than 20 halvings resolve at this
    // tolerance: panels are accepted because the depth ran out.
    const auto r = integrate_with_error<double>([](double x) { return std::sin(1.0 / (x + 1e-6)); }, 0.0, 1.0, 1e-14);
    CHECK(r.status == IntegralStatus::DepthCapped);
    CHECK(std::isfinite(r.value));
    CHECK(r.error_estimate > 1e-14);
}

TEST_CASE("The rate of convergence flags what the rule cannot resolve, and the limit of that is stated", "[calculus]") {
    // sin(20x) on [0, 1]: five samples at the root cannot resolve twenty
    // radians, so the estimate falls between the first two levels by far less
    // than the factor 8 of a smooth integrand, though every panel goes on to
    // meet its tolerance. Suspicious -- here the true error turns out small, but
    // the rule could not know that, and says so.
    const auto wiggly = integrate_with_error<double>([](double x) { return std::sin(20.0 * x); }, 0.0, 1.0, 1e-6);
    CHECK(wiggly.status == IntegralStatus::Suspicious);
    CHECK_FALSE(wiggly.trusted());

    // A smooth integrand at the same tolerance is Converged.
    const auto smooth = integrate_with_error<double>([](double x) { return std::sin(x); }, 0.0, std::numbers::pi, 1e-6);
    CHECK(smooth.status == IntegralStatus::Converged);

    // sqrt(x) at the default tolerance runs into the depth limit; the estimate
    // is large and says so: DepthCapped is the honest case.
    CHECK(integrate_with_error<double>([](double x) { return std::sqrt(x); }, 0.0, 1.0).status == IntegralStatus::DepthCapped);

    // The limit. x^0.7 at a loose tolerance is accepted as Converged, and its
    // true error is several times its estimate: one rule's estimate can be
    // small and wrong, and no status of that rule sees it. It takes a second
    // witness that samples differently (Gauss-Kronrod, tanh-sinh), and a
    // disagreement between witnesses is the signal.
    const auto rough = integrate_with_error<double>([](double x) { return std::pow(x, 0.7); }, 0.0, 1.0, 1e-3);
    CHECK(rough.status == IntegralStatus::Converged);
    CHECK(std::abs(rough.value - 1.0 / 1.7) > 3.0 * rough.error_estimate);

    // A polynomial Simpson integrates exactly has an estimate of zero at every
    // level: nothing to fall, nothing to doubt.
    const auto q = integrate_with_error<double>([](double x) { return x * x * x; }, 0.0, 2.0);
    CHECK(q.trusted());
    CHECK_THAT(q.value, WithinAbs(4.0, 1e-13));
}

TEST_CASE("The estimate and the derivative both pass through a Dual", "[calculus][dual]") {
    using D = Dual<double>;
    const D a = D::variable(2.0);
    const auto r = integrate_with_error<D>([](D x) { return x * x; }, D(0.0), a);
    CHECK(r.trusted());
    CHECK_THAT(r.value.value, WithinAbs(8.0 / 3.0, 1e-12));
    CHECK_THAT(r.value.deriv, WithinAbs(4.0, 1e-9));         // d/da of the integral is f(a) = a^2
}

TEST_CASE("riemannian_minimize on a space with intrinsic coordinates agrees with the Frechet mean", "[calculus][symmetry]") {
    // A space whose points are coordinates of the manifold itself (the kappa family, a chart)
    // has no normal: the raised gradient is already tangent, and there is nothing to project.
    // Two paths to the minimiser of the summed squared distances that share nothing: the
    // gradient of that sum through a Dual, and the fixed point of the mean of the logs.
    for (const double kappa : {-1.0, 0.0, 1.0}) {
        const ConstantCurvature<2> space{kappa};
        const std::vector<Vec<double, 2>> pts{{0.10, 0.20}, {-0.30, 0.05}, {0.25, -0.35}};

        const auto sum_of_squares = [kappa, &pts](const auto& x) {
            using Sc = std::decay_t<decltype(x[0])>;
            const ConstantCurvature<2, Sc> sp{Sc(kappa)};
            Sc acc{0};
            for (const auto& p : pts) {
                const auto d = sp.distance(x, Vec<Sc, 2>{Sc(p[0]), Sc(p[1])});
                acc = acc + d * d;
            }
            return acc;
        };
        const auto by_gradient = riemannian_minimize(space, sum_of_squares, Vec<double, 2>{0.0, 0.0});
        const auto by_logs = frechet_mean(space, pts, Vec<double, 2>{0.0, 0.0});
        INFO("kappa " << kappa);
        CHECK_THAT(by_gradient[0], WithinAbs(by_logs[0], 1e-6));
        CHECK_THAT(by_gradient[1], WithinAbs(by_logs[1], 1e-6));
    }
}
