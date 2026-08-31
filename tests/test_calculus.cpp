#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/algebra/calculus.hpp>
#include <spatium/spaces/hyperbolic.hpp>
#include <spatium/spaces/sphere.hpp>
#include <cmath>
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
