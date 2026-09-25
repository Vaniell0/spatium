// A geodesic ray against a geodesic ball, in closed form, held to a dense
// walk of the space's own distance along the ray -- an oracle that uses
// nothing of the formula it checks, only exp_map and distance.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <spatium/spatial/ball_tree.hpp>
#include <spatium/spatial/geodesic_ball.hpp>

#include <cmath>
#include <format>
#include <random>

using namespace spatium;
using namespace spatium::spatial;

namespace {

template<std::size_t N, typename T>
Vec<T, N> gaussian(std::mt19937_64& rng) {
    std::normal_distribution<T> g(0, 1);
    Vec<T, N> v{};
    for (std::size_t i = 0; i < N; ++i) v[i] = g(rng);
    return v;
}

// A point, a unit tangent there, and the rest of what a space needs to
// generate a case, one overload set per space.
template<std::size_t N>
struct Euc {
    using S = Euclidean<N, double>;
    using P = Vec<double, N>;
    S space{};
    P point(std::mt19937_64& rng) const { return P{gaussian<N, double>(rng) * 2.0}; }
    P unit_tangent(const P&, std::mt19937_64& rng) const {
        P v = gaussian<N, double>(rng);
        return P{v * (1.0 / v.norm())};
    }
    // A unit tangent at p perpendicular to u.
    P normal_to(const P&, const P& u, std::mt19937_64& rng) const {
        P w = gaussian<N, double>(rng);
        w = P{w - u * w.dot(u)};
        return P{w * (1.0 / w.norm())};
    }
    double reach() const { return 8.0; }
};

template<std::size_t N>
struct Sph {
    using S = Sphere<N, double>;
    using P = Vec<double, N + 1>;
    S space{2.0};
    P point(std::mt19937_64& rng) const {
        P v = gaussian<N + 1, double>(rng);
        return P{v * (space.radius / v.norm())};
    }
    P unit_tangent(const P& p, std::mt19937_64& rng) const {
        P w = gaussian<N + 1, double>(rng);
        const P n{p * (1.0 / p.norm())};
        w = P{w - n * w.dot(n)};
        return P{w * (1.0 / w.norm())};
    }
    P normal_to(const P& p, const P& u, std::mt19937_64& rng) const {
        P w = unit_tangent(p, rng);
        w = P{w - u * w.dot(u)};
        return P{w * (1.0 / w.norm())};
    }
    // Three turns, so the periodic case -- a second pass of the same ball --
    // is exercised and the first entry must still be the one returned.
    double reach() const { return 3.0 * 2.0 * std::numbers::pi * space.radius; }
};

template<std::size_t N>
struct Hyp {
    using S = Hyperbolic<N, double>;
    using P = Vec<double, N + 1>;
    S space{};
    P origin() const { P o{}; o[0] = 1.0; return o; }
    P tangent_part(const P& p, P w) const {
        w = P{w + p * S::minkowski(w, p)};
        return P{w * (1.0 / std::sqrt(S::minkowski(w, w)))};
    }
    P point(std::mt19937_64& rng) const {
        const P o = origin();
        P w = gaussian<N + 1, double>(rng);
        w[0] = 0.0;
        w = P{w * (1.0 / w.norm())};
        return space.exp_map(o, w, std::uniform_real_distribution<double>(0.0, 2.5)(rng));
    }
    P unit_tangent(const P& p, std::mt19937_64& rng) const {
        return tangent_part(p, gaussian<N + 1, double>(rng));
    }
    P normal_to(const P& p, const P& u, std::mt19937_64& rng) const {
        P w = unit_tangent(p, rng);
        w = P{w - u * S::minkowski(w, u)};
        return P{w * (1.0 / std::sqrt(S::minkowski(w, w)))};
    }
    double reach() const { return 6.0; }
};

// The first t at which the ray is inside the ball, by walking it.
template<typename G>
std::pair<double, double> walk(const G& g, const typename G::P& p, const typename G::P& v,
                               const GeodesicBall<typename G::S>& b, double tmax, int steps) {
    double prev = 0.0;
    for (int k = 0; k <= steps; ++k) {
        const double t = tmax * k / steps;
        if (g.space.distance(g.space.exp_map(p, v, t), b.c) <= b.r) return {prev, t};
        prev = t;
    }
    return {-1.0, -1.0};
}

template<typename G>
void check_space(const G& g, const char* name) {
    std::mt19937_64 rng(20260925);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const int kCases = 300, kSteps = 200000;
    int entered = 0, grazes = 0;
    for (int i = 0; i < kCases; ++i) {
        const auto p = g.point(rng);
        const auto v = g.unit_tangent(p, rng);
        const double tmax = g.reach() * (0.2 + 0.8 * unit(rng));
        GeodesicBall<typename G::S> b;
        if (i % 2 == 0) {
            b = {g.point(rng), 0.05 + 1.0 * unit(rng)};
        } else {
            // Tangent to the ray at t*, off by a relative 1e-9..1e-3 either way.
            const double tstar = tmax * unit(rng), r = 0.1 + 0.8 * unit(rng);
            const auto at = g.space.exp_map(p, v, tstar);
            // the ray's direction at t*, by a tiny finite step, then a normal to it
            const auto ahead = g.space.exp_map(p, v, tstar + 1e-6);
            auto u = typename G::P{(ahead - at) * 1e6};
            u = g.unit_tangent(at, rng) * 0.0 + u;   // (u is already tangent up to 1e-6)
            const auto n = g.normal_to(at, typename G::P{u * (1.0 / std::sqrt(std::abs(u.dot(u))))}, rng);
            const double off = std::pow(10.0, -9.0 + 6.0 * unit(rng)) * (unit(rng) < 0.5 ? -1.0 : 1.0);
            b = {g.space.exp_map(at, n, r * (1.0 + off)), r};
            ++grazes;
        }
        const auto got = ray_interval(g.space, p, v, b, tmax);
        const auto [before, first] = walk(g, p, v, b, tmax, kSteps);
        INFO(std::format("{} case {}: walk first {:.9f} (after {:.9f}), closed form {}", name, i, first,
                         before, got ? std::format("[{:.9f}, {:.9f}]", got->first, got->second) : "miss"));
        if (first >= 0.0) {
            ++entered;
            REQUIRE(got);
            // The true entry lies between the last sample outside and the
            // first inside; the closed form must land there.
            CHECK(got->first <= first + 1e-9);
            CHECK(got->first >= before - 1e-9);
        } else if (got) {
            // Between two samples the ray can graze in and out; then the
            // formula's entry must really be at the ball's edge.
            const double d = g.space.distance(g.space.exp_map(p, v, got->first), b.c);
            CHECK(d <= b.r * (1.0 + 1e-6) + 1e-9);
        }
    }
    INFO(std::format("{}: {} entered, {} built as grazes", name, entered, grazes));
    CHECK(entered > kCases / 5);
}

}  // namespace

TEST_CASE("A geodesic ray meets a ball where walking the ray does", "[spatial][geodesic_ball]") {
    SECTION("Euclidean 3") { check_space(Euc<3>{}, "E3"); }
    SECTION("Sphere 2") { check_space(Sph<2>{}, "S2"); }
    SECTION("Sphere 3") { check_space(Sph<3>{}, "S3"); }
    SECTION("Hyperbolic 2") { check_space(Hyp<2>{}, "H2"); }
    SECTION("Hyperbolic 3") { check_space(Hyp<3>{}, "H3"); }
}

// A merged ball holds every point of both children, sampled, on the sphere
// -- where past the injectivity radius the merge falls back to the triangle
// inequality -- and in hyperbolic space.
template<typename G>
static void check_merge(const G& g) {
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    for (int i = 0; i < 200; ++i) {
        const GeodesicBall<typename G::S> a{g.point(rng), 0.05 + unit(rng)}, b{g.point(rng), 0.05 + unit(rng)};
        const auto m = merge(g.space, a, b);
        for (const auto& child : {a, b})
            for (int k = 0; k < 50; ++k) {
                const auto x = g.space.exp_map(child.c, g.unit_tangent(child.c, rng), child.r * unit(rng));
                INFO(std::format("case {}: d {:.12f} r {:.12f}", i, g.space.distance(m.c, x), m.r));
                REQUIRE(g.space.distance(m.c, x) <= m.r);
            }
        CHECK(lower_distance(g.space, m, a.c) == 0.0);
    }
}

TEST_CASE("A merged geodesic ball contains both of its children", "[spatial][geodesic_ball]") {
    SECTION("Sphere 2") { check_merge(Sph<2>{}); }
    SECTION("Hyperbolic 3") { check_merge(Hyp<3>{}); }
}

// The tree answers what checking every item answers: the same first item a
// ray enters, at the same parameter, and the same nearest item -- on a
// sphere, where merged nodes can pass the injectivity radius, in
// hyperbolic space, and in the flat case.
template<typename G>
static void check_tree(const G& g, const char* name) {
    std::mt19937_64 rng(11);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::vector<GeodesicBall<typename G::S>> items;
    for (int i = 0; i < 600; ++i) items.push_back({g.point(rng), 0.02 + 0.1 * unit(rng)});
    const auto tree = GeodesicBallTree<typename G::S>::build(g.space, items);
    int hits = 0;
    for (int k = 0; k < 1500; ++k) {
        const auto p = g.point(rng);
        const auto v = g.unit_tangent(p, rng);
        const double tmax = g.reach() * unit(rng);
        std::optional<std::pair<std::size_t, double>> brute;
        for (std::size_t i = 0; i < items.size(); ++i)
            if (auto in = ray_interval(g.space, p, v, items[i], tmax); in && (!brute || in->first < brute->second))
                brute = std::pair{i, in->first};
        const auto got = tree.ray_cast(p, v, tmax);
        INFO(std::format("{} ray {}", name, k));
        REQUIRE(got.has_value() == brute.has_value());
        if (got) {
            ++hits;
            CHECK(got->t == brute->second);
        }
        double best = std::numeric_limits<double>::infinity();
        for (const auto& b : items) best = std::min(best, lower_distance(g.space, b, p));
        const auto near = tree.nearest(p);
        REQUIRE(near);
        CHECK(near->distance == best);
    }
    INFO(std::format("{}: {} rays hit", name, hits));
    CHECK(hits > 100);
}

TEST_CASE("A geodesic ball tree answers as checking every ball does", "[spatial][geodesic_ball]") {
    SECTION("Euclidean 3") { check_tree(Euc<3>{}, "E3"); }
    SECTION("Sphere 2") { check_tree(Sph<2>{}, "S2"); }
    SECTION("Hyperbolic 3") { check_tree(Hyp<3>{}, "H3"); }
}
