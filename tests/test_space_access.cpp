// Customization points (core/access.hpp): members, ADL, derivation.
//
// Two things are held here. Every space in spaces/ answers through the
// customization points exactly as through its members, and its closed-form
// distance agrees with the one derived from its own exp/log/metric -- a
// check of each space against itself that nothing asked for before. And a
// type with no members at all, which the member-based concepts reject, is a
// Riemannian space through ADL, and gets distance, geodesics and a Frechet
// mean it never wrote.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/core/access.hpp>
#include <spatium/algebra/calculus.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <spatium/spaces/hyperbolic.hpp>
#include <spatium/spaces/product.hpp>
#include <spatium/spaces/spd.hpp>
#include <spatium/spaces/sphere.hpp>
#include <cmath>
#include <numbers>
#include <vector>

using namespace spatium;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// ── A space somebody else wrote ───────────────────────────────
// A flat cylinder R x S^1 of radius r, points (z, theta). No nested types,
// no members: the author of this type never heard of spatium. It becomes a
// space through a traits specialisation and three free functions.

namespace third_party {
struct Cylinder { double radius = 1.0; };

using P = spatium::Vec<double, 2>;

inline double wrap(double a) {
    constexpr double pi = std::numbers::pi;
    a = std::fmod(a + pi, 2 * pi);
    if (a < 0) a += 2 * pi;
    return a - pi;
}
inline P exp_map(const Cylinder&, const P& p, const P& v, double t) {
    return P{p[0] + t * v[0], wrap(p[1] + t * v[1])};
}
inline P log_map(const Cylinder&, const P& p, const P& q) {
    return P{q[0] - p[0], wrap(q[1] - p[1])};
}
inline double metric_at(const Cylinder& c, const P&, const P& u, const P& v) {
    return u[0] * v[0] + c.radius * c.radius * u[1] * v[1];
}
} // namespace third_party

template<>
struct spatium::spaces::space_traits<third_party::Cylinder> {
    using point_type   = third_party::P;
    using tangent_type = third_party::P;
    using scalar_type  = double;
};

TEST_CASE("A space with no members is a space through ADL", "[access]") {
    using third_party::Cylinder;
    static_assert(spaces::Riemannian<Cylinder>);
    static_assert(spaces::Distanced<Cylinder>);        // derived, never written
    static_assert(!RiemannianManifold<Cylinder>);      // the member-based concept says no

    const Cylinder c{.radius = 2.0};
    const third_party::P p{0.0, 0.0}, q{3.0, std::numbers::pi / 2};
    // sqrt(dz^2 + r^2 dtheta^2)
    CHECK_THAT(spaces::distance(c, p, q),
               WithinRel(std::sqrt(9.0 + std::numbers::pi * std::numbers::pi), 1e-14));

    // Across the seam: theta = 3 and -3 are 2pi - 6 apart, not 6.
    const third_party::P a{0.0, 3.0}, b{0.0, -3.0};
    CHECK_THAT(spaces::distance(c, a, b), WithinRel(2.0 * (2 * std::numbers::pi - 6.0), 1e-12));

    const auto m = spaces::midpoint(c, a, b);
    CHECK_THAT(spaces::distance(c, a, m), WithinRel(spaces::distance(c, m, b), 1e-12));
    CHECK_THAT(std::abs(m[1]), WithinAbs(std::numbers::pi, 1e-12));   // on the seam

    // A Frechet mean the type never implemented, across the seam as well.
    const std::vector<third_party::P> pts{{1.0, 3.0}, {1.0, -3.0}, {2.0, 3.1}, {0.0, -3.1}};
    const auto mean = algebra::frechet_mean(c, pts);
    CHECK_THAT(mean[0], WithinAbs(1.0, 1e-9));
    CHECK_THAT(std::abs(mean[1]), WithinAbs(std::numbers::pi, 1e-9));
}

// ── Every space in spaces/, against itself ─────────────────────

namespace {
// Distinct pairs only. At p == q the closed forms are the weaker side:
// measured when this test was written, Sphere's acos and Hyperbolic's acosh
// give d(p, p) of 1.5e-8 to 2.1e-8 (the derivative of acos is infinite at
// 1), while |log_p p|_g is exactly 0. That is a precision limit of those
// closed forms, not a disagreement about the geometry.
template<class S>
void check_against_itself(const S& s, const std::vector<spaces::point_t<S>>& pts, double tol) {
    static_assert(spaces::Riemannian<S>);
    for (std::size_t i = 0; i < pts.size(); ++i)
        for (std::size_t j = i + 1; j < pts.size(); ++j) {
            const auto& p = pts[i];
            const auto& q = pts[j];
            // The customization point finds the member.
            CHECK(spaces::distance(s, p, q) == s.distance(p, q));
            // The closed form agrees with |log_p q|_g, derived from the
            // space's own log and metric.
            const auto v = spaces::tangent_t<S>{spaces::log_map(s, p, q)};
            CHECK_THAT(spaces::norm_at(s, p, v), WithinAbs(s.distance(p, q), tol));
            // The midpoint is equidistant, at half the distance.
            const auto m = spaces::midpoint(s, p, q);
            CHECK_THAT(s.distance(p, m), WithinAbs(s.distance(p, q) / 2, tol));
            CHECK_THAT(s.distance(m, q), WithinAbs(s.distance(p, q) / 2, tol));
        }
}
} // namespace

TEST_CASE("Closed-form distances agree with |log|_g", "[access]") {
    SECTION("Euclidean") {
        Euclidean<3> e;
        check_against_itself(e, {Vec3{0, 0, 0}, Vec3{1, 2, 3}, Vec3{-2, 0.5, 1}}, 1e-12);
    }
    SECTION("Sphere") {
        Sphere<2> s{.radius = 2.0};
        const Vec3 n{0, 0, 2};
        check_against_itself(s, {n, s.exp_map(n, Vec3{0.7, -0.3, 0}, 1.0),
                                 s.exp_map(n, Vec3{-0.4, 1.1, 0}, 1.0)}, 1e-10);
    }
    SECTION("Hyperbolic") {
        Hyperbolic<2> h;
        const auto o = Hyperbolic<2>::origin();
        check_against_itself(h, {o, h.exp_map(o, Vec3{0, 0.8, -0.3}, 1.0),
                                 h.exp_map(o, Vec3{0, -0.2, 1.1}, 1.0)}, 1e-9);
    }
    SECTION("Product S2 x E1") {
        ProductSpace<Sphere<2>, Euclidean<1>> ps;
        const Vec<double, 4> b{0, 0, 1, 0};
        check_against_itself(ps, {b, ps.exp_map(b, Vec<double, 4>{0.5, 0.1, 0, 2.0}, 1.0),
                                  ps.exp_map(b, Vec<double, 4>{-0.3, 0.6, 0, -1.0}, 1.0)}, 1e-10);
    }
    SECTION("SPD, affine-invariant") {
        SPDAffineInvariant<2> spd;
        using M = SPDAffineInvariant<2>::PointType;
        M a = M::identity(), b = M::identity(), c = M::identity();
        b(0, 0) = 2.0; b(0, 1) = b(1, 0) = 0.3; b(1, 1) = 0.7;
        c(0, 0) = 0.5; c(0, 1) = c(1, 0) = -0.2; c(1, 1) = 1.8;
        check_against_itself(spd, {a, b, c}, 1e-9);
    }
}
