// Minkowski spacetime in geometric units (spaces/minkowski.hpp), and what makes a
// Lorentzian space a different concept from a Riemannian one.
//
// The checks are pairs of paths that share nothing: relativistic velocity addition
// against the sum of hyperbolic distances against the composition of boost matrices;
// the Thomas-Wigner rotation of two boosts against the angular defect of a triangle
// in Hyperbolic<2> and against a closed form; the post-Newtonian expansion read from
// a Series against its known coefficients.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/core/verify.hpp>
#include <spatium/algebra/series.hpp>
#include <spatium/spaces/hyperbolic.hpp>
#include <spatium/spaces/minkowski.hpp>
#include <cmath>
#include <numbers>
#include <vector>

using namespace spatium;
using Catch::Matchers::WithinAbs;

namespace {

using M3 = Minkowski<3>;
using E4 = Vec<double, 4>;

// A "spacetime" with two time directions: signature (-, -, +). It is not Lorentzian,
// and the reverse triangle inequality must fail on it -- the control for the verifier.
struct TwoTimes {
    using ScalarType = double;
    using PointType = Vec<double, 3>;
    using TangentVector = Vec<double, 3>;
    static constexpr std::size_t dimension = 3;
    bool contains(const PointType&) const { return true; }
    double interval(const PointType& a, const PointType& b) const {
        const PointType d = b - a;
        return -d[0] * d[0] - d[1] * d[1] + d[2] * d[2];
    }
    Causal causal(const PointType& a, const PointType& b) const {
        const double s = interval(a, b);
        return s < -1e-12 ? Causal::Timelike : (s > 1e-12 ? Causal::Spacelike : Causal::Null);
    }
    double proper_time(const PointType& a, const PointType& b) const { return std::sqrt(-interval(a, b)); }
    bool precedes(const PointType& a, const PointType& b) const { return causal(a, b) != Causal::Spacelike && b[0] > a[0]; }
    PointType exp_map(const PointType& p, const TangentVector& v, double t) const { return PointType{p + v * t}; }
    TangentVector log_map(const PointType& p, const PointType& q) const { return TangentVector{q - p}; }
    double metric_at(const PointType&, const TangentVector& u, const TangentVector& v) const {
        return -u[0] * v[0] - u[1] * v[1] + u[2] * v[2];
    }
};
static_assert(LorentzianManifold<TwoTimes>);

}  // namespace

TEST_CASE("a Lorentzian space is its own concept: no distance", "[minkowski]") {
    static_assert(LorentzianManifold<Minkowski<3>>);
    static_assert(!MetricSpace<Minkowski<3>>);
    static_assert(!RiemannianManifold<Minkowski<3>>);
    static_assert(RiemannianManifold<Hyperbolic<3>> && !LorentzianManifold<Hyperbolic<3>>);
    CHECK(Minkowski<3>::dimension == 4);
}

TEST_CASE("the interval classifies events and gives proper time", "[minkowski]") {
    const M3 m;
    const E4 origin{0, 0, 0, 0};
    const E4 inside{2, 1, 0, 0}, outside{1, 2, 0, 0}, on_cone{1, 0.6, 0.8, 0};
    CHECK(m.causal(origin, inside) == Causal::Timelike);
    CHECK(m.causal(origin, outside) == Causal::Spacelike);
    CHECK(m.causal(origin, on_cone) == Causal::Null);
    CHECK_THAT(m.interval(origin, inside), WithinAbs(-3.0, 1e-14));
    CHECK_THAT(m.proper_time(origin, inside), WithinAbs(std::sqrt(3.0), 1e-14));
    CHECK_THAT(m.proper_length(origin, outside), WithinAbs(std::sqrt(3.0), 1e-14));
    CHECK(std::isnan(m.proper_time(origin, outside)));         // spacelike: no clock, NaN not a number
    CHECK(std::isnan(m.proper_length(origin, inside)));
    CHECK(m.precedes(origin, inside));
    CHECK(!m.precedes(inside, origin));
    CHECK(!m.precedes(origin, outside));
    // geodesics are straight lines
    const E4 v{1, 0.2, 0, 0};
    const E4 p = m.exp_map(origin, v, 3.0);
    CHECK_THAT(p[0], WithinAbs(3.0, 1e-14));
    CHECK_THAT(m.log_map(origin, p)[1], WithinAbs(0.6, 1e-14));
}

TEST_CASE("the reverse triangle inequality holds, and fails with two times", "[minkowski][symmetry]") {
    const std::vector<E4> events{{0, 0, 0, 0}, {1, 0.2, 0, 0}, {2.5, 0.5, 0.3, 0}, {4, 1, 0.5, 0.2},
                                 {5, 0, 0, 0}, {0.5, 3, 0, 0}, {3, -0.4, 0.6, 0.1}};
    const auto ok = verify_lorentzian(M3{}, events);
    INFO(ok.message());
    CHECK(ok.passed);

    // The twins: a straight clock reads longer than one that turned back.
    const M3 m;
    const E4 a{0, 0, 0, 0}, far{3, 2, 0, 0}, b{6, 0, 0, 0};
    CHECK(m.proper_time(a, b) > m.proper_time(a, far) + m.proper_time(far, b));

    // control: the same axiom must fail where the signature is wrong
    const std::vector<Vec<double, 3>> pts{{0, 0, 0}, {1, 0, 0}, {1.5, 1, 0}};
    const auto bad = verify_lorentzian(TwoTimes{}, pts);
    CHECK(!bad.passed);
}

TEST_CASE("a boost preserves the interval and proper time", "[minkowski]") {
    const M3 m;
    const auto L = M3::boost(Vec<double, 3>{0.3, -0.4, 0.2});
    const E4 a{0.3, 1.0, -0.5, 0.7}, b{2.0, 0.4, 0.1, -0.3};
    const E4 La = L * a, Lb = L * b;
    CHECK_THAT(m.interval(La, Lb), WithinAbs(m.interval(a, b), 1e-13));
    // Lambda^T eta Lambda = eta
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j) {
            double s = 0;
            for (std::size_t k = 0; k < 4; ++k) s += L(k, i) * (k == 0 ? -1.0 : 1.0) * L(k, j);
            CHECK_THAT(s, WithinAbs(i == j ? (i == 0 ? -1.0 : 1.0) : 0.0, 1e-13));
        }
    // the zero boost is the identity, with no 0/0 in it
    const auto I = M3::boost(Vec<double, 3>{0, 0, 0});
    CHECK(I(0, 0) == 1.0);
    CHECK(I(2, 2) == 1.0);
    CHECK(I(1, 2) == 0.0);
    // light stays light
    CHECK(std::isnan(M3::gamma(Vec<double, 3>{1.0, 0, 0})));
}

TEST_CASE("the four-velocities are Hyperbolic<N>: a rapidity is a distance", "[minkowski][symmetry]") {
    const Hyperbolic<3> h;
    const Vec<double, 3> v{0.3, -0.4, 0.2};
    const E4 u = M3::four_velocity(v);
    CHECK(h.contains(u));
    // rapidity = atanh(|v|) = the hyperbolic distance from rest
    CHECK_THAT(h.distance(Hyperbolic<3>::origin(), u), WithinAbs(std::atanh(std::sqrt(v.norm_squared())), 1e-13));
}

TEST_CASE("collinear velocities add by three independent paths", "[minkowski][symmetry]") {
    const double v1 = 0.6, v2 = 0.7;
    const Hyperbolic<1> h;
    // 1. the relativistic addition formula
    const double added = (v1 + v2) / (1 + v1 * v2);
    // 2. rapidities add: atanh v1 + atanh v2, and it is a hyperbolic distance
    const double rapidity = std::atanh(v1) + std::atanh(v2);
    // 3. the composition of the two boosts as matrices
    const auto L = Minkowski<1>::boost(Vec<double, 1>{v1}) * Minkowski<1>::boost(Vec<double, 1>{v2});
    const double from_matrices = L(1, 0) / L(0, 0);
    CHECK_THAT(std::atanh(added), WithinAbs(rapidity, 1e-13));
    CHECK_THAT(from_matrices, WithinAbs(added, 1e-13));
    // and the composite four-velocity is that far from rest in Hyperbolic<1>
    Vec<double, 2> w{L(0, 0), L(1, 0)};
    CHECK_THAT(h.distance(Hyperbolic<1>::origin(), w), WithinAbs(rapidity, 1e-12));
}

namespace {
// The angle at p between the geodesics to q and to r, in Hyperbolic<2>.
double angle_at(const Hyperbolic<2>& h, const Vec<double, 3>& p, const Vec<double, 3>& q, const Vec<double, 3>& r) {
    const auto a = h.log_map(p, q), b = h.log_map(p, r);
    const double c = Hyperbolic<2>::minkowski(a, b) /
                     std::sqrt(Hyperbolic<2>::minkowski(a, a) * Hyperbolic<2>::minkowski(b, b));
    return std::acos(std::clamp(c, -1.0, 1.0));
}
}  // namespace

TEST_CASE("the Thomas-Wigner rotation is the defect of a hyperbolic triangle", "[minkowski][symmetry]") {
    // Two boosts at right angles, u along x and then v along y. Their composite is a boost
    // and a rotation; the rotation is what the composition of boosts leaves behind.
    for (const auto& [u, v] : {std::pair{0.6, 0.5}, std::pair{0.9, 0.8}, std::pair{0.3, 0.2}}) {
        const auto Bu = Minkowski<2>::boost(Vec<double, 2>{u, 0.0});
        const auto Bv = Minkowski<2>::boost(Vec<double, 2>{0.0, v});
        const auto L = Bu * Bv;

        // path 1: strip the boost off the composite and read the rotation that is left
        const double w0 = L(0, 0);
        const Vec<double, 2> w{L(1, 0) / w0, L(2, 0) / w0};
        const auto R = Minkowski<2>::boost(Vec<double, 2>{-w[0], -w[1]}) * L;
        CHECK_THAT(R(0, 0), WithinAbs(1.0, 1e-12));                    // purely spatial
        const double angle = std::abs(std::atan2(R(2, 1), R(1, 1)));

        // path 2: the triangle (rest, u, composite) in the hyperbolic plane of velocities
        const Hyperbolic<2> h;
        const auto O = Hyperbolic<2>::origin();
        const auto U = Minkowski<2>::four_velocity(Vec<double, 2>{u, 0.0});
        const Vec<double, 3> W{L(0, 0), L(1, 0), L(2, 0)};
        const double defect = std::numbers::pi - angle_at(h, O, U, W) - angle_at(h, U, O, W) - angle_at(h, W, O, U);

        // path 3: the closed form for perpendicular boosts, cos = (g_u + g_v) / (1 + g_u g_v)
        const double gu = 1 / std::sqrt(1 - u * u), gv = 1 / std::sqrt(1 - v * v);
        const double closed = std::acos((gu + gv) / (1 + gu * gv));

        INFO("u " << u << " v " << v);
        CHECK(angle > 1e-3);                                            // there is a rotation to find
        CHECK_THAT(defect, WithinAbs(angle, 1e-10));
        CHECK_THAT(closed, WithinAbs(angle, 1e-10));
    }
}

TEST_CASE("the Newtonian limit and its post-Newtonian corrections, read from a Series in c", "[minkowski][series]") {
    // E = (gamma - 1) m c^2 as an ordinary function of c. With c = 1/eps, a Series
    // expands it in 1/c^2 exactly: 1/2 m v^2 + 3/8 m v^4/c^2 + 5/16 m v^6/c^4 + ...
    using S = Series<double, 12>;
    const double m = 1.7, v = 0.8;
    const S c = S::at_infinity();
    const S E = kinetic_energy(S(m), S(v), c);
    REQUIRE(E.ok());
    CHECK_THAT(E.coefficient(0), WithinAbs(0.5 * m * v * v, 1e-14));
    CHECK_THAT(E.coefficient(2), WithinAbs(3.0 / 8.0 * m * std::pow(v, 4), 1e-14));
    CHECK_THAT(E.coefficient(4), WithinAbs(5.0 / 16.0 * m * std::pow(v, 6), 1e-14));
    CHECK(E.coefficient(1) == 0.0);                                  // only even powers of 1/c

    // the same, as a limit
    const auto lim = limit([m, v](auto cc) { return kinetic_energy(decltype(cc)(m), decltype(cc)(v), cc); }, AtInfinity{});
    CHECK(lim.kind == LimitKind::Finite);
    CHECK_THAT(lim.value, WithinAbs(0.5 * m * v * v, 1e-14));
    CHECK(lim.exponent == 0);

    // a double at large c is a difference of two numbers that are the same
    const double e_double = kinetic_energy(m, v, 1e8);
    CHECK(std::abs(e_double - 0.5 * m * v * v) > 0.1 * 0.5 * m * v * v);
}
