// A metric as ten fields: the Boyer-Lindquist Kerr metric written that way
// must be KerrMetric entry for entry, derivatives included, and the
// Kerr-Schild form -- which exists nowhere else in the library to compare
// with -- is held to what any Kerr metric must satisfy: it solves the
// vacuum equations, its determinant is that of Minkowski space, and on the
// equator its g_tt is the Boyer-Lindquist one.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <spatium/physics/relativity/geodesic.hpp>
#include <spatium/physics/relativity/kerr.hpp>
#include <spatium/physics/relativity/metric_field.hpp>

#include <cmath>
#include <format>
#include <random>

using namespace spatium;
using namespace spatium::physics::relativity;
using Catch::Matchers::WithinAbs;
using V4 = Vec<double, 4>;

namespace {

// R_{mu nu} from the Christoffel symbols, their derivatives by central
// differences of step h:
//   R_{mu nu} = d_l G^l_{mu nu} - d_nu G^l_{mu l} + G^l_{l s} G^s_{mu nu} - G^l_{nu s} G^s_{mu l}
template<typename Metric>
Matrix<double, 4, 4> ricci(const Metric& m, const V4& x, double h) {
    const auto G = christoffel(m, x);
    std::array<std::array<Matrix<double, 4, 4>, 4>, 4> dG{};   // dG[k][l](mu,nu) = d_k G^l_{mu nu}
    for (std::size_t k = 0; k < 4; ++k) {
        V4 xp = x, xm = x;
        xp[k] += h;
        xm[k] -= h;
        const auto Gp = christoffel(m, xp), Gm = christoffel(m, xm);
        for (std::size_t l = 0; l < 4; ++l)
            for (std::size_t mu = 0; mu < 4; ++mu)
                for (std::size_t nu = 0; nu < 4; ++nu)
                    dG[k][l](mu, nu) = (Gp[l](mu, nu) - Gm[l](mu, nu)) / (2 * h);
    }
    Matrix<double, 4, 4> R{};
    for (std::size_t mu = 0; mu < 4; ++mu)
        for (std::size_t nu = 0; nu < 4; ++nu) {
            double v = 0;
            for (std::size_t l = 0; l < 4; ++l) {
                v += dG[l][l](mu, nu) - dG[nu][l](mu, l);
                for (std::size_t s = 0; s < 4; ++s)
                    v += G[l](l, s) * G[s](mu, nu) - G[l](nu, s) * G[s](mu, l);
            }
            R(mu, nu) = v;
        }
    return R;
}

double max_abs(const Matrix<double, 4, 4>& m) {
    double v = 0;
    for (std::size_t i = 0; i < 4; ++i)
        for (std::size_t j = 0; j < 4; ++j) v = std::max(v, std::abs(m(i, j)));
    return v;
}

double det4(const Matrix<double, 4, 4>& m) {
    // Laplace expansion along the first row; small and exact enough here.
    auto minor3 = [&](int skip) {
        int c[3], k = 0;
        for (int j = 0; j < 4; ++j) if (j != skip) c[k++] = j;
        auto e = [&](int i, int j) { return m(i + 1, c[j]); };
        return e(0, 0) * (e(1, 1) * e(2, 2) - e(1, 2) * e(2, 1)) -
               e(0, 1) * (e(1, 0) * e(2, 2) - e(1, 2) * e(2, 0)) +
               e(0, 2) * (e(1, 0) * e(2, 1) - e(1, 1) * e(2, 0));
    };
    double d = 0;
    for (int j = 0; j < 4; ++j) d += (j % 2 ? -1 : 1) * m(0, j) * minor3(j);
    return d;
}

}  // namespace

TEST_CASE("The Kerr metric written as fields is KerrMetric", "[relativity][metric_field]") {
    const double M = 1.0, a = 0.9;
    const auto field = kerr_boyer_lindquist(M, a);
    const KerrMetric<double> kerr{M, a};
    std::mt19937_64 rng(5);
    std::uniform_real_distribution<double> r(2.2, 30.0), th(0.1, 3.0), ph(0.0, 6.28);
    for (int i = 0; i < 200; ++i) {
        const V4 x{0.3, r(rng), th(rng), ph(rng)};
        const auto a_md = metric_derivatives(field, x), b_md = metric_derivatives(kerr, x);
        for (std::size_t mu = 0; mu < 4; ++mu)
            for (std::size_t nu = 0; nu < 4; ++nu) {
                CHECK_THAT(a_md.g(mu, nu), WithinAbs(b_md.g(mu, nu), 1e-13 * (1 + std::abs(b_md.g(mu, nu)))));
                for (std::size_t k = 0; k < 4; ++k)
                    CHECK_THAT(a_md.dg[k](mu, nu),
                               WithinAbs(b_md.dg[k](mu, nu), 1e-12 * (1 + std::abs(b_md.dg[k](mu, nu)))));
            }
    }
}

TEST_CASE("A metric that reads what a spacetime point has not is refused", "[relativity][metric_field]") {
    using F = io::build::Field<double>;
    std::array<F, 10> g;
    for (auto& e : g) e = F(1.0);
    g[4] = F::coord(1) + F::u();   // g_rr reading a surface parameter
    const auto m = MetricField<double>::make(g);
    REQUIRE_FALSE(m);
    CHECK(m.error().message.find("op 1 is U") != std::string::npos);
}

TEST_CASE("Both forms of the Kerr metric solve the vacuum equations", "[relativity][metric_field]") {
    const double M = 1.0, a = 0.7;
    const auto bl = kerr_boyer_lindquist(M, a);
    const auto ks = kerr_schild(M, a);
    std::mt19937_64 rng(9);
    std::uniform_real_distribution<double> u(-1.0, 1.0);

    // Off the horizon and away from the ring singularity, in both charts.
    double worst_bl = 0, worst_ks = 0;
    for (int i = 0; i < 40; ++i) {
        const double rr = 3.0 + 6.0 * (u(rng) + 1.0) / 2.0;
        const double th = 0.3 + 2.5 * (u(rng) + 1.0) / 2.0, ph = 3.0 * u(rng);
        const V4 xbl{0.0, rr, th, ph};
        const V4 xks{0.0, rr * std::sin(th) * std::cos(ph), rr * std::sin(th) * std::sin(ph), rr * std::cos(th)};
        INFO(std::format("case {}", i));
        const double rbl = max_abs(ricci(bl, xbl, 1e-4)), rks = max_abs(ricci(ks, xks, 1e-4));
        worst_bl = std::max(worst_bl, rbl);
        worst_ks = std::max(worst_ks, rks);
        // Central differences of step 1e-4 leave O(h^2) of truncation;
        // measured, the worst over these points is 1.2e-6, against the
        // control below at more than 1e-2.
        CHECK(rbl < 1e-5);
        CHECK(rks < 1e-5);
        // Kerr-Schild keeps Minkowski's determinant, -1, everywhere.
        CHECK_THAT(det4(ks(xks)), WithinAbs(-1.0, 1e-12));
    }

    // The control: the same check must fail on a metric that is not a
    // vacuum solution, or passing it says nothing.
    using F = io::build::Field<double>;
    std::array<F, 10> g;
    for (auto& e : g) e = F(0.0);
    // eta + a bump in g_tt: a static metric with a source.
    const F x = F::coord(1), y = F::coord(2), z = F::coord(3);
    g[0] = F(-1.0) + F(0.3) / (F(1.0) + x * x + y * y + z * z);
    g[4] = F(1.0);
    g[7] = F(1.0);
    g[9] = F(1.0);
    const auto bumped = *MetricField<double>::make(g);
    const double control = max_abs(ricci(bumped, V4{0.0, 1.0, 0.5, 0.2}, 1e-4));
    UNSCOPED_INFO(std::format("vacuum residual worst: BL {:.2e}, KS {:.2e}; control {:.2e}", worst_bl, worst_ks, control));
    CHECK(control > 1e-2);
    CHECK(worst_bl < 1e-5);
}

TEST_CASE("Kerr-Schild and Boyer-Lindquist agree on the equator", "[relativity][metric_field]") {
    // t is the same Killing time in both, so g_tt at the same event is the
    // same number; on the equator the Kerr-Schild radius is sqrt(r^2 + a^2).
    const double M = 1.0, a = 0.9;
    const auto bl = kerr_boyer_lindquist(M, a);
    const auto ks = kerr_schild(M, a);
    for (double r : {2.5, 4.0, 10.0, 100.0}) {
        const double rho = std::sqrt(r * r + a * a);
        CHECK_THAT(ks(V4{0.0, rho, 0.0, 0.0})(0, 0), WithinAbs(bl(V4{0.0, r, std::numbers::pi / 2, 0.0})(0, 0), 1e-13));
    }
}

// ── One clock ────────────────────────────────────────────────────
//
// A motion reads time as its first parameter, a metric as coordinate 0;
// `in_spacetime()` is the same field on the metric's clock. It must mean
// the same thing on both -- through the field and through its lowered
// form -- and it is what lets a hole's path, written as a motion, enter
// the metric that hole makes.
TEST_CASE("A motion field reads the same on the metric's clock", "[relativity][metric_field][clock]") {
    using F = io::build::Field<double>;
    const F t = F::t();
    const F e = min(max(t * F(0.5), F(0.0)), F(1.0));
    const F motion = F(3.0) * cos(F(0.7) * t) + e * e * (F(3.0) - F(2.0) * e) + sqrt(t * t + F(1.0));
    const auto clock = motion.in_spacetime();
    REQUIRE(clock);
    const auto pod = io::build::lower(motion), pod_clock = io::build::lower(*clock);
    REQUIRE(pod);
    REQUIRE(pod_clock);
    std::vector<double> s0(pod->ops.size()), s1(pod_clock->ops.size());
    for (double time : {-2.0, 0.0, 0.37, 1.5, 4.0, 11.0}) {
        io::build::FieldInputs<double> as_motion{};
        as_motion.u = time;
        io::build::FieldInputs<double> as_event{};
        as_event.x = V4{time, 5.0, -2.0, 7.0};   // the spatial coordinates must not matter
        const double a = io::build::interpret(pod->ops.data(), static_cast<std::uint32_t>(pod->ops.size()),
                                              pod->tables.data(), pod->points.data(), s0.data(), as_motion);
        const double b = io::build::interpret(pod_clock->ops.data(), static_cast<std::uint32_t>(pod_clock->ops.size()),
                                              pod_clock->tables.data(), pod_clock->points.data(), s1.data(), as_event);
        CHECK(a == b);
    }
}

TEST_CASE("A field that reads what an event has not cannot join the metric's clock",
          "[relativity][metric_field][clock]") {
    using F = io::build::Field<double>;
    const auto surface = (F::t() + F::v()).in_spacetime();
    REQUIRE_FALSE(surface);
    CHECK(surface.error().message.find("is V") != std::string::npos);
    CHECK_FALSE((F::t() * F::point(0)).in_spacetime());
    CHECK_FALSE(F([](double u, double) { return u; }).in_spacetime());
}

// A non-spinning hole on a circular path, its path written as a motion:
// at every event its metric is the static hole's metric at the offset from
// where the hole is at that event's time.
TEST_CASE("A hole whose path is a motion field makes the metric of a hole at that point",
          "[relativity][metric_field][clock]") {
    using F = io::build::Field<double>;
    const double M = 1.0, R = 12.0, w = 0.02;
    const F cx = *(F(R) * cos(F(w) * F::t())).in_spacetime();
    const F cy = *(F(R) * sin(F(w) * F::t())).in_spacetime();
    const F x = F::coord(1) - cx, y = F::coord(2) - cy, z = F::coord(3);
    const F r = sqrt(x * x + y * y + z * z);
    const std::array<F, 4> l{F(1.0), x / r, y / r, z / r};
    std::array<F, 10> g;
    for (std::size_t e = 0; e < 10; ++e) {
        const auto [i, j] = MetricField<double>::kEntries[e];
        g[e] = F(i != j ? 0.0 : (i == 0 ? -1.0 : 1.0)) + F(2.0 * M) / r * l[i] * l[j];
    }
    const auto moving = *MetricField<double>::make(g);
    const auto at_rest = kerr_schild(M, 0.0);
    for (double time : {0.0, 40.0, 157.0, 300.0})
        for (const V4 ev : {V4{time, 3.0, 4.0, 1.0}, V4{time, -20.0, 2.0, -3.0}}) {
            const V4 offset{time, ev[1] - R * std::cos(w * time), ev[2] - R * std::sin(w * time), ev[3]};
            const auto a = moving(ev), b = at_rest(offset);
            for (std::size_t i = 0; i < 4; ++i)
                for (std::size_t j = 0; j < 4; ++j)
                    CHECK_THAT(a(i, j), WithinAbs(b(i, j), 1e-12));
        }
}
