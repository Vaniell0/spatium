// Multi-source gravity (gravity_field.hpp) and its Barnes-Hut tree
// (gravity_tree.hpp).
//
// The periodic field is checked against an independent brute-force lattice
// sum, not against itself: images of every source inside a ball around the
// point, weighted by a radial step smoothed by a Gaussian of width delta
// (0.5 erfc((r - R0) / delta), R0 = 6 delta), the uniform background
// cancelling exactly because the ball is centred on the point. By Poisson
// summation the truncation error is the transform of the smoothed tail at the
// reciprocal lattice, of order exp(-(pi delta / L)^2), so delta = 1.5 L is
// about 1e-10 relative; the test shows the oracle converged by comparing two
// widths.
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/physics/mechanics/body.hpp>
#include <spatium/physics/mechanics/force.hpp>
#include <spatium/physics/mechanics/gravity_field.hpp>
#include <spatium/physics/mechanics/gravity_tree.hpp>
#include <spatium/physics/mechanics/integrator.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <vector>

using namespace spatium;
using namespace spatium::physics::mechanics;
using Catch::Approx;

namespace {

template<std::size_t N>
std::vector<GravitySource<N>> random_sources(std::size_t n, double L, unsigned seed,
                                             bool varied_mass = true) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0.0, L), m(0.5, 1.5);
    std::vector<GravitySource<N>> s(n);
    for (auto& b : s) {
        for (std::size_t i = 0; i < N; ++i) b.position[i] = u(rng);
        b.mass = varied_mass ? m(rng) : 1.0;
    }
    return s;
}

double taper(double u, double delta) {   // 1 up to 6 delta, 0 from 12 delta
    return 0.5 * std::erfc((u - 6 * delta) / delta);
}

// Acceleration at x from the periodic images of `src` inside the smoothed ball around x. G = 1.
template<std::size_t N>
Vec<double, N> lattice_oracle(const std::vector<GravitySource<N>>& src, const Vec<double, N>& x,
                              double L, unsigned p, double eps, double delta) {
    const double R = 12 * delta;
    const int M = static_cast<int>(std::ceil(R / L)) + 1;
    Vec<double, N> a{};
    for (const auto& s : src) {
        std::array<int, N> n;
        n.fill(-M);
        for (;;) {
            Vec<double, N> d;
            double r2 = 0;
            for (std::size_t i = 0; i < N; ++i) {
                d[i] = x[i] - s.position[i] - n[i] * L;
                r2 += d[i] * d[i];
            }
            const double r = std::sqrt(r2);
            if (r < R) {
                const double q = 0.5 * (p + 1);
                const double g = taper(r, delta) * s.mass * std::pow(r2 + eps * eps, -q);
                for (std::size_t i = 0; i < N; ++i) a[i] -= d[i] * g;
            }
            std::size_t i = 0;
            for (; i < N; ++i) { if (++n[i] <= M) break; n[i] = -M; }
            if (i == N) break;
        }
    }
    return a;
}

template<std::size_t N>
double rel_diff(const Vec<double, N>& a, const Vec<double, N>& b) {
    Vec<double, N> d{a - b};
    return d.norm() / std::max(b.norm(), 1e-300);
}

template<std::size_t N>
void check_against_lattice(unsigned p, double eps_over_L, std::size_t nsrc, double delta_over_L,
                           double tol) {
    const double L = 1.3;
    const double eps = eps_over_L * L;
    Vec<double, N> box;
    for (std::size_t i = 0; i < N; ++i) box[i] = L;
    GravityLaw<N> law{p, 1.0, eps, 0.0};
    GravityField<N> field(law, box);
    auto src = random_sources<N>(nsrc, L, 7 + p);
    field.set_sources(src);
    std::mt19937 rng(11);
    std::uniform_real_distribution<double> u(0.0, L);
    for (int t = 0; t < 6; ++t) {
        Vec<double, N> x;
        for (std::size_t i = 0; i < N; ++i) x[i] = u(rng);
        const auto want = lattice_oracle<N>(src, x, L, p, eps, delta_over_L * L);
        const auto want2 = lattice_oracle<N>(src, x, L, p, eps, 1.4 * delta_over_L * L);
        // the oracle itself has converged: a 30% wider taper gives the same answer
        REQUIRE(rel_diff<N>(want, want2) < 1e-9);
        REQUIRE(rel_diff<N>(field.acceleration(x), want) < tol);
    }
}

} // namespace

TEST_CASE("periodic field matches the brute-force lattice sum", "[gravity][periodic]") {
    // (dimension, exponent, softening) -- the Newtonian field in 2D and 3D, and the 1/r law in 3D.
    SECTION("2D, 1/r") { check_against_lattice<2>(1, 0.0, 5, 1.5, 1e-8); }
    SECTION("2D, 1/r, softened") { check_against_lattice<2>(1, 0.05, 5, 1.5, 1e-8); }
    SECTION("3D, 1/r^2") { check_against_lattice<3>(2, 0.0, 4, 1.5, 1e-8); }
    SECTION("3D, 1/r^2, softened") { check_against_lattice<3>(2, 0.08, 4, 1.5, 1e-8); }
    SECTION("3D, 1/r") { check_against_lattice<3>(1, 0.0, 4, 1.5, 1e-8); }
    SECTION("3D, 1/r, softened") { check_against_lattice<3>(1, 0.08, 4, 1.5, 1e-8); }
}

TEST_CASE("periodic field does not depend on the Ewald split", "[gravity][periodic]") {
    const double L = 2.0;
    auto src = random_sources<3>(20, L, 3);
    for (double eps : {0.0, 0.1}) {
        GravityLaw<3> law{2, 1.0, eps, 0.0};
        GravityField<3> a(law, {L, L, L}, EwaldOptions{4.6, 4.4});
        GravityField<3> b(law, {L, L, L}, EwaldOptions{8.0, 5.0});
        a.set_sources(src);
        b.set_sources(src);
        REQUIRE(a.reciprocal_vectors() < b.reciprocal_vectors());
        std::mt19937 rng(5);
        std::uniform_real_distribution<double> u(0.0, L);
        for (int t = 0; t < 10; ++t) {
            Vec<double, 3> x{u(rng), u(rng), u(rng)};
            REQUIRE(rel_diff<3>(a.acceleration(x), b.acceleration(x)) < 1e-7);
        }
    }
}

TEST_CASE("periodic field is periodic, continuous and odd", "[gravity][periodic]") {
    const double L = 1.0;
    GravityField<3> field(GravityLaw<3>::newtonian(1.0, 0.02), {L, L, L});
    // A configuration symmetric under x -> -x (mod L): sources in +-pairs.
    std::vector<GravitySource<3>> src;
    std::mt19937 rng(9);
    std::uniform_real_distribution<double> u(0.0, L), m(0.5, 1.5);
    for (int i = 0; i < 6; ++i) {
        Vec<double, 3> p{u(rng), u(rng), u(rng)};
        const double mass = m(rng);
        src.push_back({p, mass});
        src.push_back({Vec<double, 3>{-p[0], -p[1], -p[2]}, mass});
    }
    field.set_sources(src);

    for (int t = 0; t < 8; ++t) {
        Vec<double, 3> x{u(rng), u(rng), u(rng)};
        const auto a = field.acceleration(x);
        // Odd about the origin.
        const auto am = field.acceleration(Vec<double, 3>{-x[0], -x[1], -x[2]});
        REQUIRE(Vec<double, 3>{a + am}.norm() < 1e-10 * std::max(1.0, a.norm()));
        // Periodic: one box over, the same.
        const auto ap = field.acceleration(Vec<double, 3>{x[0] + L, x[1] - 2 * L, x[2]});
        REQUIRE(Vec<double, 3>{a - ap}.norm() < 1e-10 * std::max(1.0, a.norm()));
    }
    // Continuous across the face: x just below L against x just above 0.
    const double h = 1e-9;
    const auto lo = field.acceleration(Vec<double, 3>{L - h, 0.3137, 0.7211});
    const auto hi = field.acceleration(Vec<double, 3>{0.0 + h - 0.0, 0.3137, 0.7211});
    const auto ref = field.acceleration(Vec<double, 3>{-h, 0.3137, 0.7211});
    REQUIRE(Vec<double, 3>{lo - ref}.norm() < 1e-6 * std::max(1.0, ref.norm()));
    REQUIRE(Vec<double, 3>{hi - field.acceleration(Vec<double, 3>{h, 0.3137, 0.7211})}.norm() == 0.0);
}

TEST_CASE("an isolated periodic system exerts no net force", "[gravity][periodic]") {
    // Third law: the sum of m_i a_i vanishes, with different couplings too (a source's
    // charge is coupling * inertia, and the sum of inertia * a is what is conserved).
    const double L = 1.0;
    GravityField<3> field(GravityLaw<3>::newtonian(1.0, 0.03), {L, L, L});
    auto src = random_sources<3>(60, L, 21);
    field.set_sources(src);
    Vec<double, 3> total{};
    double scale = 0;
    for (std::size_t i = 0; i < src.size(); ++i) {
        const auto a = field.acceleration(src[i].position, i);
        total = Vec<double, 3>{total + a * src[i].mass};
        scale += src[i].mass * a.norm();
    }
    REQUIRE(total.norm() < 1e-9 * scale);
}

TEST_CASE("potential is minus the integral of the field", "[gravity][periodic]") {
    // Force against a centred difference of the potential, periodic and free.
    for (unsigned p : {1u, 2u}) {
        for (double eps : {0.0, 0.07}) {
            GravityLaw<3> law{p, 1.0, eps, 0.0};
            GravityField<3> per(law, {1.0, 1.0, 1.0});
            GravityField<3> fre(law);
            auto src = random_sources<3>(6, 1.0, 4 + p);
            per.set_sources(src);
            fre.set_sources(src);
            for (auto* f : {&per, &fre}) {
                Vec<double, 3> x{0.31, 0.62, 0.17};
                const double h = 1e-5;
                Vec<double, 3> g{};
                for (int i = 0; i < 3; ++i) {
                    Vec<double, 3> xp = x, xm = x;
                    xp[i] += h; xm[i] -= h;
                    g[i] = (f->potential(xp) - f->potential(xm)) / (2 * h);
                }
                const auto a = f->acceleration(x);
                Vec<double, 3> want{g * (-1.0)};
                REQUIRE(Vec<double, 3>{a - want}.norm() < 1e-6 * std::max(1.0, a.norm()));
            }
        }
    }
}

TEST_CASE("energy and momentum are conserved by a symplectic integrator on a torus",
          "[gravity][periodic]") {
    const double L = 1.0;
    const std::size_t n = 24;
    GravityField<3> field(GravityLaw<3>::newtonian(1.0, 0.08), {L, L, L});
    std::mt19937 rng(2);
    std::uniform_real_distribution<double> u(0.0, L), v(-0.3, 0.3);
    std::vector<PointMass<3>> b(n);
    for (auto& x : b) {
        x = PointMass<3>(1.0 / n, Vec<double, 3>{u(rng), u(rng), u(rng)},
                         Vec<double, 3>{v(rng), v(rng), v(rng)});
    }
    auto load = [&] {
        std::vector<GravitySource<3>> s;
        for (const auto& x : b) s.push_back({x.state.position, x.mass});
        field.set_sources(std::move(s));
    };
    auto energy = [&] { return total_kinetic_energy(b) + field.potential_energy(); };
    auto accel = [&] {
        std::vector<Vec<double, 3>> a(n);
        for (std::size_t i = 0; i < n; ++i) a[i] = field.acceleration(b[i].state.position, i);
        return a;
    };
    load();
    const double e0 = energy();
    const auto p0 = total_momentum(b);
    double worst = 0;
    const double dt = 2e-3;
    auto a = accel();
    for (int step = 0; step < 400; ++step) {
        for (std::size_t i = 0; i < n; ++i) {
            b[i].state.velocity = Vec<double, 3>{b[i].state.velocity + a[i] * (0.5 * dt)};
            b[i].state.position = Vec<double, 3>{b[i].state.position + b[i].state.velocity * dt};
        }
        load();
        a = accel();
        for (std::size_t i = 0; i < n; ++i)
            b[i].state.velocity = Vec<double, 3>{b[i].state.velocity + a[i] * (0.5 * dt)};
        worst = std::max(worst, std::abs(energy() - e0) / std::abs(e0));
    }
    REQUIRE(worst < 2e-3);
    REQUIRE(Vec<double, 3>{total_momentum(b) - p0}.norm() < 1e-10);
}

TEST_CASE("a body in a fixed field takes the single-body integrators", "[gravity]") {
    // GravityPull is a force of the force.hpp signature: verlet_step and
    // composite_force take it as they take PointGravity, and the coupling scales the pull.
    GravityField<3> field(GravityLaw<3>::newtonian(1.0));
    field.set_sources({{Vec<double, 3>{0, 0, 0}, 1.0}});
    PointGravity<3> point(Vec<double, 3>{0, 0, 0}, 1.0);

    PointMass<3> a(2.0, Vec<double, 3>{1, 0, 0}, Vec<double, 3>{0, 1, 0});
    PointMass<3> b = a, c = a;
    GravityPull<GravityField<3>> pull{&field, 1.0};
    GravityPull<GravityField<3>> pull2{&field, 2.0};
    const auto f1 = pull(a, 0.0), f2 = pull2(a, 0.0), fp = point(a, 0.0);
    REQUIRE(Vec<double, 3>{f1 - fp}.norm() < 1e-14);
    REQUIRE(Vec<double, 3>{f2 - fp * 2.0}.norm() < 1e-14);

    double e0 = a.kinetic_energy() - 2.0 * 1.0 / 1.0;
    auto composite = composite_force(pull, Damper<3>(0.0));
    for (int i = 0; i < 6000; ++i) {
        verlet_step(a, pull, 1e-3);
        verlet_step(b, composite, 1e-3);
        verlet_step(c, point, 1e-3);
    }
    REQUIRE(Vec<double, 3>{a.state.position - c.state.position}.norm() < 1e-12);
    REQUIRE(Vec<double, 3>{a.state.position - b.state.position}.norm() < 1e-12);
    const double e1 = a.kinetic_energy() - 2.0 / a.state.position.norm();
    REQUIRE(std::abs(e1 - e0) < 1e-5);
}

TEST_CASE("degenerate source sets", "[gravity]") {
    SECTION("no sources") {
        GravityField<3> f(GravityLaw<3>::newtonian());
        REQUIRE(f.acceleration(Vec<double, 3>{1, 2, 3}).norm() == 0.0);
        GravityField<3> g(GravityLaw<3>::newtonian(), {1.0, 1.0, 1.0});
        REQUIRE(g.acceleration(Vec<double, 3>{.1, .2, .3}).norm() == 0.0);
        GravityTree<3> t(f);
        REQUIRE(t.acceleration(Vec<double, 3>{1, 2, 3}).norm() == 0.0);
    }
    SECTION("one body feels nothing from itself, free or periodic") {
        GravityField<3> f(GravityLaw<3>::newtonian(1.0, 0.1));
        f.set_sources({{Vec<double, 3>{.5, .5, .5}, 3.0}});
        REQUIRE(f.acceleration(Vec<double, 3>{.5, .5, .5}, 0).norm() == 0.0);
        GravityField<3> g(GravityLaw<3>::newtonian(1.0, 0.1), {1.0, 1.0, 1.0});
        g.set_sources({{Vec<double, 3>{.5, .5, .5}, 3.0}});
        REQUIRE(g.acceleration(Vec<double, 3>{.5, .5, .5}, 0).norm() < 1e-12);
        GravityTree<3> t(g);
        REQUIRE(t.acceleration(Vec<double, 3>{.5, .5, .5}, 0).norm() < 1e-12);
    }
    SECTION("two bodies in free space: Newton, in 3D and 2D") {
        GravityField<3> f(GravityLaw<3>::newtonian(2.0));
        f.set_sources({{Vec<double, 3>{0, 0, 0}, 3.0}, {Vec<double, 3>{2, 0, 0}, 5.0}});
        const auto a0 = f.acceleration(Vec<double, 3>{0, 0, 0}, 0);
        REQUIRE(a0[0] == Approx(2.0 * 5.0 / 4.0));
        const auto a1 = f.acceleration(Vec<double, 3>{2, 0, 0}, 1);
        REQUIRE(a1[0] == Approx(-2.0 * 3.0 / 4.0));
        GravityField<2> g(GravityLaw<2>::newtonian(2.0));
        g.set_sources({{Vec<double, 2>{0, 0}, 3.0}, {Vec<double, 2>{2, 0}, 5.0}});
        REQUIRE(g.acceleration(Vec<double, 2>{0, 0}, 0)[0] == Approx(2.0 * 5.0 / 2.0));
        // the softer law in 3D
        GravityField<3> h(GravityLaw<3>::inverse_distance(2.0));
        h.set_sources({{Vec<double, 3>{2, 0, 0}, 5.0}});
        REQUIRE(h.acceleration(Vec<double, 3>{0, 0, 0})[0] == Approx(2.0 * 5.0 / 2.0));
    }
    SECTION("coincident bodies: softening keeps it finite, none gives zero") {
        std::vector<GravitySource<3>> s(5, GravitySource<3>{Vec<double, 3>{.3, .3, .3}, 1.0});
        GravityField<3> soft(GravityLaw<3>::newtonian(1.0, 0.05));
        soft.set_sources(s);
        REQUIRE(soft.acceleration(Vec<double, 3>{.3, .3, .3}, 0).norm() == 0.0);
        const auto a = soft.acceleration(Vec<double, 3>{.3 + 0.01, .3, .3}, 0);
        REQUIRE(std::isfinite(a.norm()));
        REQUIRE(a.norm() < 5.0 / (0.05 * 0.05));
        GravityField<3> hard(GravityLaw<3>::newtonian());
        hard.set_sources(s);
        REQUIRE(hard.acceleration(Vec<double, 3>{.3, .3, .3}).norm() == 0.0);
        // the tree builds over them and agrees
        GravityTree<3> t(soft);
        REQUIRE(Vec<double, 3>{t.acceleration(Vec<double, 3>{.3 + 0.01, .3, .3}, 0) - a}.norm()
                < 1e-12 * a.norm());
    }
    SECTION("the clamp caps the length") {
        GravityLaw<3> law = GravityLaw<3>::newtonian(1.0, 0.0);
        law.max_acceleration = 2.0;
        GravityField<3> f(law);
        f.set_sources({{Vec<double, 3>{0, 0, 0}, 10.0}});
        REQUIRE(f.acceleration(Vec<double, 3>{1, 0, 0}).norm() == Approx(2.0));
        REQUIRE(f.acceleration(Vec<double, 3>{100, 0, 0}).norm() == Approx(10.0 / 1e4));
        REQUIRE(f.acceleration_unclamped(Vec<double, 3>{1, 0, 0}).norm() == Approx(10.0));
    }
    SECTION("a periodic law outside 1..N-1 is refused") {
        REQUIRE_THROWS((GravityField<3>(GravityLaw<3>{3, 1.0, 0.0, 0.0}, {1.0, 1.0, 1.0})));
        REQUIRE_THROWS((GravityField<2>(GravityLaw<2>{2, 1.0, 0.0, 0.0}, {1.0, 1.0})));
        REQUIRE_NOTHROW((GravityField<3>(GravityLaw<3>{1, 1.0, 0.0, 0.0}, {1.0, 1.0, 1.0})));
    }
}

namespace {

struct TreeError { double mean, max_rel_rms, max_rel; };

template<std::size_t N>
TreeError tree_error(const GravityField<N>& field, double theta, std::size_t n) {
    GravityTree<N> tree(field, theta);
    const auto& src = field.sources();
    std::vector<double> err(n), mag(n);
    double rms = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const auto d = field.acceleration(src[i].position, i);
        const auto t = tree.acceleration(src[i].position, theta, i);
        err[i] = Vec<double, N>{d - t}.norm();
        mag[i] = d.norm();
        rms += mag[i] * mag[i];
    }
    rms = std::sqrt(rms / n);
    TreeError e{0, 0, 0};
    for (std::size_t i = 0; i < n; ++i) {
        e.mean += err[i] / mag[i] / n;
        e.max_rel = std::max(e.max_rel, err[i] / mag[i]);
        e.max_rel_rms = std::max(e.max_rel_rms, err[i] / rms);
    }
    return e;
}

} // namespace

TEST_CASE("tree against direct sum", "[gravity][tree]") {
    const double L = 1.0;
    SECTION("theta = 0 opens everything and is the direct sum") {
        GravityField<3> f(GravityLaw<3>::newtonian(1.0, 0.02), {L, L, L});
        f.set_sources(random_sources<3>(300, L, 1));
        REQUIRE(tree_error<3>(f, 0.0, 300).max_rel < 1e-10);
        REQUIRE(tree_error<3>(f, GravityTree<3>::default_theta, 300).mean < 3e-2);
        GravityField<2> g(GravityLaw<2>::newtonian(1.0, 0.02));
        g.set_sources(random_sources<2>(300, L, 1));
        REQUIRE(tree_error<2>(g, 0.0, 300).max_rel < 1e-10);
    }
    SECTION("error falls with theta, and the default is accurate") {
        GravityField<3> f(GravityLaw<3>::newtonian(1.0, 0.02), {L, L, L});
        f.set_sources(random_sources<3>(1500, L, 2));
        const auto e3 = tree_error<3>(f, 0.3, 1500);
        const auto e5 = tree_error<3>(f, 0.5, 1500);
        const auto e9 = tree_error<3>(f, 0.9, 1500);
        REQUIRE(e3.mean < e5.mean);
        REQUIRE(e5.mean < e9.mean);
        // On a torus the damped short-range kernel makes the monopole worse than in free space
        // at the same theta; the default keeps the mean near 1% of the force on a body.
        REQUIRE(e5.mean < 3e-2);
        REQUIRE(e5.max_rel_rms < 8e-2);
    }
    SECTION("free space, 2D torus and the 1/r law") {
        GravityField<3> f(GravityLaw<3>::newtonian(1.0, 0.02));
        f.set_sources(random_sources<3>(800, L, 3));
        REQUIRE(tree_error<3>(f, 0.5, 800).mean < 8e-3);
        GravityField<2> g(GravityLaw<2>::newtonian(1.0, 0.02), {L, L});
        g.set_sources(random_sources<2>(800, L, 3));
        REQUIRE(tree_error<2>(g, 0.5, 800).mean < 8e-2);
        GravityField<3> h(GravityLaw<3>::inverse_distance(1.0, 0.02), {L, L, L});
        h.set_sources(random_sources<3>(800, L, 3));
        REQUIRE(tree_error<3>(h, 0.5, 800).mean < 6e-2);
    }
    SECTION("moments: the root holds the total mass and the centre of mass") {
        GravityField<3> f(GravityLaw<3>::newtonian());
        auto s = random_sources<3>(200, L, 4);
        f.set_sources(s);
        GravityTree<3> t(f);
        double m = 0;
        Vec<double, 3> c{};
        for (const auto& x : s) { m += x.mass; c = Vec<double, 3>{c + x.position * x.mass}; }
        c = Vec<double, 3>{c / m};
        REQUIRE(t.moments()[0].mass == Approx(m));
        REQUIRE(Vec<double, 3>{t.moments()[0].centre - c}.norm() < 1e-12);
    }
}
