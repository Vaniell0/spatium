// Metric-agnostic geodesic engine: exact Christoffel symbols via Dual<T>,
// RK4 integration of the full 4-coordinate geodesic equation.
// Validation, not just unit coverage:
//   - Cross-check against the plane-confined u(phi) equation examples/
//     blackhole_demo.cpp integrates directly (Verlet) -- that ODE is
//     exact for Schwarzschild given a correct impact parameter b, so
//     agreement here validates the new 4D/Dual<T> machinery against a
//     known-good, independently-implemented method.
//   - Energy/angular-momentum conservation drift over a real trajectory,
//     same style as physics/mechanics's Kepler-orbit and LGVI checks.
//   - Closed-form checks: capture/escape flips across the photon-sphere
//     critical impact parameter b_crit = 3*sqrt(3)*M; weak-field
//     deflection converges toward 4M/b as b/M grows (verified via a
//     standalone probe before locking tolerances here -- 26%/11%/6.4%
//     relative error at b/M = 15/30/50, consistent with the known
//     next-order O(M/b) correction to the leading-order formula).
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/physics/relativity/schwarzschild.hpp>
#include <spatium/physics/relativity/kerr.hpp>
#include <spatium/physics/relativity/geodesic.hpp>
#include <spatium/physics/relativity/accretion_disk.hpp>
#include <christoffel_closed_form.hpp>
#include <cmath>
#include <numbers>

using Catch::Matchers::WithinAbs;

using namespace spatium;
using namespace spatium::physics::relativity;

namespace {

// Mirrors examples/blackhole_demo.cpp's trace_geodesic(): velocity-Verlet
// on the plane-confined photon equation u'' = -u + 3*M*u^2, u = 1/r. Not
// the code under test -- an independent reference to cross-check against.
double old_u_phi_delta_phi(double M, double D, double b, double R_escape) {
    constexpr double DPHI = 0.0005;
    constexpr long MAX_STEPS = 4'000'000;
    double u = 1.0 / D;
    double up = std::sqrt(std::max(0.0, 1.0 / (b * b) - u * u + 2.0 * M * u * u * u));
    double phi = 0.0;
    const double U_CAPTURE = 1.0 / (2.0 * M * 1.001);
    const double U_ESCAPE = 1.0 / R_escape;
    auto upp_of = [M](double uu) { return -uu + 3.0 * M * uu * uu; };
    double upp = upp_of(u);
    for (long i = 0; i < MAX_STEPS; ++i) {
        double u_new = u + DPHI * up + 0.5 * DPHI * DPHI * upp;
        double upp_new = upp_of(u_new);
        double up_new = up + 0.5 * DPHI * (upp + upp_new);
        u = u_new; up = up_new; upp = upp_new;
        phi += DPHI;
        if (u >= U_CAPTURE) return -1.0;  // captured
        if (u <= U_ESCAPE && up < 0.0) return phi;
    }
    return phi;
}

// Sets up an inward-aimed photon at radius D with impact parameter b
// (E=1, L=b -- an affine-parameter gauge choice, doesn't affect the
// trajectory's shape) and integrates via the new engine until capture
// (r <= horizon*1.001) or escape (r >= R_escape, moving outward).
// Returns final phi, or -1.0 if captured.
double new_geodesic_delta_phi(double M, double D, double b, double R_escape, double dlambda) {
    constexpr long MAX_STEPS = 5'000'000;
    SchwarzschildMetric<double> metric{M};
    double f0 = 1.0 - 2.0 * M / D;
    double E = 1.0, L = b;
    double ur0 = std::sqrt(std::max(0.0, E * E - f0 * L * L / (D * D)));
    Vec<double, 8> state{0.0, D, std::numbers::pi / 2.0, 0.0, E / f0, -ur0, 0.0, L / (D * D)};
    const double R_capture = 2.0 * M * 1.001;
    for (long i = 0; i < MAX_STEPS; ++i) {
        state = geodesic_step(metric, state, dlambda);
        double r = state[1];
        if (r <= R_capture) return -1.0;  // captured
        if (r >= R_escape && state[5] > 0.0) return state[3];
    }
    return state[3];
}

// Generic (T, Metric)-agnostic version of new_geodesic_delta_phi() above,
// for the fp32-vs-fp64 precision check below: is a T=float geodesic
// integration close enough to T=double to be worth running on hardware
// (a T4 GPU) where fp64 throughput is 1/32 of fp32? This is the decision
// gate for whether the planned CUDA port (see project memory) can use
// fp32 at all, including right at the numerically hardest region (near
// b_crit, the photon-sphere critical impact parameter) -- not just
// asserted, measured, including AT that boundary.
template<Scalar T, typename Metric>
struct PrecisionTrace { bool captured; T phi; };

template<Scalar T, typename Metric>
PrecisionTrace<T, Metric> trace_for_precision(const Metric& metric, T M, T D, T b, T R_capture,
                                               T R_escape, T dlambda, long max_steps) {
    using std::sqrt;
    T f0 = T(1) - T(2) * M / D;
    T E = T(1), L = b;
    T ur0 = sqrt(std::max(T(0), E * E - f0 * L * L / (D * D)));
    Vec<T, 8> state{T(0), D, std::numbers::pi_v<T> / T(2), T(0), E / f0, -ur0, T(0), L / (D * D)};
    for (long i = 0; i < max_steps; ++i) {
        state = geodesic_step(metric, state, dlambda);
        T r = state[1];
        if (!std::isfinite(r) || r <= R_capture) return {true, state[3]};
        if (r >= R_escape && state[5] > T(0)) return {false, state[3]};
    }
    return {true, state[3]};
}

} // namespace

TEST_CASE("Schwarzschild geodesic integration: fp32 matches fp64 capture/escape classification "
          "even near b_crit, and delta_phi agreement stays within ~1e-3 relative",
          "[relativity][precision]") {
    double M = 1.0, D = 30.0, R_escape = 50.0;
    double b_crit = schwarzschild_critical_impact_parameter(M);
    SchwarzschildMetric<double> metric_d{M};
    SchwarzschildMetric<float> metric_f{float(M)};
    double R_capture_d = 2.0 * M * 1.15;
    float R_capture_f = float(R_capture_d);

    // Ratios spanning safely-escaping to safely-captured, including
    // several right at the critical curve where the trajectory winds
    // chaotically many times -- exactly the region a GPU port's
    // precision has to survive, not just the easy cases.
    for (double ratio : {0.5, 0.9, 0.99, 0.999, 0.9999, 1.0001, 1.001, 1.01, 1.1, 2.0, 5.0}) {
        double b = ratio * b_crit;
        auto rd = trace_for_precision<double>(metric_d, M, D, b, R_capture_d, R_escape, 0.01, 200000);
        auto rf = trace_for_precision<float>(metric_f, float(M), float(D), float(b), R_capture_f,
                                              float(R_escape), 0.01f, 200000);
        INFO("ratio=" << ratio << " b=" << b);
        CHECK(rd.captured == rf.captured);
        if (!rd.captured && !rf.captured) {
            double rel_err = std::abs(rd.phi - double(rf.phi)) / std::max(1e-9, std::abs(rd.phi));
            CHECK(rel_err < 1e-3);
        }
    }
}

TEST_CASE("Kerr geodesic integration: fp32 matches fp64 capture/escape classification, "
          "delta_phi agreement stays within ~1e-3 relative",
          "[relativity][kerr][precision]") {
    double M = 1.0, a = 0.9, D = 30.0, R_escape = 50.0;
    KerrMetric<double> metric_d{M, a};
    KerrMetric<float> metric_f{float(M), float(a)};
    double R_capture_d = kerr_outer_horizon_radius(M, a) * 1.15;
    float R_capture_f = float(R_capture_d);

    for (double b : {2.0, 3.0, 3.5, 3.6, 3.65, 3.7, 3.8, 4.0, 5.0, 8.0}) {
        auto rd = trace_for_precision<double>(metric_d, M, D, b, R_capture_d, R_escape, 0.01, 200000);
        auto rf = trace_for_precision<float>(metric_f, float(M), float(D), float(b), R_capture_f,
                                              float(R_escape), 0.01f, 200000);
        INFO("b=" << b);
        CHECK(rd.captured == rf.captured);
        if (!rd.captured && !rf.captured) {
            double rel_err = std::abs(rd.phi - double(rf.phi)) / std::max(1e-9, std::abs(rd.phi));
            CHECK(rel_err < 1e-3);
        }
    }
}

TEST_CASE("geodesic engine matches the plane-confined u(phi) reference", "[relativity]") {
    double M = 1.0, D = 30.0, b = 6.0, R_escape = 50.0;
    double phi_old = old_u_phi_delta_phi(M, D, b, R_escape);
    double phi_new = new_geodesic_delta_phi(M, D, b, R_escape, 0.01);
    REQUIRE(phi_old > 0.0);
    REQUIRE(phi_new > 0.0);
    double rel_diff = std::abs(phi_old - phi_new) / phi_old;
    REQUIRE(rel_diff < 1e-3);
}

TEST_CASE("geodesic engine conserves energy and angular momentum", "[relativity]") {
    double M = 1.0, D = 30.0, b = 6.0;
    SchwarzschildMetric<double> metric{M};
    double f0 = 1.0 - 2.0 * M / D;
    double E = 1.0, L = b;
    double ur0 = std::sqrt(std::max(0.0, E * E - f0 * L * L / (D * D)));
    Vec<double, 8> state{0.0, D, std::numbers::pi / 2.0, 0.0, E / f0, -ur0, 0.0, L / (D * D)};

    double E0 = killing_energy(metric, state);
    double L0 = killing_angular_momentum(metric, state);

    constexpr double dlambda = 0.01;
    for (int i = 0; i < 20000; ++i)
        state = geodesic_step(metric, state, dlambda);

    double Ef = killing_energy(metric, state);
    double Lf = killing_angular_momentum(metric, state);

    REQUIRE(std::abs(Ef - E0) < 1e-9);
    REQUIRE(std::abs(Lf - L0) < 1e-9);
}

TEST_CASE("capture/escape flips across the photon-sphere critical impact parameter", "[relativity]") {
    double M = 1.0, D = 30.0, R_escape = 50.0;
    double b_crit = schwarzschild_critical_impact_parameter(M);

    double phi_below = new_geodesic_delta_phi(M, D, 0.99 * b_crit, R_escape, 0.01);
    double phi_above = new_geodesic_delta_phi(M, D, 1.01 * b_crit, R_escape, 0.01);

    REQUIRE(phi_below < 0.0);  // captured
    REQUIRE(phi_above > 0.0);  // escaped
}

TEST_CASE("weak-field deflection converges toward 4M/b", "[relativity]") {
    double M = 1.0, b = 50.0, D = 200.0, R_escape = 2000.0;
    double phi_gr = new_geodesic_delta_phi(M, D, b, R_escape, 0.02);
    double phi_flat = new_geodesic_delta_phi(1e-12, D, b, R_escape, 0.02);
    REQUIRE(phi_gr > 0.0);
    REQUIRE(phi_flat > 0.0);

    double deflection = phi_gr - phi_flat;
    double predicted = 4.0 * M / b;
    double rel_err = std::abs(deflection - predicted) / predicted;
    REQUIRE(rel_err < 0.10);  // leading-order formula; ~6% measured at this b/M
}

TEST_CASE("circular_orbit_ut matches the closed-form value at the ISCO", "[relativity][disk]") {
    double M = 1.0;
    double r_isco = schwarzschild_isco_radius(M);
    REQUIRE(r_isco == 6.0);

    // u^t = 1/sqrt(1 - 3M/r); at r=6M this is 1/sqrt(1-0.5) = sqrt(2).
    double ut = circular_orbit_ut(M, r_isco);
    REQUIRE_THAT(ut, WithinAbs(std::numbers::sqrt2, 1e-12));
}

TEST_CASE("disk_redshift_factor reduces to the pure gravitational case at b=0", "[relativity][disk]") {
    double M = 1.0, r = 10.0;
    double ut = circular_orbit_ut(M, r);
    REQUIRE_THAT(disk_redshift_factor(M, r, 0.0), WithinAbs(ut, 1e-12));
}

TEST_CASE("disk_redshift_factor is asymmetric between the approaching and receding limbs",
          "[relativity][disk]") {
    double M = 1.0, r = 10.0, b = 3.0;
    double ut = circular_orbit_ut(M, r);
    double z_prograde = disk_redshift_factor(M, r, b);    // photon co-rotating with the disk
    double z_retrograde = disk_redshift_factor(M, r, -b); // photon counter-rotating

    REQUIRE(z_prograde < ut);       // approaching side: Doppler blueshift reduces (1+z)
    REQUIRE(z_retrograde > ut);     // receding side: Doppler redshift adds to it
    REQUIRE(z_prograde < z_retrograde);
}

TEST_CASE("KerrMetric(a=0) matches SchwarzschildMetric exactly", "[relativity][kerr]") {
    double M = 1.0;
    SchwarzschildMetric<double> schw{M};
    KerrMetric<double> kerr{M, 0.0};
    Vec<double, 4> x{0.0, 8.0, 1.1, 0.4};  // generic off-equatorial point

    Matrix<double, 4, 4> gs = schw(x);
    Matrix<double, 4, 4> gk = kerr(x);
    for (std::size_t mu = 0; mu < 4; ++mu)
        for (std::size_t nu = 0; nu < 4; ++nu)
            CHECK_THAT(gk(mu, nu), WithinAbs(gs(mu, nu), 1e-12));
}

TEST_CASE("kerr_outer_horizon_radius(a=0) is 2M", "[relativity][kerr]") {
    REQUIRE_THAT(kerr_outer_horizon_radius(1.0, 0.0), WithinAbs(2.0, 1e-12));
}

TEST_CASE("kerr_photon_orbit_radius and kerr_isco_radius reduce to Schwarzschild at a=0",
          "[relativity][kerr]") {
    double M = 1.0;
    CHECK_THAT(kerr_photon_orbit_radius(M, 0.0, true), WithinAbs(3.0, 1e-9));
    CHECK_THAT(kerr_photon_orbit_radius(M, 0.0, false), WithinAbs(3.0, 1e-9));
    CHECK_THAT(kerr_isco_radius(M, 0.0, true), WithinAbs(6.0, 1e-9));
    CHECK_THAT(kerr_isco_radius(M, 0.0, false), WithinAbs(6.0, 1e-9));
}

TEST_CASE("kerr_photon_orbit_radius and kerr_isco_radius match the known extremal (a=M) values",
          "[relativity][kerr]") {
    // Textbook extremal-Kerr results (a=M): prograde photon orbit and
    // ISCO both collapse onto the horizon at r=M; retrograde photon
    // orbit sits at 4M, retrograde ISCO at 9M. kerr_isco_radius's z1
    // term involves cbrt(1-a_star^2), which converges to its a=M limit
    // only as the CUBE ROOT of (1-a_star) -- an a=M-1e-9 probe measured
    // ~1.6e-3 residual on the prograde ISCO, consistent with
    // (1e-9)^(1/3) ~ 1e-3, not a bug. 1e-2 tolerance at a=M-1e-9 gives
    // real margin over that expected residual while still being a tight
    // check against the wrong (non-extremal) value.
    double M = 1.0, a = 1.0 - 1e-9;  // just inside extremal to avoid the a=M singular endpoint
    CHECK_THAT(kerr_photon_orbit_radius(M, a, true), WithinAbs(1.0, 1e-4));
    CHECK_THAT(kerr_photon_orbit_radius(M, a, false), WithinAbs(4.0, 1e-4));
    CHECK_THAT(kerr_isco_radius(M, a, true), WithinAbs(1.0, 1e-2));
    CHECK_THAT(kerr_isco_radius(M, a, false), WithinAbs(9.0, 1e-2));
}

TEST_CASE("kerr_equatorial_omega/ut reduce to Schwarzschild's circular-orbit values at a=0",
          "[relativity][kerr]") {
    double M = 1.0, r = 10.0;
    CHECK_THAT(kerr_equatorial_omega(M, 0.0, r, true), WithinAbs(keplerian_omega(M, r), 1e-12));
    CHECK_THAT(kerr_equatorial_ut(M, 0.0, r, true), WithinAbs(circular_orbit_ut(M, r), 1e-9));
    CHECK_THAT(kerr_equatorial_ut(M, 0.0, r, false), WithinAbs(circular_orbit_ut(M, r), 1e-9));
}

TEST_CASE("kerr_disk_redshift_factor reduces to disk_redshift_factor at a=0",
          "[relativity][kerr][disk]") {
    double M = 1.0, r = 10.0, b = 3.0;
    CHECK_THAT(kerr_disk_redshift_factor(M, 0.0, r, b, true),
               WithinAbs(disk_redshift_factor(M, r, b), 1e-9));
    CHECK_THAT(kerr_disk_redshift_factor(M, 0.0, r, -b, true),
               WithinAbs(disk_redshift_factor(M, r, -b), 1e-9));
}

// gpu/christoffel_closed_form.hpp's hand/symbolically-derived, straight-
// line expressions are meant to replace geodesic.hpp's generic Dual<T> +
// general invert() machinery inside the GPU kernel (Phase 2 of the CUDA
// port) -- SIMT execution punishes the indirection/branching that
// genericity costs, and a kernel that will only ever run these two
// published metrics doesn't need to pay for supporting a hypothetical
// future one. Before any CUDA is written at all, the closed-form
// version has to reproduce the already-validated generic engine
// component-by-component, not just "look right" against a textbook.
// Both are exact analytic evaluations (no discretization on either
// side, unlike the fp32-vs-fp64 precision cases above), so agreement is
// expected to near machine precision everywhere except right next to
// genuine coordinate features (Delta -> 0 at the horizon) where the
// two independently-structured expressions can differ in how they
// accumulate rounding -- checked with a relative tolerance that scales
// with each component's own magnitude, not a single absolute bound.
namespace {
void check_christoffel_match(const std::array<Matrix<double, 4, 4>, 4>& ref,
                              const double (&cf)[4][4][4], double a, double r, double theta) {
    for (std::size_t lam = 0; lam < 4; ++lam)
        for (std::size_t mu = 0; mu < 4; ++mu)
            for (std::size_t nu = 0; nu < 4; ++nu) {
                double got = cf[lam][mu][nu];
                double want = ref[lam](mu, nu);
                double scale = std::max(1.0, std::abs(want));
                INFO("a=" << a << " r=" << r << " theta=" << theta << " lambda=" << lam
                          << " mu=" << mu << " nu=" << nu << " got=" << got << " want=" << want);
                CHECK(std::abs(got - want) <= 1e-6 * scale);
            }
}
}  // namespace

TEST_CASE("schwarzschild_christoffel_closed_form matches the generic Dual<T> christoffel()",
          "[relativity][gpu][closed-form]") {
    double M = 1.0;
    SchwarzschildMetric<double> metric{M};
    for (double r : {2.2, 2.5, 3.0, 4.0, 6.0, 10.0, 25.0}) {
        for (double theta : {0.3, 1.0, std::numbers::pi / 2.0, 2.0, 2.8}) {
            Vec<double, 4> x{0.0, r, theta, 0.0};
            auto Gamma_ref = christoffel(metric, x);
            double Gamma_cf[4][4][4];
            spatium::gpu::schwarzschild_christoffel_closed_form(M, r, theta, Gamma_cf);
            check_christoffel_match(Gamma_ref, Gamma_cf, 0.0, r, theta);
        }
    }
}

TEST_CASE("kerr_christoffel_closed_form matches the generic Dual<T> christoffel()",
          "[relativity][kerr][gpu][closed-form]") {
    double M = 1.0;
    for (double a : {0.0, 0.5, 0.9, 0.99}) {
        KerrMetric<double> metric{M, a};
        double r_horizon = kerr_outer_horizon_radius(M, a);
        for (double r_over_horizon : {1.2, 1.5, 2.0, 3.0, 6.0, 15.0}) {
            double r = r_horizon * r_over_horizon;
            for (double theta : {0.3, 1.0, std::numbers::pi / 2.0, 2.0, 2.8}) {
                Vec<double, 4> x{0.0, r, theta, 0.0};
                auto Gamma_ref = christoffel(metric, x);
                double Gamma_cf[4][4][4];
                spatium::gpu::kerr_christoffel_closed_form(M, a, r, theta, Gamma_cf);
                check_christoffel_match(Gamma_ref, Gamma_cf, a, r, theta);
            }
        }
    }
}
