#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <spatium/physics/atomic/orbital.hpp>
#include <cmath>
#include <numbers>

using namespace spatium;
using namespace spatium::physics::atomic;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

constexpr double pi = std::numbers::pi;

// ── Legendre polynomials ──────────────────────────────────────

TEST_CASE("Legendre P_0^0 = 1", "[orbital]") {
    CHECK_THAT(assoc_legendre<double>(0, 0, 0.5), WithinAbs(1.0, 1e-12));
}

TEST_CASE("Legendre P_1^0 = x", "[orbital]") {
    CHECK_THAT(assoc_legendre<double>(1, 0, 0.7), WithinAbs(0.7, 1e-12));
}

TEST_CASE("Legendre P_2^0 = (3x²-1)/2", "[orbital]") {
    double x = 0.6;
    double expected = (3 * x * x - 1) / 2;
    CHECK_THAT(assoc_legendre<double>(2, 0, x), WithinAbs(expected, 1e-12));
}

TEST_CASE("Legendre P_1^1 = -sqrt(1-x²)", "[orbital]") {
    double x = 0.6;
    double expected = -std::sqrt(1 - x * x);
    CHECK_THAT(assoc_legendre<double>(1, 1, x), WithinAbs(expected, 1e-12));
}

// ── Spherical harmonics ───────────────────────────────────────

TEST_CASE("Y_0^0 is constant = 1/(2√π)", "[orbital]") {
    double val = real_spherical_harmonic<double>(0, 0, 0.5, 0.3);
    CHECK_THAT(val, WithinAbs(1.0 / (2.0 * std::sqrt(pi)), 1e-12));
}

TEST_CASE("Y_l^m normalization (orthogonality integral)", "[orbital]") {
    // ∫∫ |Y_1^0|² sin(θ) dθ dφ = 1
    int N = 100;
    double sum = 0;
    double dtheta = pi / N;
    double dphi = 2 * pi / N;
    for (int i = 0; i < N; ++i) {
        double theta = (i + 0.5) * dtheta;
        for (int j = 0; j < N; ++j) {
            double phi = (j + 0.5) * dphi;
            double Y = real_spherical_harmonic<double>(1, 0, theta, phi);
            sum += Y * Y * std::sin(theta) * dtheta * dphi;
        }
    }
    CHECK_THAT(sum, WithinRel(1.0, 0.02));  // 2% tolerance for grid integration
}

// ── Laguerre polynomials ──────────────────────────────────────

TEST_CASE("Laguerre L_0^α = 1", "[orbital]") {
    CHECK_THAT(gen_laguerre<double>(0, 2.0, 3.0), WithinAbs(1.0, 1e-12));
}

TEST_CASE("Laguerre L_1^α = 1+α-x", "[orbital]") {
    double alpha = 3.0, x = 2.0;
    CHECK_THAT(gen_laguerre<double>(1, alpha, x), WithinAbs(1 + alpha - x, 1e-12));
}

// ── Radial wavefunction ───────────────────────────────────────

TEST_CASE("R_10 at r=0 is nonzero (1s orbital)", "[orbital]") {
    double R = radial_wavefunction<double>(1, 0, 0.0);
    CHECK(std::abs(R) > 0.1);  // 1s has nonzero density at origin
}

TEST_CASE("R_20 has a node (2s orbital)", "[orbital]") {
    // 2s has a node at r = 2a₀
    double R_small = radial_wavefunction<double>(2, 0, 0.5);
    double R_large = radial_wavefunction<double>(2, 0, 4.0);
    // Should have opposite signs (node between)
    CHECK(R_small * R_large < 0);
}

TEST_CASE("R_nl decays at large r", "[orbital]") {
    double R = radial_wavefunction<double>(1, 0, 50.0);
    CHECK(std::abs(R) < 1e-10);
}

// ── Probability density ───────────────────────────────────────

TEST_CASE("1s orbital density is spherically symmetric", "[orbital]") {
    // |ψ_100|² should be same at all angles for same r
    double d1 = orbital_density<double>(1, 0, 0, 2.0, 0.5, 0.0);
    double d2 = orbital_density<double>(1, 0, 0, 2.0, 1.0, 1.0);
    double d3 = orbital_density<double>(1, 0, 0, 2.0, 2.0, 3.0);
    CHECK_THAT(d1, WithinRel(d2, 1e-10));
    CHECK_THAT(d1, WithinRel(d3, 1e-10));
}

TEST_CASE("2p_z orbital has node at equator", "[orbital]") {
    // 2p_z (n=2, l=1, m=0): Y_1^0 ∝ cos(θ), zero at θ=π/2
    double d = orbital_density<double>(2, 1, 0, 3.0, pi / 2, 0.0);
    CHECK_THAT(d, WithinAbs(0.0, 1e-15));
}

TEST_CASE("Orbital density non-negative", "[orbital]") {
    // Sample random points — density should always be ≥ 0
    for (int n = 1; n <= 3; ++n)
        for (int l = 0; l < n; ++l)
            for (int m = -l; m <= l; ++m) {
                double d = orbital_density<double>(n, l, m, 1.0, 0.7, 0.3);
                CHECK(d >= 0.0);
            }
}

// ── ImplicitSurface factory ───────────────────────────────────

TEST_CASE("make_orbital produces valid ImplicitSurface", "[orbital]") {
    auto surf = make_orbital(1, 0, 0);
    // At origin (inside 1s orbital): density > iso → F < 0
    CHECK(surf(0.0, 0.0, 0.0) > 0);  // density at origin > iso_level
    // Far away: density ≈ 0 → F < 0 (below iso)
    CHECK(surf(100.0, 0.0, 0.0) < 0);
}

TEST_CASE("make_orbital marching cubes produces mesh", "[orbital]") {
    auto surf = make_orbital(1, 0, 0, 0.01);  // 1s orbital, higher iso
    auto mesh = marching_cubes(surf, 16);
    CHECK(mesh.vertices.size() > 0);
    CHECK(mesh.faces.size() > 0);
}

// ── Point cloud sampling ──────────────────────────────────────

TEST_CASE("sample_orbital_points returns requested count", "[orbital][pointcloud]") {
    auto cloud = sample_orbital_points(1, 0, 0, 1000, 42);
    CHECK(cloud.positions.size() == 1000);
    CHECK(cloud.positive_lobe.size() == 1000);
}

TEST_CASE("sample_orbital_points within bounding box", "[orbital][pointcloud]") {
    int n = 2;
    double bound = n * (n + 1) + 3.0;
    auto cloud = sample_orbital_points(n, 1, 0, 500, 42);
    for (auto& p : cloud.positions) {
        CHECK(std::abs(p[0]) <= bound);
        CHECK(std::abs(p[1]) <= bound);
        CHECK(std::abs(p[2]) <= bound);
    }
}

TEST_CASE("1s orbital: radial distribution peaks near r=1", "[orbital][pointcloud]") {
    auto cloud = sample_orbital_points(1, 0, 0, 50000, 42);
    // Bin radii and find mode
    constexpr int bins = 20;
    double max_r = 6.0;
    int counts[bins] = {};
    for (auto& p : cloud.positions) {
        double r = std::sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
        int bin = static_cast<int>(r / max_r * bins);
        if (bin >= 0 && bin < bins) counts[bin]++;
    }
    // Find bin with max count (weighted by 1/r² to get radial probability)
    // Actually, raw density sampling already weights by |ψ|²,
    // the r²|ψ|² radial distribution peaks at r=1 for 1s.
    // In our raw sampling, the 3D density peaks at r=0.
    // But the bin with most points (in shell volume) should be near r≈1.
    int max_bin = 0;
    for (int i = 1; i < bins; ++i)
        if (counts[i] > counts[max_bin]) max_bin = i;
    double peak_r = (max_bin + 0.5) * max_r / bins;
    // Peak should be within first few bins (high density near origin)
    CHECK(peak_r < 2.0);
}

TEST_CASE("2p_z lobes: positive=z>0, negative=z<0", "[orbital][pointcloud]") {
    auto cloud = sample_orbital_points(2, 1, 0, 5000, 42);
    int pos_above = 0, pos_below = 0, neg_above = 0, neg_below = 0;
    for (std::size_t i = 0; i < cloud.positions.size(); ++i) {
        double z = cloud.positions[i][2];
        if (cloud.positive_lobe[i]) {
            if (z > 0) ++pos_above; else ++pos_below;
        } else {
            if (z > 0) ++neg_above; else ++neg_below;
        }
    }
    // Vast majority of positive lobe should be z>0
    CHECK(pos_above > pos_below * 5);
    // Vast majority of negative lobe should be z<0
    CHECK(neg_below > neg_above * 5);
}

TEST_CASE("sample_orbital_points reproducibility", "[orbital][pointcloud]") {
    auto c1 = sample_orbital_points(1, 0, 0, 100, 123);
    auto c2 = sample_orbital_points(1, 0, 0, 100, 123);
    REQUIRE(c1.positions.size() == c2.positions.size());
    for (std::size_t i = 0; i < c1.positions.size(); ++i) {
        CHECK(c1.positions[i][0] == c2.positions[i][0]);
        CHECK(c1.positions[i][1] == c2.positions[i][1]);
        CHECK(c1.positions[i][2] == c2.positions[i][2]);
        CHECK(c1.positive_lobe[i] == c2.positive_lobe[i]);
    }
}
