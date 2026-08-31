// Atom upgrade: Slater Z_eff, STO radial, Hund spin occupancy, multi-threaded
// rejection sampling. Exercised through the legacy header path (unchanged
// consumer API).
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <spatium/physics/elements.hpp>
#include <spatium/physics/atomic/orbital.hpp>
#include <cmath>

using namespace spatium::physics::atomic;
using Catch::Approx;

TEST_CASE("slater_z_eff: hydrogen", "[atom][slater]") {
    auto cfg = element(1).electron_config();
    // 1s electron of H: no other electrons → σ = 0, Z_eff = 1.
    REQUIRE(slater_z_eff(1, 1, 0, cfg) == Approx(1.0).margin(1e-9));
}

TEST_CASE("slater_z_eff: helium 1s", "[atom][slater]") {
    auto cfg = element(2).electron_config();
    // 1s² He: σ from the other 1s electron = 0.30, Z_eff = 2 - 0.30 = 1.70.
    REQUIRE(slater_z_eff(2, 1, 0, cfg) == Approx(1.70).margin(1e-9));
}

TEST_CASE("slater_z_eff: lithium 2s", "[atom][slater]") {
    auto cfg = element(3).electron_config();
    // Li 1s² 2s¹: outer 2s sees 2×0.85 from 1s → σ = 1.70, Z_eff = 1.30.
    REQUIRE(slater_z_eff(3, 2, 0, cfg) == Approx(1.30).margin(1e-9));
}

TEST_CASE("slater_z_eff: neon 2p", "[atom][slater]") {
    auto cfg = element(10).electron_config();
    // Ne 1s² 2s² 2p⁶: a 2p electron sees 4 same-group neighbors (other 2s/2p
    // share the (2s,2p) group, count 6 - 1 = 5 excluding self? No — 2s²+2p⁶=8
    // minus self = 7, all at 0.35) plus 2×0.85 from 1s.
    // σ = 7·0.35 + 2·0.85 = 2.45 + 1.70 = 4.15. Z_eff = 10 - 4.15 = 5.85.
    REQUIRE(slater_z_eff(10, 2, 1, cfg) == Approx(5.85).margin(1e-9));
}

TEST_CASE("slater_n_star values", "[atom][slater]") {
    REQUIRE(slater_n_star(1) == Approx(1.0));
    REQUIRE(slater_n_star(3) == Approx(3.0));
    REQUIRE(slater_n_star(4) == Approx(3.7));
    REQUIRE(slater_n_star(5) == Approx(4.0));
    REQUIRE(slater_n_star(6) == Approx(4.2));
}

TEST_CASE("slater_radial_wavefunction: finite at r=0 for s, zero for l>0", "[atom][slater]") {
    // STO with n*=1, ζ=1 at r=0: N · 0^0 · e^0 = N (finite).
    auto s_val = slater_radial_wavefunction<double>(1.0, 1.0, 0.0);
    REQUIRE(std::isfinite(s_val));
    REQUIRE(s_val > 0.0);

    // With n*=2, at r=0: r^(n*-1) = r^1 = 0 → node at origin (p-like).
    auto p_val = slater_radial_wavefunction<double>(2.0, 1.0, 0.0);
    REQUIRE(p_val == Approx(0.0).margin(1e-12));
}

TEST_CASE("fill_subshell_spins: p³ is maximally unpaired (Hund)", "[atom][spin]") {
    std::vector<SpinOccupancy> occ;
    fill_subshell_spins(1 /* l=p */, 3 /* half-filled */, occ);
    REQUIRE(occ.size() == 3);
    // Each m gets one ↑, none paired.
    for (const auto& o : occ) {
        REQUIRE(o.up == 1);
        REQUIRE(o.down == 0);
        REQUIRE(o.unpaired());
    }
}

TEST_CASE("fill_subshell_spins: p⁵ has one unpaired", "[atom][spin]") {
    std::vector<SpinOccupancy> occ;
    fill_subshell_spins(1, 5, occ);
    int paired = 0, unpaired = 0;
    for (const auto& o : occ) { if (o.paired()) ++paired; else if (o.unpaired()) ++unpaired; }
    REQUIRE(paired == 2);
    REQUIRE(unpaired == 1);
}

TEST_CASE("unpaired_electrons: Fe has 4 unpaired (paramagnetic)", "[atom][spin]") {
    // Fe: [Ar] 3d⁶ 4s². 4s² is paired. 3d⁶: 5 m's, Hund gives 5↑ + 1↓ = 4 unpaired.
    auto cfg = element(26).electron_config();
    REQUIRE(unpaired_electrons(cfg) == 4);
}

TEST_CASE("unpaired_electrons: Cr anomaly gives 6 unpaired", "[atom][spin]") {
    // Cr: [Ar] 3d⁵ 4s¹ (half-filled d + single s = 6 parallel spins).
    auto cfg = element(24).electron_config();
    REQUIRE(unpaired_electrons(cfg) == 6);
}

TEST_CASE("unpaired_electrons: noble gas = 0", "[atom][spin]") {
    REQUIRE(unpaired_electrons(element(2).electron_config()) == 0);   // He
    REQUIRE(unpaired_electrons(element(10).electron_config()) == 0);  // Ne
    REQUIRE(unpaired_electrons(element(18).electron_config()) == 0);  // Ar
}

TEST_CASE("sample_orbital_points: threaded vs single produce same count", "[atom][mt]") {
    // 1s orbital of H, same seed must yield num_points regardless of thread count.
    auto c1 = sample_orbital_points<double>(1, 0, 0, 1000, 42, 1);
    auto c2 = sample_orbital_points<double>(1, 0, 0, 1000, 42, 4);
    REQUIRE(c1.positions.size() == 1000);
    REQUIRE(c2.positions.size() == 1000);
}
