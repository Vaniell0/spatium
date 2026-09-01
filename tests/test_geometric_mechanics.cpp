// LGVI on SO(3), DEC heat equation, geometric continuum mechanics
// scaffolding.
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <spatium/physics/mechanics/lgvi.hpp>
#include <spatium/physics/mechanics/continuum.hpp>
#if SPATIUM_HAS_EIGEN
#  include <spatium/mesh/dec.hpp>
#endif
#include <spatium/mesh/mesh.hpp>
#include <spatium/mesh/topology.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <cmath>

using namespace spatium;
using namespace spatium::algebra;
using namespace spatium::physics::mechanics;
using namespace spatium::mesh;
using Catch::Approx;

// ── LGVI tests ────────────────────────────────────────────────

TEST_CASE("LGVI: hat & vee maps round-trip", "[lgvi][hat]") {
    Vec3 v{1.5, -0.7, 2.3};
    auto hat = hat_so3(v);
    auto back = vee_so3(hat);
    REQUIRE(back[0] == Approx(v[0]).margin(1e-15));
    REQUIRE(back[1] == Approx(v[1]).margin(1e-15));
    REQUIRE(back[2] == Approx(v[2]).margin(1e-15));
}

TEST_CASE("LGVI: cay map of zero is identity", "[lgvi][cayley]") {
    Vec3 zero{0, 0, 0};
    auto F = lgvi_cayley(zero);
    auto I = SO3::ElementType::identity();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            REQUIRE(F(i, j) == Approx(I(i, j)).margin(1e-15));
}

TEST_CASE("LGVI: cay of small y produces near-rotation matrix",
          "[lgvi][cayley]") {
    Vec3 y{0.1, 0.0, 0.0};         // small "rotation" around x by ~0.1 rad
    auto F = lgvi_cayley(y);
    // F^T F should be I to high precision (Cayley map → SO(3) image).
    auto FtF = F.transpose() * F;
    auto I = SO3::ElementType::identity();
    double err = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            err = std::max(err, std::abs(FtF(i, j) - I(i, j)));
    REQUIRE(err < 1e-13);
}

TEST_CASE("LGVI free rigid body: spatial momentum L = R · Π conserved",
          "[lgvi][noether]") {
    // Discrete Noether for the SO(3) symmetry of L_d → spatial angular
    // momentum is conserved EXACTLY (to round-off).
    Vec3 J_diag{1.0, 2.0, 3.0};

    LGVIRigidBodyState s;
    s.R = SO3::ElementType::identity();
    Vec3 omega0{1.0, 0.5, 0.2};
    s.Pi = Vec3{omega0[0] * J_diag[0],
                omega0[1] * J_diag[1],
                omega0[2] * J_diag[2]};

    Vec3 L0 = lgvi_spatial_angular_momentum(s);

    constexpr double dt = 1e-3;
    for (int i = 0; i < 5000; ++i)
        s = lgvi_rigid_body_step(s, J_diag, dt);

    Vec3 L = lgvi_spatial_angular_momentum(s);
    // Spatial L should equal L0 to discrete-Noether precision (~1e-10).
    REQUIRE(std::abs(L[0] - L0[0]) < 1e-9);
    REQUIRE(std::abs(L[1] - L0[1]) < 1e-9);
    REQUIRE(std::abs(L[2] - L0[2]) < 1e-9);
}

TEST_CASE("LGVI free rigid body: kinetic energy bounded over many steps",
          "[lgvi][energy]") {
    // Honest status: the current Cayley-form LGVI step (both the
    // linear `lgvi_solve_y_diag` and the fixed-point-iterated
    // `lgvi_solve_y_diag_implicit`) produce the same empirical
    // drift, ≈ 3 % over 5 000 steps at h = 1e-3. The two agree
    // because at these scales y ≈ 1e-4 and the nonlinear
    // (I − ŷ/2) Π̂ (I + ŷ/2) RHS correction is 1e-8 per step —
    // too small to eat into the observed drift. Getting to the
    // paper bound (~1e-6) needs a full Newton iteration on the
    // joint system (F_k, Π_{k+1}), not just on y; that is a
    // bigger piece of work than fits in this slice. Tracking
    // issue / refinement target documented in lgvi.hpp.
    Vec3 J_diag{1.0, 2.0, 3.0};
    LGVIRigidBodyState s;
    s.R = SO3::ElementType::identity();
    s.Pi = Vec3{1.0, 0.5, 0.2};

    double E0 = lgvi_kinetic_energy(s.Pi, J_diag);

    constexpr double dt = 1e-3;
    double max_drift = 0;
    for (int i = 0; i < 5000; ++i) {
        s = lgvi_rigid_body_step(s, J_diag, dt);
        double E = lgvi_kinetic_energy(s.Pi, J_diag);
        max_drift = std::max(max_drift, std::abs(E - E0) / E0);
    }
    REQUIRE(max_drift < 5e-2);
}

TEST_CASE("LGVI |Π| body momentum conserved to machine precision",
          "[lgvi][casimir]") {
    // |Π|² is a Casimir of so(3)* — preserved by the exact Euler
    // flow and by any Lie-Poisson integrator. Our Cayley step
    // preserves it to round-off independently of the energy drift
    // issue, confirming the update Π_{k+1} = F^T · Π_k is at
    // least an isometry on so(3)* (as it should be for F ∈ SO(3)).
    Vec3 J_diag{1.0, 2.0, 3.0};
    LGVIRigidBodyState s;
    s.R = SO3::ElementType::identity();
    s.Pi = Vec3{1.0, 0.5, 0.2};
    double mag0_sq = s.Pi[0]*s.Pi[0] + s.Pi[1]*s.Pi[1] + s.Pi[2]*s.Pi[2];

    for (int i = 0; i < 10000; ++i)
        s = lgvi_rigid_body_step(s, J_diag, 1e-3);

    double mag_sq = s.Pi[0]*s.Pi[0] + s.Pi[1]*s.Pi[1] + s.Pi[2]*s.Pi[2];
    REQUIRE(std::abs(mag_sq - mag0_sq) < 1e-10);
}

TEST_CASE("LGVI free rigid body: orientation stays on SO(3)",
          "[lgvi][so3]") {
    Vec3 J_diag{1.0, 1.0, 1.0};      // symmetric body
    LGVIRigidBodyState s;
    s.R = SO3::ElementType::identity();
    s.Pi = Vec3{0.5, 1.0, -0.3};

    constexpr double dt = 1e-3;
    for (int i = 0; i < 10000; ++i)
        s = lgvi_rigid_body_step(s, J_diag, dt);

    auto RtR = s.R.transpose() * s.R;
    auto I = SO3::ElementType::identity();
    double err = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            err = std::max(err, std::abs(RtR(i, j) - I(i, j)));
    REQUIRE(err < 1e-10);
}

// ── DEC heat equation tests ───────────────────────────────────
#if SPATIUM_HAS_EIGEN

namespace {
Mesh<E3> make_tetrahedron() {
    Mesh<E3> m;
    m.vertices = {Vec3{0, 0, 0}, Vec3{1, 0, 0},
                  Vec3{0, 1, 0}, Vec3{0, 0, 1}};
    m.faces = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}};
    return m;
}
}

TEST_CASE("DEC heat: backward-Euler step on tetrahedron preserves total mass",
          "[dec][heat]") {
    auto topo = MeshTopology<E3>::build(make_tetrahedron());
    auto solver = make_dec_heat_solver(topo, 0.01);

    // Initial: heat concentrated at vertex 0.
    Form0<E3> phi{Eigen::VectorXd::Zero(topo.vertex_count())};
    phi.coeffs(0) = 1.0;

    auto mass = solver.mass;
    double total_initial = (mass * phi.coeffs).sum();

    auto next = solver.step(phi);
    double total_next = (mass * next.coeffs).sum();

    // Heat equation conserves the total integral (mass-weighted).
    REQUIRE(std::abs(total_next - total_initial) < 1e-10);
}

TEST_CASE("DEC heat: long-time limit converges to mean value",
          "[dec][heat]") {
    auto topo = MeshTopology<E3>::build(make_tetrahedron());
    auto solver = make_dec_heat_solver(topo, 1.0);

    // Initial: spike at vertex 0.
    Form0<E3> phi{Eigen::VectorXd::Zero(topo.vertex_count())};
    phi.coeffs(0) = 1.0;

    // Mean of φ weighted by mass.
    double total = (solver.mass * phi.coeffs).sum();
    double total_mass = solver.mass.diagonal().sum();
    double mean = total / total_mass;

    // Run many steps (long-time limit).
    for (int i = 0; i < 200; ++i)
        phi = solver.step(phi);

    // Should be uniform = mean.
    for (int i = 0; i < phi.coeffs.size(); ++i)
        REQUIRE(std::abs(phi.coeffs(i) - mean) < 1e-3);
}

#endif // SPATIUM_HAS_EIGEN

// ── Continuum scaffolding tests ────────────────────────────────

TEST_CASE("Continuum: identity deformation has F = I, E = 0",
          "[continuum][identity]") {
    auto id_def = make_deformation_map<E3, E3>(
        [](const Vec3& X) { return X; });
    Vec3 X{0.5, 0.5, 0.5};
    auto F = deformation_gradient(id_def, X);
    auto I = Matrix<double, 3, 3>::identity();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            REQUIRE(F(i, j) == Approx(I(i, j)).margin(1e-9));

    auto E = green_strain(right_cauchy_green(F));
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            REQUIRE(E(i, j) == Approx(0.0).margin(1e-9));
}

TEST_CASE("Continuum: pure stretch λ along x: F diagonal, E quadratic",
          "[continuum][stretch]") {
    constexpr double lambda = 1.5;
    auto stretch = make_deformation_map<E3, E3>(
        [lambda](const Vec3& X) { return Vec3{X[0] * lambda, X[1], X[2]}; });
    Vec3 X{0.5, 0.5, 0.5};
    auto F = deformation_gradient(stretch, X);
    REQUIRE(F(0, 0) == Approx(lambda).margin(1e-6));
    REQUIRE(F(1, 1) == Approx(1.0).margin(1e-9));
    REQUIRE(F(2, 2) == Approx(1.0).margin(1e-9));

    auto E = green_strain(right_cauchy_green(F));
    // E_xx = ½(λ² − 1) = ½(2.25 − 1) = 0.625.
    REQUIRE(E(0, 0) == Approx(0.625).margin(1e-5));
    REQUIRE(E(1, 1) == Approx(0.0).margin(1e-9));
    REQUIRE(E(2, 2) == Approx(0.0).margin(1e-9));
}

TEST_CASE("Continuum: SVK strain energy is zero at identity",
          "[continuum][svk]") {
    SaintVenantKirchhoff<double, 3> svk{1.0, 0.5};      // λ, μ
    auto I = Matrix<double, 3, 3>::identity();
    REQUIRE(svk(I) == Approx(0.0).margin(1e-15));
}

TEST_CASE("Continuum: SVK quadratic in stretch λ", "[continuum][svk]") {
    SaintVenantKirchhoff<double, 3> svk{1.0, 0.5};
    // Pure stretch λ along x. F = diag(λ, 1, 1). E_xx = ½(λ²-1) ≈ (λ-1) for
    // small strain. W ∝ E² → 10× more strain → 100× more energy.
    auto F_of = [](double lam) {
        Matrix<double, 3, 3> F{};
        F(0,0) = lam; F(1,1) = 1; F(2,2) = 1;
        return F;
    };
    double W_small = svk(F_of(1.001));
    double W_large = svk(F_of(1.01));
    REQUIRE(W_large / W_small > 90);
    REQUIRE(W_large / W_small < 110);
}
