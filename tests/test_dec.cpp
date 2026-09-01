// DEC primitives. Validate against existing differential.hpp.
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#if SPATIUM_HAS_EIGEN

#include <spatium/mesh/mesh.hpp>
#include <spatium/mesh/topology.hpp>
#include <spatium/mesh/dec.hpp>
#include <spatium/mesh/differential.hpp>
#include <spatium/spaces/euclidean.hpp>

using namespace spatium;
using namespace spatium::mesh;
using Catch::Approx;

namespace {
// Unit tetrahedron — small but covers all DEC operators (V=4, E=6, F=4).
Mesh<E3> make_tetrahedron() {
    Mesh<E3> m;
    m.vertices = {Vec3{0, 0, 0}, Vec3{1, 0, 0},
                  Vec3{0, 1, 0}, Vec3{0, 0, 1}};
    m.faces = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}};
    return m;
}
}

TEST_CASE("DEC: tetrahedron has 6 edges + 4 faces", "[dec]") {
    auto topo = MeshTopology<E3>::build(make_tetrahedron());
    REQUIRE(topo.vertex_count() == 4);
    REQUIRE(topo.edge_count() == 6);
    REQUIRE(topo.mesh().faces.size() == 4);
}

TEST_CASE("DEC: d_0 has signed incidence pattern", "[dec][d0]") {
    auto topo = MeshTopology<E3>::build(make_tetrahedron());
    auto d0 = exterior_derivative_0(topo);
    REQUIRE(d0.rows() == topo.edge_count());
    REQUIRE(d0.cols() == topo.vertex_count());
    // Each edge row has exactly one -1 and one +1 → row sum = 0.
    Eigen::MatrixXd dense = Eigen::MatrixXd(d0);
    for (int e = 0; e < dense.rows(); ++e) {
        REQUIRE(dense.row(e).sum() == Approx(0.0).margin(1e-12));
        // Two non-zero entries per edge.
        int nonzero = 0;
        for (int v = 0; v < dense.cols(); ++v)
            if (std::abs(dense(e, v)) > 1e-12) ++nonzero;
        REQUIRE(nonzero == 2);
    }
}

TEST_CASE("DEC: d_1 ∘ d_0 = 0 (closure)", "[dec][closure]") {
    auto topo = MeshTopology<E3>::build(make_tetrahedron());
    double residual = d1_d0_closure_residual(topo);
    REQUIRE(residual == Approx(0.0).margin(1e-12));
}

TEST_CASE("DEC: hodge_star_0 = lumped mass", "[dec][hodge]") {
    auto topo = MeshTopology<E3>::build(make_tetrahedron());
    auto s0 = hodge_star_0(topo);
    auto m  = build_mass_matrix(topo);
    REQUIRE(s0.rows() == m.rows());
    REQUIRE(s0.cols() == m.cols());
    Eigen::MatrixXd diff = Eigen::MatrixXd(s0) - Eigen::MatrixXd(m);
    REQUIRE(diff.cwiseAbs().maxCoeff() == Approx(0.0).margin(1e-15));
}

TEST_CASE("DEC: hodge_star_2 is positive on a non-degenerate mesh", "[dec][hodge]") {
    auto topo = MeshTopology<E3>::build(make_tetrahedron());
    auto s2 = hodge_star_2(topo);
    REQUIRE(s2.rows() == 4);
    for (int f = 0; f < 4; ++f) {
        double v = Eigen::MatrixXd(s2)(f, f);
        REQUIRE(v > 0.0);  // 1/area should be finite & positive
    }
}

TEST_CASE("DEC: laplace_beltrami matches build_laplacian", "[dec][laplacian]") {
    auto topo = MeshTopology<E3>::build(make_tetrahedron());

    auto L_dec = laplace_beltrami_dec(topo);
    auto L_old = build_laplacian(topo);

    // Both formulations land on the SAME positive-on-diagonal cotan Laplacian
    // (the comment in differential.hpp claiming negative-semidefinite is
    // misleading: it's PSD by construction `L(i,i) = -sum_j L(i,j)` with
    // L(i,j) ≤ 0). DEC reproduces it exactly through `d₀ᵀ ⋆₁ d₀`.
    Eigen::MatrixXd dense_dec = Eigen::MatrixXd(L_dec.weak);
    Eigen::MatrixXd dense_old = Eigen::MatrixXd(L_old);
    Eigen::MatrixXd diff = dense_dec - dense_old;
    REQUIRE(diff.cwiseAbs().maxCoeff() == Approx(0.0).margin(1e-10));
}

TEST_CASE("DEC: Form arithmetic", "[dec][form]") {
    auto topo = MeshTopology<E3>::build(make_tetrahedron());
    Form0<E3> a{Eigen::VectorXd::Constant(topo.vertex_count(), 2.0)};
    Form0<E3> b{Eigen::VectorXd::Constant(topo.vertex_count(), 3.0)};
    auto sum = a + b;
    REQUIRE(sum.coeffs.size() == 4);
    REQUIRE(sum.coeffs(0) == 5.0);
    auto scaled = a * 0.5;
    REQUIRE(scaled.coeffs(0) == 1.0);
}

#endif // SPATIUM_HAS_EIGEN
