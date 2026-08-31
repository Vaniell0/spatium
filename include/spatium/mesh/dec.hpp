#pragma once

// Discrete Exterior Calculus (DEC) primitives on triangle meshes.
//
// Wraps the operators in `mesh/differential.hpp` (cotangent Laplacian, lumped
// mass, face gradients, integrated divergence) into the canonical DEC algebra:
//
//   d_k     : Form<k>      → Form<k+1>      (exterior derivative)
//   star_k  : Form<k>      → Form<n-k>      (Hodge dual; ⋆₀ = mass, ⋆₁ = cot)
//   delta_k : Form<k>      → Form<k-1>      (codifferential = ⋆⁻¹ d^T ⋆)
//   Δ       = δd + dδ                       (Laplace-Beltrami)
//
// On a triangulated 2-manifold:
//   Form<0> ↔ vertex-valued field (n_v components)
//   Form<1> ↔ edge-valued one-form (n_e components)
//   Form<2> ↔ face-valued two-form (n_f components)
//
// Validation: `laplace_beltrami(topo)` reproduces the existing
// `differential::build_laplacian(topo)` (sign-aware) — DEC is a typed view
// over the same numerics, not a reimplementation.
//
// Eigen is required (Form storage = `Eigen::VectorXd`, operators = `Eigen::SparseMatrix`).
// References:
//   - Hirani, A. N. *Discrete Exterior Calculus* (PhD, Caltech 2003).
//   - Crane, K. *Discrete Differential Geometry: An Applied Introduction*.
//   - de Goes, F. et al. CAGD 2024 (polygonal generalisation).

#include <spatium/_export_macro.hpp>

#if defined(SPATIUM_HAS_EIGEN) && SPATIUM_HAS_EIGEN

#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/mesh/topology.hpp>
#  include <spatium/mesh/differential.hpp>
#  include <Eigen/Sparse>
#  include <Eigen/SparseCholesky>
#  include <cmath>
#  include <memory>
#endif

SPATIUM_EXPORT namespace spatium::mesh {

// ── Strongly-typed differential form ──────────────────────────
// `K` is the form degree (0, 1, or 2 on a triangle mesh). The wrapper carries
// the ambient topology so we can bound-check + chain operators ergonomically.

template<int K, Surface S>
struct DiscreteForm {
    using ScalarT = typename S::ScalarType;
    static constexpr int degree = K;

    Eigen::Matrix<ScalarT, Eigen::Dynamic, 1> coeffs;

    DiscreteForm() = default;
    explicit DiscreteForm(Eigen::Matrix<ScalarT, Eigen::Dynamic, 1> c)
        : coeffs(std::move(c)) {}

    [[nodiscard]] Eigen::Index size() const { return coeffs.size(); }

    DiscreteForm operator+(const DiscreteForm& rhs) const {
        return DiscreteForm{coeffs + rhs.coeffs};
    }
    DiscreteForm operator-(const DiscreteForm& rhs) const {
        return DiscreteForm{coeffs - rhs.coeffs};
    }
    DiscreteForm operator*(ScalarT s) const {
        return DiscreteForm{coeffs * s};
    }
};

template<Surface S> using Form0 = DiscreteForm<0, S>;
template<Surface S> using Form1 = DiscreteForm<1, S>;
template<Surface S> using Form2 = DiscreteForm<2, S>;

// ── Exterior derivative d_0: Form<0> → Form<1> ────────────────
// Edge-vertex incidence matrix. Sign convention: each edge with sorted
// vertices (v0, v1, v0 < v1) contributes -1 at v0 and +1 at v1, so for a
// scalar field φ the row gives φ(v1) - φ(v0) — the discrete gradient
// integrated along the edge.

template<Surface S>
Eigen::SparseMatrix<typename S::ScalarType> exterior_derivative_0(
    const MeshTopology<S>& topo)
{
    using T = typename S::ScalarType;
    const auto nv = topo.vertex_count();
    const auto ne = topo.edge_count();
    Eigen::SparseMatrix<T> D(static_cast<Eigen::Index>(ne),
                             static_cast<Eigen::Index>(nv));
    std::vector<Eigen::Triplet<T>> trips;
    trips.reserve(ne * 2);
    for (std::uint32_t e = 0; e < ne; ++e) {
        const auto& E = topo.edge(e);
        trips.emplace_back(static_cast<Eigen::Index>(e),
                           static_cast<Eigen::Index>(E.v0), T{-1});
        trips.emplace_back(static_cast<Eigen::Index>(e),
                           static_cast<Eigen::Index>(E.v1), T{1});
    }
    D.setFromTriplets(trips.begin(), trips.end());
    return D;
}

// ── Exterior derivative d_1: Form<1> → Form<2> ────────────────
// Face-edge incidence with orientation (face vertices ordered (a,b,c) imply
// boundary edges ab, bc, ca; sign +1 if the topology edge stores them in the
// same order, -1 if reversed).

template<Surface S>
Eigen::SparseMatrix<typename S::ScalarType> exterior_derivative_1(
    const MeshTopology<S>& topo)
{
    using T = typename S::ScalarType;
    const auto nf = static_cast<Eigen::Index>(topo.mesh().faces.size());
    const auto ne = static_cast<Eigen::Index>(topo.edge_count());
    Eigen::SparseMatrix<T> D(nf, ne);
    std::vector<Eigen::Triplet<T>> trips;
    trips.reserve(static_cast<std::size_t>(nf) * 3);
    for (std::uint32_t f = 0; f < topo.mesh().faces.size(); ++f) {
        auto [a, b, c] = topo.mesh().faces[f];
        std::array<std::pair<std::uint32_t, std::uint32_t>, 3> face_oriented{{
            {a, b}, {b, c}, {c, a}
        }};
        const auto& fe = topo.face_edges(f);
        for (int k = 0; k < 3; ++k) {
            const auto& E = topo.edge(fe[k]);
            const auto& fo = face_oriented[k];
            T sign = (E.v0 == fo.first && E.v1 == fo.second) ? T{1} : T{-1};
            trips.emplace_back(static_cast<Eigen::Index>(f),
                               static_cast<Eigen::Index>(fe[k]), sign);
        }
    }
    D.setFromTriplets(trips.begin(), trips.end());
    return D;
}

// ── Hodge star ⋆₀: Form<0> → Form<2>  (vertex masses, lumped) ─
// Returns the diagonal lumped mass matrix from differential.hpp — already
// the standard ⋆₀ in the primal-vertex / dual-cell DEC.

template<Surface S>
Eigen::SparseMatrix<typename S::ScalarType> hodge_star_0(
    const MeshTopology<S>& topo)
{
    return build_mass_matrix(topo);
}

// ── Hodge star ⋆₁: Form<1> → Form<1>  (cotan-weighted edge stars) ─
// ⋆₁_ee = ½(cot α + cot β), where α and β are the angles opposite edge e in
// its two adjacent triangles. Boundary edges (one face) keep just one term.

template<Surface S>
Eigen::SparseMatrix<typename S::ScalarType> hodge_star_1(
    const MeshTopology<S>& topo)
{
    using T = typename S::ScalarType;
    const auto ne = topo.edge_count();
    Eigen::SparseMatrix<T> H(static_cast<Eigen::Index>(ne),
                             static_cast<Eigen::Index>(ne));
    std::vector<Eigen::Triplet<T>> trips;
    trips.reserve(ne);

    const auto& mesh = topo.mesh();
    for (std::uint32_t e = 0; e < ne; ++e) {
        const auto& E = topo.edge(e);
        T weight{0};
        for (auto fi : E.adj_faces) {
            if (fi == no_face) continue;
            auto [a, b, c] = mesh.faces[fi];
            // The third vertex of the triangle — the one opposite the edge —
            // is `a + b + c − E.v0 − E.v1`. The two on-edge vertices cancel,
            // leaving the off-edge index. Avoids a 3-way `if` chain.
            std::uint32_t opp = a + b + c - E.v0 - E.v1;
            const auto& pi = mesh.vertices[E.v0];
            const auto& pj = mesh.vertices[E.v1];
            const auto& pk = mesh.vertices[opp];
            weight += T{0.5} * cotangent_weight<T, S::PointType::size>(pi, pj, pk);
        }
        trips.emplace_back(static_cast<Eigen::Index>(e),
                           static_cast<Eigen::Index>(e), weight);
    }
    H.setFromTriplets(trips.begin(), trips.end());
    return H;
}

// ── Hodge star ⋆₂: Form<2> → Form<0>  (1 / face area, diagonal) ─
// On the dual side, a primal face maps to a dual vertex; ⋆₂ scales by the
// reciprocal of the primal area so that integrating ⋆₂ω over the face
// recovers ω evaluated at the dual vertex.

template<Surface S>
Eigen::SparseMatrix<typename S::ScalarType> hodge_star_2(
    const MeshTopology<S>& topo)
{
    using T = typename S::ScalarType;
    const auto nf = static_cast<Eigen::Index>(topo.mesh().faces.size());
    Eigen::SparseMatrix<T> H(nf, nf);
    std::vector<Eigen::Triplet<T>> trips;
    trips.reserve(static_cast<std::size_t>(nf));

    const auto& mesh = topo.mesh();
    for (std::uint32_t f = 0; f < mesh.faces.size(); ++f) {
        auto [a, b, c] = mesh.faces[f];
        T area;
        if constexpr (S::PointType::size >= 3)
            area = face_area_3d(mesh.vertices[a], mesh.vertices[b], mesh.vertices[c]);
        else
            area = face_area_2d(mesh.vertices[a], mesh.vertices[b], mesh.vertices[c]);
        T value = (area > T{0}) ? (T{1} / area) : T{0};
        trips.emplace_back(static_cast<Eigen::Index>(f),
                           static_cast<Eigen::Index>(f), value);
    }
    H.setFromTriplets(trips.begin(), trips.end());
    return H;
}

// ── Laplace-Beltrami via DEC: Δ = ⋆₀⁻¹ d₀ᵀ ⋆₁ d₀ ───────────────
// Returns the symmetric weak Laplacian (`d₀ᵀ ⋆₁ d₀`) plus its mass matrix
// so consumers can solve generalised eigenproblems / heat flow.
// Sign convention: positive semi-definite. The historical
// `build_laplacian` returns the negative of this convention; both are
// reconcilable by flipping a sign.

template<Surface S>
struct LaplaceBeltrami {
    Eigen::SparseMatrix<typename S::ScalarType> weak;   // d₀ᵀ ⋆₁ d₀  (PSD)
    Eigen::SparseMatrix<typename S::ScalarType> mass;   // ⋆₀ (diagonal)
};

template<Surface S>
LaplaceBeltrami<S> laplace_beltrami_dec(const MeshTopology<S>& topo) {
    auto d0 = exterior_derivative_0(topo);
    auto s1 = hodge_star_1(topo);
    return { d0.transpose() * s1 * d0, hodge_star_0(topo) };
}

// ── Closure check: d₁ ∘ d₀ ≡ 0 (de Rham complex) ──────────────
// Returns the maximum absolute coefficient of the product matrix; should be
// ≤ epsilon for any valid simplicial complex. Useful in debug.

template<Surface S>
typename S::ScalarType d1_d0_closure_residual(const MeshTopology<S>& topo) {
    auto d0 = exterior_derivative_0(topo);
    auto d1 = exterior_derivative_1(topo);
    Eigen::SparseMatrix<typename S::ScalarType> dd = d1 * d0;
    return Eigen::Matrix<typename S::ScalarType, Eigen::Dynamic, 1>{
        Eigen::Map<const Eigen::Matrix<typename S::ScalarType, Eigen::Dynamic, 1>>(
            dd.valuePtr(), dd.nonZeros())
    }.cwiseAbs().maxCoeff();
}

// ── Heat equation via DEC (backward Euler) ────────────────────
// Solves ∂φ/∂t = -Δ φ on the mesh implicitly:
//
//     (M + dt · L) · φ_{n+1} = M · φ_n
//
// where M = ⋆₀ (lumped mass) and L = d₀ᵀ ⋆₁ d₀ (weak Laplacian, PSD).
// Backward Euler is unconditionally stable for any dt > 0 — no CFL bound.
// One sparse Cholesky factorisation up front, cheap solves per step.

template<Surface S>
struct DecHeatSolver {
    using ScalarT = typename S::ScalarType;
    using SparseM = Eigen::SparseMatrix<ScalarT>;
    using LDLT    = Eigen::SimplicialLDLT<SparseM>;

    SparseM mass;
    std::unique_ptr<LDLT> solver;          // Eigen solvers are non-copyable

    // Take one backward-Euler step: φ ← (M + dt L)⁻¹ · M · φ.
    Form0<S> step(const Form0<S>& phi) const {
        Eigen::Matrix<ScalarT, Eigen::Dynamic, 1> rhs = mass * phi.coeffs;
        Eigen::Matrix<ScalarT, Eigen::Dynamic, 1> next = solver->solve(rhs);
        return Form0<S>{next};
    }
};

template<Surface S>
DecHeatSolver<S> make_dec_heat_solver(const MeshTopology<S>& topo,
                                      typename S::ScalarType dt)
{
    using SparseM = typename DecHeatSolver<S>::SparseM;
    // Use `build_laplacian` (negative semi-definite cotan operator) +
    // `build_mass_matrix` directly — `laplace_beltrami_dec` would
    // assemble the same operator via d₀ᵀ ⋆₁ d₀, paying an extra
    // sparse multiply. That DEC-typed path stays available as the
    // canonical reference checked by test_dec.cpp.
    SparseM mass  = build_mass_matrix(topo);
    SparseM L_neg = build_laplacian(topo);            // negative semi-definite
    SparseM A     = mass - L_neg * dt;                // = mass + L_pos · dt
    DecHeatSolver<S> h;
    h.mass = mass;
    h.solver = std::make_unique<typename DecHeatSolver<S>::LDLT>();
    h.solver->compute(A);
    return h;
}

} // namespace spatium::mesh

#endif  // SPATIUM_HAS_EIGEN
