#pragma once

// Discrete differential operators on triangle meshes.
// Requires Eigen (SPATIUM_HAS_EIGEN) for sparse matrices.

#include <spatium/_export_macro.hpp>

#if defined(SPATIUM_HAS_EIGEN) && SPATIUM_HAS_EIGEN

#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/mesh/topology.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <Eigen/Sparse>
#  include <algorithm>
#  include <cmath>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::mesh {

using spatium::Vec;
using spatium::Scalar;
using spatium::epsilon;

// Cotangent of the angle at vertex pk in triangle (pi, pj, pk).
// Returns cot(angle_k) = dot(pi-pk, pj-pk) / |cross(pi-pk, pj-pk)|.
// Clamped to avoid infinity on degenerate triangles.
template<Scalar T, std::size_t N>
T cotangent_weight(const Vec<T, N>& pi, const Vec<T, N>& pj, const Vec<T, N>& pk) {
    auto u = pi - pk;
    auto v = pj - pk;
    auto cross_norm = u.cross(v).norm();
    auto dot_val = u.dot(v);
    constexpr T max_cot = T{1e6};
    if (cross_norm < epsilon<T>())
        return dot_val >= T{0} ? max_cot : -max_cot;
    return std::clamp(dot_val / cross_norm, -max_cot, max_cot);
}

// Build the cotangent Laplacian matrix L (n x n, symmetric, negative semi-definite).
// L(i,j) = -0.5 * sum of cotangent weights from adjacent triangles.
// L(i,i) = -sum_j L(i,j).
template<Surface S>
Eigen::SparseMatrix<typename S::ScalarType> build_laplacian(const MeshTopology<S>& topo) {
    using T = typename S::ScalarType;
    auto n = static_cast<int>(topo.vertex_count());
    auto& m = topo.mesh();

    using Triplet = Eigen::Triplet<T>;
    std::vector<Triplet> trips;
    trips.reserve(m.faces.size() * 6);

    for (std::size_t fi = 0; fi < m.faces.size(); ++fi) {
        auto [a, b, c] = m.faces[fi];
        const auto& pa = m.vertices[a];
        const auto& pb = m.vertices[b];
        const auto& pc = m.vertices[c];

        // Cotangent weight for edge (b,c) opposite vertex a
        T w_a = cotangent_weight(pb, pc, pa);
        // Edge (a,c) opposite vertex b
        T w_b = cotangent_weight(pa, pc, pb);
        // Edge (a,b) opposite vertex c
        T w_c = cotangent_weight(pa, pb, pc);

        // Off-diagonal: L(i,j) += -0.5 * cot(opposite angle)
        auto add_entry = [&](uint32_t i, uint32_t j, T w) {
            T val = T{-0.5} * w;
            trips.emplace_back(static_cast<int>(i), static_cast<int>(j), val);
            trips.emplace_back(static_cast<int>(j), static_cast<int>(i), val);
        };

        add_entry(b, c, w_a);
        add_entry(a, c, w_b);
        add_entry(a, b, w_c);
    }

    Eigen::SparseMatrix<T> L(n, n);
    L.setFromTriplets(trips.begin(), trips.end());

    // Diagonal: L(i,i) = -sum_j L(i,j)
    for (int i = 0; i < n; ++i) {
        T row_sum = T{0};
        for (typename Eigen::SparseMatrix<T>::InnerIterator it(L, i); it; ++it) {
            if (it.row() != i)
                row_sum += it.value();
        }
        L.coeffRef(i, i) = -row_sum;
    }

    return L;
}

// Build the lumped mass matrix A (diagonal, n x n).
// A(i,i) = (1/3) * sum of areas of triangles incident to vertex i.
template<Surface S>
Eigen::SparseMatrix<typename S::ScalarType> build_mass_matrix(const MeshTopology<S>& topo) {
    using T = typename S::ScalarType;
    auto n = static_cast<int>(topo.vertex_count());
    auto& m = topo.mesh();

    Eigen::Matrix<T, Eigen::Dynamic, 1> diag =
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);

    for (std::size_t fi = 0; fi < m.faces.size(); ++fi) {
        auto [a, b, c] = m.faces[fi];
        T area = face_area_3d(m.vertices[a], m.vertices[b], m.vertices[c]);
        T third = area / T{3};
        diag(static_cast<int>(a)) += third;
        diag(static_cast<int>(b)) += third;
        diag(static_cast<int>(c)) += third;
    }

    using Triplet = Eigen::Triplet<T>;
    std::vector<Triplet> trips;
    trips.reserve(n);
    for (int i = 0; i < n; ++i)
        trips.emplace_back(i, i, std::max(diag(i), epsilon<T>()));

    Eigen::SparseMatrix<T> A(n, n);
    A.setFromTriplets(trips.begin(), trips.end());
    return A;
}

// Compute per-face gradient of a scalar field u defined at vertices.
// For triangle (i,j,k) with area A_f and normal n_f:
//   grad_f(u) = (1/(2*A_f)) * sum_l u_l * (e_opposite_l x n_f)
template<Surface S>
std::vector<Vec<typename S::ScalarType, 3>> face_gradients(
    const MeshTopology<S>& topo,
    const Eigen::Matrix<typename S::ScalarType, Eigen::Dynamic, 1>& u)
{
    using T = typename S::ScalarType;
    auto& m = topo.mesh();
    std::vector<Vec<T, 3>> grads(m.faces.size());

    for (std::size_t fi = 0; fi < m.faces.size(); ++fi) {
        auto [a, b, c] = m.faces[fi];
        const auto& pa = m.vertices[a];
        const auto& pb = m.vertices[b];
        const auto& pc = m.vertices[c];

        auto e_ab = pb - pa;  // edge opposite c... no, e_ab = b-a
        auto e_bc = pc - pb;
        auto e_ca = pa - pc;
        auto normal = e_ab.cross(e_bc);  // not normalized, magnitude = 2*area
        T twice_area = normal.norm();

        if (twice_area < epsilon<T>()) {
            grads[fi] = Vec<T, 3>{};
            continue;
        }

        auto n = normal / twice_area;  // unit normal

        // grad = (1/(2A)) * (u_a * (n x e_bc) + u_b * (n x e_ca) + u_c * (n x e_ab))
        auto ga = n.cross(e_bc) * u(static_cast<int>(a));
        auto gb = n.cross(e_ca) * u(static_cast<int>(b));
        auto gc = n.cross(e_ab) * u(static_cast<int>(c));

        grads[fi] = Vec<T, 3>{(ga + gb + gc) / twice_area};
    }

    return grads;
}

// Compute integrated divergence from per-face vector field X.
// div_i = (1/2) * sum over incident faces f of:
//   cot(angle_j) * (p_i - p_k) . X_f  +  cot(angle_k) * (p_i - p_j) . X_f
// where j,k are the other two vertices of face f.
template<Surface S>
Eigen::Matrix<typename S::ScalarType, Eigen::Dynamic, 1> integrated_divergence(
    const MeshTopology<S>& topo,
    const std::vector<Vec<typename S::ScalarType, 3>>& X)
{
    using T = typename S::ScalarType;
    auto n = static_cast<int>(topo.vertex_count());
    auto& m = topo.mesh();

    Eigen::Matrix<T, Eigen::Dynamic, 1> div =
        Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);

    for (std::size_t fi = 0; fi < m.faces.size(); ++fi) {
        auto [a, b, c] = m.faces[fi];
        const auto& pa = m.vertices[a];
        const auto& pb = m.vertices[b];
        const auto& pc = m.vertices[c];

        T cot_a = cotangent_weight(pb, pc, pa);
        T cot_b = cotangent_weight(pa, pc, pb);
        T cot_c = cotangent_weight(pa, pb, pc);

        const auto& xf = X[fi];

        // For vertex a: cot_b * (pa - pc) . X + cot_c * (pa - pb) . X
        div(static_cast<int>(a)) += T{0.5} * (cot_b * (pa - pc).dot(xf) + cot_c * (pa - pb).dot(xf));
        // For vertex b: cot_a * (pb - pc) . X + cot_c * (pb - pa) . X
        div(static_cast<int>(b)) += T{0.5} * (cot_a * (pb - pc).dot(xf) + cot_c * (pb - pa).dot(xf));
        // For vertex c: cot_a * (pc - pb) . X + cot_b * (pc - pa) . X
        div(static_cast<int>(c)) += T{0.5} * (cot_a * (pc - pb).dot(xf) + cot_b * (pc - pa).dot(xf));
    }

    return div;
}

// Mean edge length of the mesh.
template<Surface S>
typename S::ScalarType mean_edge_length(const MeshTopology<S>& topo) {
    using T = typename S::ScalarType;
    auto& m = topo.mesh();
    T total{0};
    auto ne = topo.edge_count();
    for (uint32_t i = 0; i < ne; ++i) {
        auto& e = topo.edge(i);
        total += (m.vertices[e.v0] - m.vertices[e.v1]).norm();
    }
    return ne > 0 ? total / static_cast<T>(ne) : T{1};
}

} // namespace spatium::mesh

#endif // SPATIUM_HAS_EIGEN
