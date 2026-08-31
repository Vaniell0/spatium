#pragma once

// Heat method geodesics (Crane et al. 2013).
// O(h^2) accuracy vs O(h) for Dijkstra.
// Requires Eigen (SPATIUM_HAS_EIGEN) for sparse linear solvers.

#include <spatium/_export_macro.hpp>

#if defined(SPATIUM_HAS_EIGEN) && SPATIUM_HAS_EIGEN

#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/mesh/geodesic_types.hpp>
#  include <spatium/mesh/differential.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <Eigen/SparseCholesky>
#  include <cassert>
#  include <cmath>
#  include <limits>
#  include <memory>
#  include <span>
#endif

SPATIUM_EXPORT namespace spatium::mesh {

using spatium::Vec;
using spatium::epsilon;

// Pre-factored heat method solver.
// Build once (expensive: two Cholesky factorizations), query many times (cheap).
// SimplicialLDLT is non-movable, so HeatSolver is heap-allocated via build().
template<Surface S>
class HeatSolver {
public:
    using T = typename S::ScalarType;
    using Solver = Eigen::SimplicialLDLT<Eigen::SparseMatrix<T>>;

    static std::unique_ptr<HeatSolver> build(const MeshTopology<S>& topo, const S& /*space*/) {
        auto solver = std::unique_ptr<HeatSolver>(new HeatSolver);
        solver->topo_ = &topo;

        auto n = static_cast<int>(topo.vertex_count());

        // Time step: t = h^2
        T h = mean_edge_length(topo);
        solver->t_ = h * h;

        // Build operators
        solver->L_ = build_laplacian(topo);
        auto A = build_mass_matrix(topo);

        // Heat solve: (A + t*L) u = delta_s
        // L is positive semi-definite (negative off-diagonal, positive diagonal).
        Eigen::SparseMatrix<T> heat_op = A + solver->t_ * solver->L_;
        solver->heat_solver_ = std::make_unique<Solver>();
        solver->heat_solver_->compute(heat_op);
        assert(solver->heat_solver_->info() == Eigen::Success);

        // Poisson solve: L * phi = div(X).
        // L is positive semi-definite. Pin vertex 0 to fix gauge.
        Eigen::SparseMatrix<T> L_pinned = solver->L_;
        for (typename Eigen::SparseMatrix<T>::InnerIterator it(L_pinned, 0); it; ++it)
            it.valueRef() = T{0};
        for (int k = 0; k < L_pinned.outerSize(); ++k) {
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(L_pinned, k); it; ++it) {
                if (it.row() == 0) it.valueRef() = T{0};
            }
        }
        L_pinned.coeffRef(0, 0) = T{1};

        solver->poisson_solver_ = std::make_unique<Solver>();
        solver->poisson_solver_->compute(L_pinned);
        assert(solver->poisson_solver_->info() == Eigen::Success);

        solver->mass_diag_.resize(n);
        for (int i = 0; i < n; ++i)
            solver->mass_diag_(i) = A.coeff(i, i);

        return solver;
    }

    DistanceField<S> distances(uint32_t source) const {
        return distances(std::span<const uint32_t>(&source, 1));
    }

    DistanceField<S> distances(std::span<const uint32_t> sources) const {
        auto n = static_cast<int>(topo_->vertex_count());

        // Step 1: RHS for heat equation (delta at sources, weighted by mass)
        Eigen::Matrix<T, Eigen::Dynamic, 1> rhs =
            Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(n);
        for (auto s : sources)
            rhs(static_cast<int>(s)) = mass_diag_(static_cast<int>(s));

        // Step 2: Solve heat diffusion
        Eigen::Matrix<T, Eigen::Dynamic, 1> u = heat_solver_->solve(rhs);

        // Step 3: Compute normalized negative gradient per face
        auto grads = face_gradients(*topo_, u);
        for (auto& g : grads) {
            T gn = g.norm();
            if (gn > epsilon<T>())
                g = Vec<T, 3>{g * (T{-1} / gn)};
            else
                g = Vec<T, 3>{};
        }

        // Step 4: Compute divergence
        auto div = integrated_divergence(*topo_, grads);
        div(0) = T{0};  // pinned vertex

        // Step 5: Solve Poisson for distance
        Eigen::Matrix<T, Eigen::Dynamic, 1> phi = poisson_solver_->solve(div);

        // Step 6: Shift so that min source distance = 0
        T min_phi = std::numeric_limits<T>::max();
        for (auto s : sources)
            min_phi = std::min(min_phi, phi(static_cast<int>(s)));
        for (int i = 0; i < n; ++i)
            phi(i) -= min_phi;

        // Package result
        DistanceField<S> field;
        field.distances.resize(n);
        field.predecessors.assign(n, no_vertex);
        for (int i = 0; i < n; ++i)
            field.distances[i] = std::abs(phi(i));

        return field;
    }

private:
    HeatSolver() = default;

    const MeshTopology<S>* topo_ = nullptr;
    Eigen::SparseMatrix<T> L_;
    std::unique_ptr<Solver> heat_solver_;
    std::unique_ptr<Solver> poisson_solver_;
    Eigen::Matrix<T, Eigen::Dynamic, 1> mass_diag_;
    T t_{};
};

// One-shot convenience function.
template<Surface S>
DistanceField<S> heat_geodesic_distances(
    const MeshTopology<S>& topo,
    const S& space,
    uint32_t source)
{
    auto solver = HeatSolver<S>::build(topo, space);
    return solver->distances(source);
}

template<Surface S>
DistanceField<S> heat_geodesic_distances(
    const MeshTopology<S>& topo,
    const S& space,
    std::span<const uint32_t> sources)
{
    auto solver = HeatSolver<S>::build(topo, space);
    return solver->distances(sources);
}

} // namespace spatium::mesh

#endif // SPATIUM_HAS_EIGEN
