#pragma once

// Domain 6: general linear-solve method dispatch (Jacobi vs. direct
// Gaussian elimination) -- the general-linear-algebra gap flagged early
// in the session's own survey of include/spatium/ ("no general N×N
// linear algebra, only 2x2/3x3 closed-form inverse()"), picked as the
// next domain over A4 (contact physics, still structurally blocked on
// an ipc-toolkit backend swap) 2026-08-26.
//
// New Spatium primitive first, same two-layer pattern every domain has
// needed: include/spatium/algebra/linear_solve.hpp. Fixed N=4 (same
// scale as domain 4's CircularOrbit state), flattened into the
// registry's span<const double> shape: in = [A (row-major, N*N=16), b
// (N=4)] = 20 doubles, out = x (N=4).

#include <registry.hpp>
#include <spatium/algebra/linear_solve.hpp>
#include <string>

namespace rsc {

inline constexpr std::size_t kLinearN = 4;
inline constexpr int kJacobiMaxIter = 20;  // matches rootfind_ops.hpp's Newton budget

// Row-major flat span -> Matrix<double,N,N> (Matrix's own storage is
// column-major internally; the conversion happens once here, not on
// every element access inside the solvers).
inline spatium::Matrix<double, kLinearN, kLinearN> unflatten_matrix(std::span<const double> flat) {
    spatium::Matrix<double, kLinearN, kLinearN> A;
    for (std::size_t i = 0; i < kLinearN; ++i)
        for (std::size_t j = 0; j < kLinearN; ++j) A(i, j) = flat[i * kLinearN + j];
    return A;
}

inline spatium::Vec<double, kLinearN> unflatten_vec(std::span<const double> flat) {
    spatium::Vec<double, kLinearN> v;
    for (std::size_t i = 0; i < kLinearN; ++i) v[i] = flat[i];
    return v;
}

namespace detail {

inline std::vector<std::string> linear_input_names() {
    std::vector<std::string> names;
    for (std::size_t i = 0; i < kLinearN; ++i)
        for (std::size_t j = 0; j < kLinearN; ++j)
            names.push_back("a" + std::to_string(i) + std::to_string(j));
    for (std::size_t i = 0; i < kLinearN; ++i) names.push_back("b" + std::to_string(i));
    return names;
}

inline std::vector<std::string> linear_output_names() {
    std::vector<std::string> names;
    for (std::size_t i = 0; i < kLinearN; ++i) names.push_back("x" + std::to_string(i));
    return names;
}

} // namespace detail

inline Registry build_linear_registry() {
    Registry reg;

    reg.add({.name = "linalg_jacobi",
             .tier = Tier::General,
             .in_size = kLinearN * kLinearN + kLinearN,
             .out_size = kLinearN,
             .input_names = detail::linear_input_names(),
             .output_names = detail::linear_output_names()},
            [](std::span<const double> in, std::span<double> out) {
                auto A = unflatten_matrix(in.subspan(0, kLinearN * kLinearN));
                auto b = unflatten_vec(in.subspan(kLinearN * kLinearN, kLinearN));
                auto x = spatium::solve_jacobi(A, b, kJacobiMaxIter);
                for (std::size_t i = 0; i < kLinearN; ++i) out[i] = x[i];
            });

    reg.add({.name = "linalg_direct",
             .tier = Tier::General,
             .in_size = kLinearN * kLinearN + kLinearN,
             .out_size = kLinearN,
             .input_names = detail::linear_input_names(),
             .output_names = detail::linear_output_names()},
            [](std::span<const double> in, std::span<double> out) {
                auto A = unflatten_matrix(in.subspan(0, kLinearN * kLinearN));
                auto b = unflatten_vec(in.subspan(kLinearN * kLinearN, kLinearN));
                // linear_task.hpp's sample_problem only ever offers a
                // matrix its own reference call has already solved
                // successfully once (see that file) -- a genuine singular
                // draw here would be a measure-zero event from a
                // continuous random distribution, the same trust level
                // rootfind_ops.hpp's bisection places in its caller-
                // guaranteed bracket.
                auto x = *spatium::solve_direct(A, b);
                for (std::size_t i = 0; i < kLinearN; ++i) out[i] = x[i];
            });

    return reg;
}

} // namespace rsc
