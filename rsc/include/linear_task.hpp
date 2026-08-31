#pragma once

// Domain 6 (general linear-solve dispatch), expressed as configuration
// over comparison_task.hpp's generic ComparisonTaskGenerator -- sixth
// instance, same shape as precision/geodesic/rootfind/ode/mesh.
//
// Ground truth is generated, not hand-labeled: candidates are checked
// cheapest-first (linalg_jacobi, index 0) against linalg_direct (last,
// the reference -- Gaussian elimination is exact regardless of
// conditioning, so its distance to itself is always 0, same shape
// precision_task.hpp's Real50 reference has). Jacobi wins whenever its
// fixed 20-iteration budget (linear_ops.hpp) actually converges; it
// measurably fails to when the matrix isn't diagonally dominant enough
// -- a genuine two-regime split, MEASURED with a standalone probe before
// writing this file (not assumed): at N=4, diag=6.0, a 50/50 mix of
// off-diagonal scale in [0.5,2.0] ("safe") vs. [4.0,7.0] ("risky") gives
// a real ~36% needs-direct fraction under the exact comparison rule
// (||jacobi - direct|| <= tolerance) this file uses.
//
// feature = log(diagonal_dominance_ratio(A) + eps) - mean, same log-and-
// center discipline every prior domain with an always-positive raw
// feature needed (geodesic's vertex_count, rootfind's |f'(x0)|, ode's
// char_freq*dt) -- diagonal_dominance_ratio() is >= 0 by construction
// (a ratio of absolute values), so this was applied up front rather than
// rediscovered as a dead-ReLU bug. feature_mean measured by actually
// sampling the real problem distribution at construction time (matches
// rootfind_task.hpp/ode_task.hpp's own warmup-pass approach).

#include <comparison_task.hpp>
#include <linear_ops.hpp>
#include <cmath>

namespace rsc {

struct LinearProblem {
    spatium::Matrix<double, kLinearN, kLinearN> A;
    spatium::Vec<double, kLinearN> b;
};

using LinearOutput = spatium::Vec<double, kLinearN>;
using LinearTaskGenerator = ComparisonTaskGenerator<LinearProblem, LinearOutput>;

namespace detail {

inline constexpr double kLinearDiag = 6.0;

inline double sample_off_scale(std::mt19937_64& rng) {
    std::bernoulli_distribution risky(0.5);
    if (risky(rng)) {
        std::uniform_real_distribution<double> d(4.0, 7.0);  // weakly/non dominant
        return d(rng);
    }
    std::uniform_real_distribution<double> d(0.5, 2.0);  // strongly dominant
    return d(rng);
}

inline LinearProblem sample_linear_problem(std::mt19937_64& rng) {
    double off_scale = sample_off_scale(rng);
    std::uniform_real_distribution<double> off(-off_scale, off_scale);
    std::uniform_real_distribution<double> b_dist(-5.0, 5.0);

    spatium::Matrix<double, kLinearN, kLinearN> A;
    for (std::size_t i = 0; i < kLinearN; ++i)
        for (std::size_t j = 0; j < kLinearN; ++j)
            A(i, j) = (i == j) ? kLinearDiag : off(rng);
    spatium::Vec<double, kLinearN> b;
    for (std::size_t i = 0; i < kLinearN; ++i) b[i] = b_dist(rng);
    return {A, b};
}

inline std::vector<double> flatten(const LinearProblem& p) {
    std::vector<double> in(kLinearN * kLinearN + kLinearN);
    for (std::size_t i = 0; i < kLinearN; ++i)
        for (std::size_t j = 0; j < kLinearN; ++j) in[i * kLinearN + j] = p.A(i, j);
    for (std::size_t i = 0; i < kLinearN; ++i) in[kLinearN * kLinearN + i] = p.b[i];
    return in;
}

} // namespace detail

inline LinearTaskGenerator build_linear_task_generator(const Registry& registry,
                                                          std::uint64_t seed = 0,
                                                          double tolerance = 1e-3) {
    auto call = [&registry](std::size_t op_index, const LinearProblem& p) -> LinearOutput {
        std::vector<double> out(kLinearN);
        registry[op_index](detail::flatten(p), out);
        return unflatten_vec(out);
    };

    std::size_t jacobi_index = registry.index_of("linalg_jacobi");
    std::size_t direct_index = registry.index_of("linalg_direct");

    std::vector<Candidate<LinearProblem, LinearOutput>> candidates{
        {"linalg_jacobi",
         [call, jacobi_index](const LinearProblem& p) { return call(jacobi_index, p); }},
        {"linalg_direct",
         [call, direct_index](const LinearProblem& p) { return call(direct_index, p); }},
    };

    auto sample_problem = [](std::mt19937_64& rng) -> LinearProblem {
        return detail::sample_linear_problem(rng);
    };

    // Measured, not guessed: draw the real problem distribution to find
    // log(diagonal_dominance_ratio)'s actual mean, same discipline as
    // rootfind_task.hpp/ode_task.hpp's own warmup passes.
    double feature_mean = 0.0;
    {
        std::mt19937_64 warmup(seed + 1);
        constexpr int kWarmupSamples = 20000;
        for (int i = 0; i < kWarmupSamples; ++i) {
            auto p = detail::sample_linear_problem(warmup);
            double ratio = spatium::diagonal_dominance_ratio(p.A);
            feature_mean += std::log(ratio + 1e-9);
        }
        feature_mean /= kWarmupSamples;
    }

    auto extract_features = [feature_mean](const LinearProblem& p) -> std::vector<double> {
        double ratio = spatium::diagonal_dominance_ratio(p.A);
        return {std::log(ratio + 1e-9) - feature_mean};
    };

    auto reference = [call, direct_index](const LinearProblem& p) { return call(direct_index, p); };

    auto distance = [](const LinearOutput& a, const LinearOutput& b) {
        return spatium::Vec<double, kLinearN>{a - b}.norm();
    };

    return LinearTaskGenerator(std::move(candidates), sample_problem, extract_features, reference,
                                distance, tolerance, seed);
}

} // namespace rsc
