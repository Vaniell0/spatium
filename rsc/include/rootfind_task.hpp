#pragma once

// Domain 3 (root-finding dispatch), expressed as configuration over
// comparison_task.hpp's generic ComparisonTaskGenerator -- same shape as
// precision_task.hpp/geodesic_task.hpp, third instance now, not a fourth
// bespoke class.
//
// Ground truth is generated, not hand-labeled, same principle as domains
// 1-2: candidates are checked cheapest-first (newton_root, index 0) against
// bisection_root (last, treated as the reference-quality candidate since a
// generous fixed bracket + 50 iterations converges it to ~1e-14 regardless
// of a). Newton wins whenever it actually converges within its 20-iteration
// budget; it fails to when x0 starts near the f'=0 inflection (see
// rootfind_ops.hpp) -- a genuine two-regime split, not assumed in advance.
//
// x0 is sampled from two mixed regimes -- "safe" (away from the
// inflection) and "risky" (close to it) -- same regime-mixing shape
// precision_task.hpp uses for root separation. feature = log(|f'(x0)|) -
// mean, both the log-scale and the centering applied up front rather than
// discovered as a bug later: geodesic_task.hpp hit exactly this failure
// mode (a raw always-positive single feature permanently dead-ReLUs about
// half the hidden units at init) with vertex_count, and |f'(x0)| = 3*x0^2
// is the same shape (always >= 0). feature_mean is measured by actually
// sampling x0's own distribution at construction time (matches
// geodesic_task.hpp's precomputed feature_mean approach), not guessed
// analytically.

#include <comparison_task.hpp>
#include <rootfind_ops.hpp>
#include <cmath>

namespace rsc {

struct RootProblem {
    double a;          // target: solve x^3 - a = 0
    double x0;         // Newton's initial guess
    double lo, hi;      // bisection bracket, always valid by construction
};

using RootOutput = double;
using RootTaskGenerator = ComparisonTaskGenerator<RootProblem, RootOutput>;

namespace detail {

inline double sample_x0(std::mt19937_64& rng) {
    std::bernoulli_distribution risky(0.5);
    if (risky(rng)) {
        std::uniform_real_distribution<double> near_zero(-0.15, 0.15);
        return near_zero(rng);
    }
    std::uniform_real_distribution<double> away(0.5, 3.0);
    std::bernoulli_distribution sign(0.5);
    double x0 = away(rng);
    return sign(rng) ? x0 : -x0;
}

} // namespace detail

inline RootTaskGenerator build_rootfind_task_generator(const Registry& registry,
                                                         std::uint64_t seed = 0,
                                                         double tolerance = 1e-6) {
    auto call = [&registry](std::size_t op_index, std::vector<double> in) {
        std::vector<double> out(1);
        registry[op_index](in, out);
        return out[0];
    };

    std::size_t newton_index = registry.index_of("newton_root");
    std::size_t bisection_index = registry.index_of("bisection_root");

    std::vector<Candidate<RootProblem, RootOutput>> candidates{
        {"newton_root",
         [call, newton_index](const RootProblem& p) { return call(newton_index, {p.x0, p.a}); }},
        {"bisection_root",
         [call, bisection_index](const RootProblem& p) {
             return call(bisection_index, {p.lo, p.hi, p.a});
         }},
    };

    auto sample_problem = [](std::mt19937_64& rng) -> RootProblem {
        std::uniform_real_distribution<double> a_dist(-8.0, 8.0);
        double a = a_dist(rng);
        double x0 = detail::sample_x0(rng);
        return {a, x0, -10.0, 10.0};
    };

    // Measured, not guessed: draw x0's own marginal distribution to find
    // the real mean of log(|f'(x0)|) -- f'(x)=3x^2 doesn't depend on a, so
    // a's distribution is irrelevant here.
    double feature_mean = 0.0;
    {
        std::mt19937_64 warmup(seed + 1);
        constexpr int kWarmupSamples = 20000;
        for (int i = 0; i < kWarmupSamples; ++i) {
            double x0 = detail::sample_x0(warmup);
            double fprime = 3.0 * x0 * x0;
            feature_mean += std::log(fprime + 1e-9);
        }
        feature_mean /= kWarmupSamples;
    }

    auto extract_features = [feature_mean](const RootProblem& p) -> std::vector<double> {
        double fprime = 3.0 * p.x0 * p.x0;
        return {std::log(fprime + 1e-9) - feature_mean};
    };

    auto reference = [](const RootProblem& p) {
        using std::cbrt;
        return cbrt(p.a);
    };

    auto distance = [](const RootOutput& a, const RootOutput& b) { return std::abs(a - b); };

    return RootTaskGenerator(std::move(candidates), sample_problem, extract_features, reference,
                              distance, tolerance, seed);
}

} // namespace rsc
