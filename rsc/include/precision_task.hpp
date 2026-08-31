#pragma once

// Domain 1 (polynomial/precision dispatch), now expressed as configuration
// over comparison_task.hpp's generic ComparisonTaskGenerator instead of a
// bespoke class -- see comparison_task.hpp's header comment for why. This
// file is now just: what a "problem" and "output" are for this domain, how
// to sample one, how to extract features, and how to compare -- the
// sample()/dispatch logic itself lives in comparison_task.hpp, shared with
// geodesic_task.hpp.
//
// Ground truth is generated, not hand-labeled: every cubic here has three
// real roots by construction (Vieta's formulas from three sampled roots),
// deliberately mixing two regimes -- widely-separated roots (easy, double
// is fine) and closely-spaced roots (delicate, casus irreducibilis's
// trig-formula cancellation can push double's answer measurably off from
// Real50's). solve_cubic_real50 doubles as both a candidate and the
// reference: it's the more expensive, presumed-adequate option, so
// comparing solve_cubic_f64 against it directly is exactly the right
// question ("is the cheap candidate close enough to the expensive one").
// solve_cubic_f64 is correct unless the two disagree beyond `tolerance` on
// any root component, in which case solve_cubic_real50 is.
//
// r0 (the root cluster's center) is sampled from [-1,1], not [-10,10] --
// measured, not assumed: a wide range let absolute root position dominate
// the raw coefficients' scale and swamp the actual separation signal, so
// a network trained on them never learned past chance. Narrowing r0
// removed that confound; raw coefficients as features, no engineered
// feature needed, reach ~0.99 held-out accuracy.

#include <comparison_task.hpp>
#include <precision_ops.hpp>
#include <algorithm>
#include <array>
#include <cmath>

namespace rsc {

using PrecisionProblem = std::array<double, 4>; // a, b, c, d
using PrecisionOutput = std::array<double, 6>;  // 3 roots, flattened (re, im) pairs

using PrecisionTaskGenerator = ComparisonTaskGenerator<PrecisionProblem, PrecisionOutput>;

inline PrecisionTaskGenerator build_precision_task_generator(const Registry& registry,
                                                               std::uint64_t seed = 0,
                                                               double tolerance = 1e-6) {
    auto call = [&registry](std::size_t op_index, const PrecisionProblem& p) {
        std::vector<double> in(p.begin(), p.end());
        std::vector<double> out(6);
        registry[op_index](in, out);
        PrecisionOutput result;
        std::copy(out.begin(), out.end(), result.begin());
        return result;
    };

    std::size_t f64_index = registry.index_of("solve_cubic_f64");
    std::size_t real50_index = registry.index_of("solve_cubic_real50");

    std::vector<Candidate<PrecisionProblem, PrecisionOutput>> candidates{
        {"solve_cubic_f64",
         [call, f64_index](const PrecisionProblem& p) { return call(f64_index, p); }},
        {"solve_cubic_real50",
         [call, real50_index](const PrecisionProblem& p) { return call(real50_index, p); }},
    };

    auto sample_problem = [](std::mt19937_64& rng) -> PrecisionProblem {
        std::uniform_real_distribution<double> r0_dist(-1.0, 1.0);
        std::bernoulli_distribution regime(0.5);
        bool delicate = regime(rng);
        std::uniform_real_distribution<double> gap_dist(delicate ? 1e-6 : 1.0,
                                                          delicate ? 1e-2 : 10.0);
        std::uniform_real_distribution<double> sign_dist(-1.0, 1.0);

        double r0 = r0_dist(rng);
        double r1 = r0 + (sign_dist(rng) < 0 ? -1.0 : 1.0) * gap_dist(rng);
        double r2 = r1 + (sign_dist(rng) < 0 ? -1.0 : 1.0) * gap_dist(rng);

        double b = -(r0 + r1 + r2);
        double c = r0 * r1 + r0 * r2 + r1 * r2;
        double d = -(r0 * r1 * r2);
        return {1.0, b, c, d};
    };

    auto extract_features = [](const PrecisionProblem& p) -> std::vector<double> {
        return {p[0], p[1], p[2], p[3]};
    };

    auto reference = [call, real50_index](const PrecisionProblem& p) {
        return call(real50_index, p);
    };

    auto distance = [](const PrecisionOutput& a, const PrecisionOutput& b) {
        double max_dev = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i) max_dev = std::max(max_dev, std::abs(a[i] - b[i]));
        return max_dev;
    };

    return PrecisionTaskGenerator(std::move(candidates), sample_problem, extract_features,
                                   reference, distance, tolerance, seed);
}

} // namespace rsc
