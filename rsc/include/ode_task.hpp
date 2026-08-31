#pragma once

// Domain 4 (Cauchy/IVP dispatch), configured over comparison_task.hpp,
// fourth instance of that shape. See ode_ops.hpp for the three test
// families (Decay/Oscillator/CircularOrbit) and why CircularOrbit -- not
// general eccentric Kepler orbits -- is what "subsumes orbital
// mechanics" means concretely here.
//
// Output is std::vector<double> rather than a fixed-size array, unlike
// precision_task.hpp/rootfind_task.hpp -- the three families have
// genuinely different state dimensions (1/2/4), and candidate vs.
// reference always compare within one family's own dimensionality per
// sample, so a dynamic size costs nothing and avoids forcing an
// artificial common fixed width.
//
// Feature = log(characteristic_frequency * dt) - mean, centered the same
// way rootfind_task.hpp's |f'(x0)| was -- applied up front, not
// discovered as a bug, now the second time this exact shape (an
// always-nonnegative product feature) has come up. Two extra one-hot
// bits (is_oscillator, is_orbit; decay is the implied [0,0] case)
// distinguish which family a given freq*dt belongs to, since the three
// families' actual error-vs-freq*dt curves aren't identical.

#include <comparison_task.hpp>
#include <ode_ops.hpp>
#include <cmath>

namespace rsc {

struct OdeProblem {
    OdeFamily family;
    double rate; // k (Decay), w (Oscillator), GM (CircularOrbit)
    double dt;
    int steps;
};

using OdeOutput = std::vector<double>;
using OdeTaskGenerator = ComparisonTaskGenerator<OdeProblem, OdeOutput>;

inline OdeTaskGenerator build_ode_task_generator(const Registry& registry, std::uint64_t seed = 0,
                                                  double tolerance = 1e-3) {
    auto call = [&registry](std::size_t op_index, const OdeProblem& p) {
        std::vector<double> in{static_cast<double>(static_cast<int>(p.family)), p.rate, p.dt,
                                static_cast<double>(p.steps)};
        std::vector<double> out(4);
        registry[op_index](in, out);
        std::size_t n = p.family == OdeFamily::Decay     ? 1
                        : p.family == OdeFamily::Oscillator ? 2
                                                             : 4;
        return OdeOutput(out.begin(), out.begin() + static_cast<long>(n));
    };

    std::size_t euler_index = registry.index_of("ode_euler");
    std::size_t rk4_index = registry.index_of("ode_rk4");

    std::vector<Candidate<OdeProblem, OdeOutput>> candidates{
        {"ode_euler", [call, euler_index](const OdeProblem& p) { return call(euler_index, p); }},
        {"ode_rk4", [call, rk4_index](const OdeProblem& p) { return call(rk4_index, p); }},
    };

    constexpr int kSteps = 20;

    auto sample_problem = [](std::mt19937_64& rng) -> OdeProblem {
        std::uniform_int_distribution<int> family_dist(0, 2);
        auto family = static_cast<OdeFamily>(family_dist(rng));

        std::uniform_real_distribution<double> rate_dist(0.5, 3.0);
        double rate = rate_dist(rng);
        double freq = ode_detail::characteristic_frequency(family, rate);

        std::bernoulli_distribution risky(0.5);
        std::uniform_real_distribution<double> risky_prod(2.0, 5.0);
        std::uniform_real_distribution<double> safe_prod(0.05, 1.0);
        double target_prod = risky(rng) ? risky_prod(rng) : safe_prod(rng);
        double dt = target_prod / freq;

        return {family, rate, dt, kSteps};
    };

    // Measured, not guessed: precompute feature_mean over the real
    // sample_problem distribution, same discipline as
    // rootfind_task.hpp's warmup pass.
    double feature_mean = 0.0;
    {
        std::mt19937_64 warmup(seed + 1);
        constexpr int kWarmupSamples = 20000;
        for (int i = 0; i < kWarmupSamples; ++i) {
            auto p = sample_problem(warmup);
            double freq = ode_detail::characteristic_frequency(p.family, p.rate);
            feature_mean += std::log(freq * p.dt);
        }
        feature_mean /= kWarmupSamples;
    }

    auto extract_features = [feature_mean](const OdeProblem& p) -> std::vector<double> {
        double freq = ode_detail::characteristic_frequency(p.family, p.rate);
        double log_prod = std::log(freq * p.dt) - feature_mean;
        double is_oscillator = p.family == OdeFamily::Oscillator ? 1.0 : 0.0;
        double is_orbit = p.family == OdeFamily::CircularOrbit ? 1.0 : 0.0;
        return {log_prod, is_oscillator, is_orbit};
    };

    auto reference = [](const OdeProblem& p) -> OdeOutput {
        double T = p.dt * p.steps;
        return ode_detail::exact(p.family, p.rate, T);
    };

    auto distance = [](const OdeOutput& a, const OdeOutput& b) {
        double m = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, std::abs(a[i] - b[i]));
        return m;
    };

    return OdeTaskGenerator(std::move(candidates), sample_problem, extract_features, reference,
                             distance, tolerance, seed);
}

} // namespace rsc
