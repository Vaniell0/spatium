#pragma once

// Rigid-body integrator dispatch, configured over comparison_task.hpp --
// same shape as ode_task.hpp, five candidates instead of two. See
// integrator_ops.hpp for the three test families and why a single
// characteristic-frequency feature (sqrt(rate)) applies uniformly across
// all of them.
//
// Feature = log(freq * dt) - mean, same centering as ode_task.hpp's, plus
// two one-hot bits (is_spring, is_point_gravity; uniform-gravity is the
// implied [0,0] case).
//
// target_prod (freq*dt) is drawn log-uniform over a single wide range
// rather than ode_task.hpp's two-tier risky/safe split: with five
// dispatch classes instead of two, a smooth spread across the full
// quality range exercises every tier, where a bimodal split tuned for a
// binary choice would leave the middle candidates (verlet, rk4)
// under-sampled.

#include <comparison_task.hpp>
#include <integrator_ops.hpp>
#include <cmath>

namespace rsc {

struct IntegratorProblem {
    IntegratorFamily family;
    double rate; // g (UniformGravity), k (Spring), GM (PointGravity)
    double dt;
    int steps;
};

using IntegratorOutput = std::vector<double>;
using IntegratorTaskGenerator = ComparisonTaskGenerator<IntegratorProblem, IntegratorOutput>;

inline IntegratorTaskGenerator build_integrator_task_generator(const Registry& registry,
                                                                 std::uint64_t seed = 0,
                                                                 double tolerance = 1e-3) {
    auto call = [&registry](std::size_t op_index, const IntegratorProblem& p) {
        std::vector<double> in{static_cast<double>(static_cast<int>(p.family)), p.rate, p.dt,
                                static_cast<double>(p.steps)};
        std::vector<double> out(4);
        registry[op_index](in, out);
        return IntegratorOutput(out.begin(), out.end());
    };

    std::size_t euler_index = registry.index_of("integrator_euler");
    std::size_t semi_implicit_index = registry.index_of("integrator_semi_implicit_euler");
    std::size_t verlet_index = registry.index_of("integrator_verlet");
    std::size_t rk4_index = registry.index_of("integrator_rk4");
    std::size_t yoshida4_index = registry.index_of("integrator_yoshida4");

    std::vector<Candidate<IntegratorProblem, IntegratorOutput>> candidates{
        {"integrator_euler",
         [call, euler_index](const IntegratorProblem& p) { return call(euler_index, p); }},
        {"integrator_semi_implicit_euler",
         [call, semi_implicit_index](const IntegratorProblem& p) {
             return call(semi_implicit_index, p);
         }},
        {"integrator_verlet",
         [call, verlet_index](const IntegratorProblem& p) { return call(verlet_index, p); }},
        {"integrator_rk4",
         [call, rk4_index](const IntegratorProblem& p) { return call(rk4_index, p); }},
        {"integrator_yoshida4",
         [call, yoshida4_index](const IntegratorProblem& p) { return call(yoshida4_index, p); }},
    };

    constexpr int kSteps = 20;
    constexpr double kLogProdMin = -4.6; // ln(0.01)
    constexpr double kLogProdMax = 2.1;  // ln(8.0)

    auto sample_problem = [](std::mt19937_64& rng) -> IntegratorProblem {
        std::uniform_int_distribution<int> family_dist(0, 2);
        auto family = static_cast<IntegratorFamily>(family_dist(rng));

        std::uniform_real_distribution<double> rate_dist(0.5, 3.0);
        double rate = rate_dist(rng);
        double freq = integrator_detail::characteristic_frequency(rate);

        std::uniform_real_distribution<double> log_prod_dist(kLogProdMin, kLogProdMax);
        double target_prod = std::exp(log_prod_dist(rng));
        double dt = target_prod / freq;

        return {family, rate, dt, kSteps};
    };

    // Measured, not guessed -- same discipline as ode_task.hpp's warmup pass.
    double feature_mean = 0.0;
    {
        std::mt19937_64 warmup(seed + 1);
        constexpr int kWarmupSamples = 20000;
        for (int i = 0; i < kWarmupSamples; ++i) {
            auto p = sample_problem(warmup);
            double freq = integrator_detail::characteristic_frequency(p.rate);
            feature_mean += std::log(freq * p.dt);
        }
        feature_mean /= kWarmupSamples;
    }

    auto extract_features = [feature_mean](const IntegratorProblem& p) -> std::vector<double> {
        double freq = integrator_detail::characteristic_frequency(p.rate);
        double log_prod = std::log(freq * p.dt) - feature_mean;
        double is_spring = p.family == IntegratorFamily::Spring ? 1.0 : 0.0;
        double is_point_gravity = p.family == IntegratorFamily::PointGravity ? 1.0 : 0.0;
        return {log_prod, is_spring, is_point_gravity};
    };

    auto reference = [](const IntegratorProblem& p) -> IntegratorOutput {
        double T = p.dt * p.steps;
        return integrator_detail::exact(p.family, p.rate, T);
    };

    auto distance = [](const IntegratorOutput& a, const IntegratorOutput& b) {
        double m = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, std::abs(a[i] - b[i]));
        return m;
    };

    return IntegratorTaskGenerator(std::move(candidates), sample_problem, extract_features,
                                    reference, distance, tolerance, seed);
}

} // namespace rsc
