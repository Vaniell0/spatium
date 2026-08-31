#pragma once

// Unified base: one combined classification head + feature space spanning
// every Eigen-independent domain built so far (Tier-1, precision,
// root-finding, Cauchy/IVP) -- "base ('болванка'): Tier1 + broad Tier2
// dispatch competence" in the design doc was always meant to be one
// combined artifact, not five separately-trained dispatchers that never
// see each other's tasks. Geodesic dispatch is deliberately NOT included
// here: it depends on SPATIUM_HAS_EIGEN, and a Dispatcher's input/output
// dimensions are fixed at construction -- conditionally changing them
// based on a compile flag is a real complication, left as a follow-up,
// not solved here. Geodesic stays its own separate specialist model.
//
// Combined feature vector (22 slots): [domain one-hot (4), tier1 (10 =
// max_in 7 + max_out 3, via features()'s existing [inputs,expected_output]
// convention), precision (4), rootfind (1), ode (3)]. Explicit one-hot
// plus non-overlapping per-domain slots, not a shared/overlapping
// layout -- avoids precision's coefficient `a` and rootfind's single
// log-feature occupying the same physical slot with the network having
// no way to tell them apart besides pattern-guessing from which other
// slots are zero.
//
// Combined label space (10 classes): tier1 keeps indices 0-3, precision
// 4-5, rootfind 6-7, ode 8-9 -- plain offsets, not a rebuilt Registry
// object. train.hpp's training/eval functions only ever compare integer
// labels against task.op_index; they never call registry[op_index]
// directly, so a literal merged Registry isn't needed for this to work.
//
// Curriculum: Tier-1 gets a fixed floor quota (0.25); the three Tier-2
// domains split the rest evenly (0.25 each). The design doc's "growing
// depth" (increasing Tier-2 share over training) is not implemented --
// fixed mixing is the simplest thing that could work; revisit only if
// training shows a real need for a schedule.

#include <dispatcher.hpp>
#include <ode_task.hpp>
#include <precision_task.hpp>
#include <registry.hpp>
#include <rootfind_task.hpp>
#include <task.hpp>
#include <tier1_ops.hpp>
#include <algorithm>
#include <array>
#include <cstdint>
#include <random>
#include <stdexcept>

namespace rsc {

enum class BaseDomain : std::size_t { Tier1 = 0, Precision = 1, Rootfind = 2, Ode = 3 };

inline constexpr std::size_t kBaseOneHotSize = 4;
inline constexpr std::size_t kBaseTier1Slots = 10;     // max_in(7) + max_out(3)
inline constexpr std::size_t kBasePrecisionSlots = 4;
inline constexpr std::size_t kBaseRootfindSlots = 1;
inline constexpr std::size_t kBaseOdeSlots = 3;
inline constexpr std::size_t kBaseFeatureDim =
    kBaseOneHotSize + kBaseTier1Slots + kBasePrecisionSlots + kBaseRootfindSlots + kBaseOdeSlots;
inline constexpr std::size_t kBaseNumOps = 4 + 2 + 2 + 2;

inline constexpr std::array<std::size_t, 4> kBaseOpOffset = {0, 4, 6, 8};
inline constexpr std::array<std::size_t, 4> kBaseFeatureOffset = {
    kBaseOneHotSize, kBaseOneHotSize + kBaseTier1Slots,
    kBaseOneHotSize + kBaseTier1Slots + kBasePrecisionSlots,
    kBaseOneHotSize + kBaseTier1Slots + kBasePrecisionSlots + kBaseRootfindSlots};

class BaseTaskGenerator {
public:
    explicit BaseTaskGenerator(std::uint64_t seed = 0)
        : tier1_registry_(build_tier1_registry()), tier1_gen_(tier1_registry_, seed),
          precision_registry_(build_precision_registry()),
          precision_gen_(build_precision_task_generator(precision_registry_, seed)),
          rootfind_registry_(build_rootfind_registry()),
          rootfind_gen_(build_rootfind_task_generator(rootfind_registry_, seed)),
          ode_registry_(build_ode_registry()),
          ode_gen_(build_ode_task_generator(ode_registry_, seed)), rng_(seed) {}

    Task sample() {
        static constexpr std::array<double, 4> kWeights = {0.25, 0.25, 0.25, 0.25};
        std::discrete_distribution<std::size_t> domain_dist(kWeights.begin(), kWeights.end());
        auto domain = static_cast<BaseDomain>(domain_dist(rng_));
        std::size_t d = static_cast<std::size_t>(domain);

        std::vector<double> x(kBaseFeatureDim, 0.0);
        x[d] = 1.0;

        std::size_t local_op = 0;
        switch (domain) {
        case BaseDomain::Tier1: {
            Task t = tier1_gen_.sample();
            auto local = features(t, tier1_registry_.max_in_size(), tier1_registry_.max_out_size());
            std::copy(local.begin(), local.end(), x.begin() + static_cast<long>(kBaseFeatureOffset[0]));
            local_op = t.op_index;
            break;
        }
        case BaseDomain::Precision: {
            Task t = precision_gen_.sample();
            std::copy(t.inputs.begin(), t.inputs.end(),
                      x.begin() + static_cast<long>(kBaseFeatureOffset[1]));
            local_op = t.op_index;
            break;
        }
        case BaseDomain::Rootfind: {
            Task t = rootfind_gen_.sample();
            std::copy(t.inputs.begin(), t.inputs.end(),
                      x.begin() + static_cast<long>(kBaseFeatureOffset[2]));
            local_op = t.op_index;
            break;
        }
        case BaseDomain::Ode: {
            Task t = ode_gen_.sample();
            std::copy(t.inputs.begin(), t.inputs.end(),
                      x.begin() + static_cast<long>(kBaseFeatureOffset[3]));
            local_op = t.op_index;
            break;
        }
        }

        std::size_t global_op = kBaseOpOffset[d] + local_op;
        return Task{global_op, std::move(x), {}};
    }

private:
    Registry tier1_registry_;
    TaskGenerator tier1_gen_;
    Registry precision_registry_;
    PrecisionTaskGenerator precision_gen_;
    Registry rootfind_registry_;
    RootTaskGenerator rootfind_gen_;
    Registry ode_registry_;
    OdeTaskGenerator ode_gen_;
    std::mt19937_64 rng_;
};

// Reads back which domain produced a task, from its one-hot slots --
// the same information sample() encoded, recovered without threading a
// separate side-channel through the generic Task type. Needed for
// per-domain REINFORCE baselines (see the "single shared baseline"
// finding in rsc/README.md's base-model section): a pooled reward
// baseline miscalibrates the advantage signal for whichever domain's
// difficulty sits far from the pooled average.
inline BaseDomain domain_of(const Task& task) {
    for (std::size_t d = 0; d < kBaseOneHotSize; ++d)
        if (task.inputs[d] == 1.0) return static_cast<BaseDomain>(d);
    throw std::logic_error("base_task::domain_of: no one-hot bit set");
}

} // namespace rsc
