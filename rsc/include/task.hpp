#pragma once

// Tier-1 task generator: the direct analog of noesis's train_chain() in
// experiments/A0_state_probe/micro_wkv.py, but against a real Registry
// instead of a two-op toy. A Task is a single-op dispatch problem — random
// operands for a randomly chosen registered op, with the correct answer
// computed by the registry itself (exact by construction, no model
// involved). Chains (composing several ops) are a later step once
// single-op dispatch works — starting there would mean debugging
// composition and dispatch at the same time.
//
// grade() is the hard pass/fail gate from the RSC design doc: wrong op
// choice or an output outside tolerance both count as failure, nothing in
// between -- this is what a trained dispatcher's reward will be built on.

#include <registry.hpp>
#include <cmath>
#include <cstdint>
#include <random>
#include <span>
#include <vector>

namespace rsc {

struct Task {
    std::size_t op_index;                // ground truth: which op solves this
    std::vector<double> inputs;           // size == registry[op_index].signature().in_size
    std::vector<double> expected_output;  // registry[op_index] applied to inputs
};

class TaskGenerator {
public:
    // registry must outlive the generator -- same convention as passing a
    // const Space& into Spatium's own concept-checked algorithms.
    TaskGenerator(const Registry& registry, std::uint64_t seed = 0)
        : registry_(registry), rng_(seed) {}

    Task sample() {
        std::uniform_int_distribution<std::size_t> op_dist(0, registry_.size() - 1);
        std::size_t op_index = op_dist(rng_);
        const auto& sig = registry_[op_index].signature();

        std::uniform_real_distribution<double> value_dist(-10.0, 10.0);
        std::vector<double> inputs(sig.in_size);
        for (auto& v : inputs) v = value_dist(rng_);

        std::vector<double> output(sig.out_size);
        registry_[op_index](inputs, output);

        return Task{op_index, std::move(inputs), std::move(output)};
    }

private:
    const Registry& registry_;
    std::mt19937_64 rng_;
};

struct Grade {
    bool op_correct;
    bool output_correct;
    bool passed() const { return op_correct && output_correct; }
};

inline Grade grade(const Task& task, std::size_t proposed_op_index,
                    std::span<const double> proposed_output, double tol = 1e-9) {
    Grade g{};
    g.op_correct = (proposed_op_index == task.op_index);
    g.output_correct = proposed_output.size() == task.expected_output.size();
    if (g.output_correct) {
        for (std::size_t i = 0; i < proposed_output.size(); ++i) {
            if (std::abs(proposed_output[i] - task.expected_output[i]) > tol) {
                g.output_correct = false;
                break;
            }
        }
    }
    return g;
}

} // namespace rsc
