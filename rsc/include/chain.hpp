#pragma once

// The "system"/chains step of the RSC build plan, and the DSL question
// resolved alongside it (see spatium/rsc/README.md's DSL section): a chain
// is a sequence of Tasks sharing one evolving accumulator, not i.i.d.
// samples like TaskGenerator. step[i].inputs = {accumulator_i, operand_i},
// step[i].op_index is the ground-truth op, step[i].expected_output =
// {accumulator_{i+1}}.
//
// This reuses Task/grade()/features()/reinforce_update() completely
// unchanged -- a chain step and an i.i.d. Tier-1 task are graded and
// trained identically, only how they're generated (a persistent
// accumulator vs. a fresh random task) differs. That was the point of
// keeping those pieces generic from the start.
//
// Teacher forcing: the accumulator always follows the chain's own
// ground-truth trajectory (registry[op_index] applied to the true
// operands), never whatever the model actually predicted -- one wrong
// guess at step i doesn't corrupt step i+1's input. Free-running (the
// model's own choice determines the next state) is a real next step, not
// attempted here -- it needs credit assignment across steps, which teacher
// forcing deliberately avoids for this first version.
//
// describe_chain() is the DSL's output side: "показывает шаги решения" --
// reuses describe() from describe.hpp per step, nothing new to parse or
// render.

#include <describe.hpp>
#include <task.hpp>
#include <train.hpp>
#include <cstdint>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace rsc {

struct ChainTask {
    double initial_state;
    std::vector<Task> steps;

    double final_result() const { return steps.back().expected_output[0]; }
};

class ChainTaskGenerator {
public:
    // registry must outlive the generator (same convention as
    // TaskGenerator). Dispatches only among registry.chainable_ops().
    ChainTaskGenerator(const Registry& registry, std::uint64_t seed = 0)
        : registry_(registry), chainable_(registry.chainable_ops()), rng_(seed) {
        if (chainable_.empty())
            throw std::invalid_argument("rsc::ChainTaskGenerator: registry has no chainable ops");
    }

    ChainTask sample(std::size_t length) {
        std::uniform_int_distribution<std::size_t> op_pick(0, chainable_.size() - 1);
        std::uniform_real_distribution<double> value_dist(-10.0, 10.0);

        ChainTask chain;
        chain.initial_state = value_dist(rng_);
        double state = chain.initial_state;

        for (std::size_t i = 0; i < length; ++i) {
            std::size_t op_index = chainable_[op_pick(rng_)];
            double operand = value_dist(rng_);

            std::vector<double> in{state, operand};
            std::vector<double> out(1);
            registry_[op_index](in, out);

            chain.steps.push_back(Task{op_index, in, out});
            state = out[0];
        }
        return chain;
    }

private:
    const Registry& registry_;
    std::vector<std::size_t> chainable_;
    std::mt19937_64 rng_;
};

// The DSL's output side: renders a solved chain as readable steps, e.g.
// "3 -> add(a=3, b=2) -> sum=5 -> multiply(a=5, b=4) -> product=20".
inline std::string describe_chain(const Registry& registry, const ChainTask& chain) {
    std::ostringstream os;
    os << chain.initial_state;
    for (const auto& step : chain.steps)
        os << ", then " << describe(registry[step.op_index].signature(), step.inputs, step.expected_output);
    return os.str();
}

// Teacher-forced REINFORCE over one chain: reinforce_update() on each step
// in order, same mechanics as a Tier-1 task, just fed the chain's own
// evolving accumulator instead of an i.i.d. sample. Returns the number of
// steps whose sampled action matched the true op.
//
// Single-sample updates (this function): measured, not assumed to work --
// on chains specifically (dispatching only between add/multiply, no
// dot3/lerp3 to make the overall number look better) accuracy stayed stuck
// at the ~0.50 chance level over 20k+ steps, never learning the underlying
// target==a+b vs. target==a*b distinction. See train_chains_batch() below,
// which fixed this by averaging gradients over a batch before applying --
// that's the one actually validated in tests/test_rsc_chain.cpp, not this.
inline std::size_t train_chain(Dispatcher& model, const ChainTask& chain, std::size_t max_in,
                                std::size_t max_out, double lr, std::mt19937_64& rng,
                                double& baseline, double baseline_lr = 0.01) {
    std::size_t correct = 0;
    for (const auto& step : chain.steps)
        if (reinforce_update(model, step, max_in, max_out, lr, rng, baseline, baseline_lr))
            ++correct;
    return correct;
}

// Minibatched REINFORCE over chains -- averages the gradient over
// batch_size whole chains (batch_size * chain_length individual steps)
// before one apply_gradients() call. This is the version that actually
// learns the add-vs-multiply distinction; see the note above.
inline double train_chains_batch(Dispatcher& model, ChainTaskGenerator& chains,
                                  std::size_t max_in, std::size_t max_out, double lr,
                                  std::size_t chain_length, std::size_t batch_size,
                                  std::size_t n_updates, std::uint64_t seed, double& baseline) {
    std::mt19937_64 rng(seed);
    std::size_t correct = 0, total = 0;
    for (std::size_t u = 0; u < n_updates; ++u) {
        Gradients sum;
        std::size_t count = 0;
        for (std::size_t b = 0; b < batch_size; ++b) {
            auto chain = chains.sample(chain_length);
            for (const auto& step : chain.steps) {
                bool ok;
                auto g = reinforce_gradient(model, step, max_in, max_out, rng, baseline, ok);
                accumulate_gradients(sum, g);
                ++count;
                if (ok) ++correct;
                ++total;
            }
        }
        scale_gradients(sum, 1.0 / static_cast<double>(count));
        model.apply_gradients(sum, lr);
    }
    return static_cast<double>(correct) / static_cast<double>(total);
}

// Held-out per-step accuracy across n_chains freshly sampled chains, no
// weight update.
inline double evaluate_chains(const Dispatcher& model, ChainTaskGenerator& chains,
                               std::size_t max_in, std::size_t max_out, std::size_t n_chains,
                               std::size_t chain_length) {
    std::size_t correct = 0, total = 0;
    for (std::size_t i = 0; i < n_chains; ++i) {
        auto chain = chains.sample(chain_length);
        for (const auto& step : chain.steps) {
            if (argmax_correct(model, step, max_in, max_out)) ++correct;
            ++total;
        }
    }
    return static_cast<double>(correct) / static_cast<double>(total);
}

} // namespace rsc
