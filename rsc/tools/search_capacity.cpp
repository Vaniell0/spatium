// How much model does a *search* need, as opposed to a classification?
//
// RSC's dispatcher today answers "which operation solves this problem",
// given the problem and its answer. A hidden layer of 64 units over a few
// dozen operations is ample for that, and the claim this experiment exists
// to test is that it is ample because the task is closer to a learned
// decision tree than to a model: each input determines its label, and
// nothing has to be held in mind across steps.
//
// Search is not that. Here the model is given a start and a target and
// must choose a sequence -- with no per-step supervision, and a reward
// only at the end, so credit for reaching the target has to be assigned
// backwards across every choice that led there. That is exactly what
// chain.hpp's teacher forcing was built to avoid, and it is the setting in
// which capacity might start to matter.
//
// Build: part of the rsc tools target.

#include "toy_search.hpp"

#include <algorithm>
#include <print>
#include <vector>

using toy::kNumActions;

namespace {
// Repeats per cell, because the first version of this experiment had one
// and could not be read. Its untrained column -- same task, same eval set,
// only the initialisation seed differing -- ranged from 2.1% to 27.6%, a
// thirteenfold spread, which is larger than any difference between hidden
// sizes it went on to report. A single number per cell was measuring the
// seed lottery and presenting it as capacity.
constexpr std::size_t kSeeds = 5;

double median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    const std::size_t n = v.size();
    return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

}  // namespace

int main() {
    constexpr std::size_t kInputDim = 6;
    constexpr std::size_t kEpisodes = 60000;
    constexpr std::size_t kBatch = 32;
    constexpr double kLearningRate = 0.05;
    constexpr std::size_t kEvalTrials = 2000;

    const std::vector<std::size_t> hidden_sizes{8, 32, 128};
    const std::vector<std::size_t> walks{2, 3, 4};

    std::println("RSC search capacity -- five integer actions, terminal reward only");
    std::println("{} episodes per cell, batch {}, lr {}, {} eval problems, {} seeds per cell",
                 kEpisodes, kBatch, kLearningRate, kEvalTrials, kSeeds);
    std::println("Paired: each seed is evaluated before and after training on the same problems.");
    std::println("");

    for (std::size_t walk : walks) {
        std::println("walk length {} (target reached by {} random actions; BFS optimum is often shorter)",
                     walk, walk);
        std::println("  {:>7} | {:>16} | {:>16} | {:>14} | {:>10}",
                     "hidden", "untrained median", "trained median", "trained range", "gain");

        for (std::size_t h : hidden_sizes) {
            std::vector<double> before_all, after_all;
            for (std::size_t s = 0; s < kSeeds; ++s) {
                rsc::Dispatcher model(kInputDim, h, kNumActions, /*seed=*/1000 * (s + 1) + h);
                before_all.push_back(toy::evaluate(model, walk, kEvalTrials, /*seed=*/99).success);
                toy::train(model, walk, kEpisodes, kLearningRate, kBatch, /*seed=*/7 + 31 * s + h);
                after_all.push_back(toy::evaluate(model, walk, kEvalTrials, /*seed=*/99).success);
            }
            const double b = median(before_all), a = median(after_all);
            const double lo = *std::min_element(after_all.begin(), after_all.end());
            const double hi = *std::max_element(after_all.begin(), after_all.end());

            std::println("  {:>7} | {:>15.1f}% | {:>15.1f}% | {:>6.1f}-{:>5.1f}% | {:>+9.1f}%",
                         h, b * 100.0, a * 100.0, lo * 100.0, hi * 100.0, (a - b) * 100.0);
        }
        std::println("");
    }

    std::println("Read the range before the median. If the ranges for two hidden sizes");
    std::println("overlap, this experiment has not distinguished them, whatever the");
    std::println("medians do.");
    return 0;
}
