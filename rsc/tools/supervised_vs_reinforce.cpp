// supervised_vs_reinforce — does the base dispatcher need REINFORCE at all?
//
// Every dispatch task carries its label: `task.op_index` is the cheapest
// correct op, found by running the candidates against a reference. The
// base was nevertheless trained by REINFORCE -- sample an action, reward 1
// if it equals the label -- which estimates, noisily, the gradient that
// cross-entropy on the same label gives exactly. The batch of 32, the
// per-domain baselines and the entropy bonus in train_base were each
// measured fixes for that noise. This runs the two side by side on the
// same generator, the same features, the same model and step size, and
// reports held-out argmax accuracy per domain as training proceeds.
//
// Configurations:
//   reinforce  hidden 64, entropy_beta 0.2, per-domain baselines (= train_base)
//   supervised hidden 64, cross-entropy on the label, nothing else
//   supervised hidden 16, the same, to see whether capacity is doing the work
//
// Usage: supervised_vs_reinforce [seeds=4] [updates=8000]

#include <base_task.hpp>
#include <train.hpp>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <print>
#include <string>
#include <thread>
#include <vector>

namespace {

enum class Method { Reinforce, Supervised };

struct Config {
    const char* name;
    Method method;
    std::size_t hidden;
};

constexpr std::array<Config, 3> kConfigs{{
    {"reinforce h64", Method::Reinforce, 64},
    {"supervised h64", Method::Supervised, 64},
    {"supervised h16", Method::Supervised, 16},
}};

constexpr std::array<int, 7> kCheckpoints{250, 500, 1000, 2000, 4000, 8000, 16000};
constexpr std::size_t kBatch = 32;
constexpr double kLr = 0.1;

// accuracy[checkpoint][domain]
using Curve = std::vector<std::array<double, 4>>;

std::array<double, 4> per_domain_accuracy(const rsc::Dispatcher& model, const std::vector<rsc::Task>& eval) {
    std::array<double, 4> hit{}, count{};
    for (const auto& t : eval) {
        const auto d = static_cast<std::size_t>(rsc::domain_of(t));
        count[d] += 1;
        hit[d] += rsc::argmax_correct(model, t, rsc::kBaseFeatureDim, 0) ? 1.0 : 0.0;
    }
    for (std::size_t d = 0; d < 4; ++d) hit[d] = count[d] > 0 ? hit[d] / count[d] : 0.0;
    return hit;
}

Curve run(const Config& c, std::uint64_t seed, int updates, const std::vector<rsc::Task>& eval) {
    rsc::Dispatcher model(rsc::kBaseFeatureDim, c.hidden, rsc::kBaseNumOps, seed);
    rsc::BaseTaskGenerator gen(1000 + seed);
    std::mt19937_64 rng(42 + seed);
    std::array<double, 4> baselines{0.5, 0.5, 0.5, 0.5};
    Curve curve;
    std::size_t next = 0;
    for (int u = 1; u <= updates; ++u) {
        rsc::Gradients sum;
        for (std::size_t b = 0; b < kBatch; ++b) {
            auto task = gen.sample();
            if (c.method == Method::Reinforce) {
                bool ok;
                auto& base = baselines[static_cast<std::size_t>(rsc::domain_of(task))];
                rsc::accumulate_gradients(
                    sum, rsc::reinforce_gradient(model, task, rsc::kBaseFeatureDim, 0, rng, base, ok, 0.01, 0.2));
            } else {
                const auto x = rsc::features(task, rsc::kBaseFeatureDim, 0);
                auto cache = model.forward_cached(x);
                rsc::accumulate_gradients(
                    sum, model.backward(cache, rsc::cross_entropy_dlogits(cache.logits, task.op_index)));
            }
        }
        rsc::scale_gradients(sum, 1.0 / static_cast<double>(kBatch));
        model.apply_gradients(sum, kLr);
        if (next < kCheckpoints.size() && u == kCheckpoints[next]) {
            curve.push_back(per_domain_accuracy(model, eval));
            ++next;
        }
    }
    return curve;
}

}  // namespace

int main(int argc, char** argv) {
    const int seeds = argc > 1 ? std::atoi(argv[1]) : 4;
    const int updates = argc > 2 ? std::atoi(argv[2]) : 8000;

    // One held-out set for every run, from a generator seed no run trains on.
    rsc::BaseTaskGenerator eval_gen(777);
    std::vector<rsc::Task> eval(8000);
    for (auto& t : eval) t = eval_gen.sample();

    std::vector<std::vector<Curve>> curves(kConfigs.size(), std::vector<Curve>(static_cast<std::size_t>(seeds)));
    std::vector<std::thread> threads;
    for (std::size_t c = 0; c < kConfigs.size(); ++c)
        for (int s = 0; s < seeds; ++s)
            threads.emplace_back([&, c, s] {
                curves[c][static_cast<std::size_t>(s)] =
                    run(kConfigs[c], static_cast<std::uint64_t>(s), updates, eval);
            });
    for (auto& t : threads) t.join();

    const char* domains[4] = {"tier1", "precision", "rootfind", "ode"};
    std::println("held-out argmax accuracy, mean over {} seeds (min..max), batch {}, lr {}", seeds, kBatch, kLr);
    for (std::size_t c = 0; c < kConfigs.size(); ++c) {
        std::println("\n{}", kConfigs[c].name);
        std::print("{:>8}", "updates");
        for (auto* d : domains) std::print("{:>22}", d);
        std::println("");
        for (std::size_t k = 0; k < curves[c][0].size(); ++k) {
            std::print("{:>8}", kCheckpoints[k]);
            for (std::size_t d = 0; d < 4; ++d) {
                double sum = 0, lo = 1, hi = 0;
                for (int s = 0; s < seeds; ++s) {
                    const double a = curves[c][static_cast<std::size_t>(s)][k][d];
                    sum += a; lo = std::min(lo, a); hi = std::max(hi, a);
                }
                std::print("{:>10.3f} ({:.2f}..{:.2f})", sum / seeds, lo, hi);
            }
            std::println("");
        }
    }
    return 0;
}
