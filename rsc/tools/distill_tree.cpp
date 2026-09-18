// Does a decision tree keep what the network learned?
//
// `dispatch_cost` established that a trained dispatcher cannot reach
// branch speed by being made smaller -- hidden 8 is still 55x an `if` --
// while a depth-4 tree costs 1.6x. So the integration question is not
// "how fast can the network be" but "how much does the tree lose", and
// that is one number nobody has taken.
//
// The method is the standard one: the network is the teacher, the tree is
// fitted to imitate *its* choices rather than the ground truth, and then
// both are scored against the ground truth on data neither saw. Fidelity
// (tree agrees with network) is reported separately from accuracy (tree is
// right), because they fail differently -- a tree can imitate a network
// faithfully into the same mistakes, and a tree can disagree with the
// network while being more right than it.
//
// Per domain, never pooled. The base curriculum's own history is the
// reason: an aggregate can look healthy while one slice sits at chance,
// which is why tests/test_rsc_base.cpp asserts a floor per domain rather
// than one number.
//
// Build: part of the rsc tools target.

#include <base_task.hpp>
#include <dispatcher.hpp>
#include <train.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <memory>
#include <numeric>
#include <print>
#include <random>
#include <vector>

namespace {

struct Sample {
    std::vector<double> x;
    std::size_t teacher;   // the network's choice -- what the tree is fitted to
    std::size_t truth;     // the registry's own answer
    std::size_t domain;
};

std::size_t argmax_of(const rsc::Dispatcher& model, const std::vector<double>& x) {
    auto c = model.forward_cached(x);
    return static_cast<std::size_t>(
        std::max_element(c.logits.begin(), c.logits.end()) - c.logits.begin());
}

// Two residuals the tree cannot build for itself. Tier-1 dispatch asks
// whether the observed output is a+b or a*b, which is a relation between
// three features, and an axis-aligned split can only approximate a
// diagonal boundary by stacking steps -- which is the suspected reason
// domain 0 is the one that resists depth. Handing the tree the residuals
// makes that boundary axis-aligned; if the suspicion is right, domain 0
// stops needing depth at all, and if it is wrong, nothing moves.
//
// The network keeps its own inputs untouched. Distillation does not
// require teacher and student to see the same features, and pretending it
// does would leave a cheap question unasked.
//
// Layout: 4 one-hot slots, then the Tier-1 block of max_in(7) + max_out(3),
// so operands are at 4 and 5 and the observed output at 11.
std::vector<double> augment(const std::vector<double>& x) {
    auto out = x;
    const double a = x[4], b = x[5], y = x[11];
    out.push_back(y - (a + b));
    out.push_back(y - (a * b));
    return out;
}

std::vector<Sample> collect(const rsc::Dispatcher& model, std::uint64_t seed, std::size_t n,
                            bool augmented) {
    rsc::BaseTaskGenerator gen(seed);
    std::vector<Sample> out;
    out.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        auto task = gen.sample();
        auto x = rsc::features(task, rsc::kBaseFeatureDim, 0);
        const std::size_t teacher = argmax_of(model, x);   // network sees the raw features
        out.push_back(Sample{augmented ? augment(x) : x, teacher, task.op_index,
                             static_cast<std::size_t>(rsc::domain_of(task))});
    }
    return out;
}

// ── A depth-limited CART over the teacher's labels ───────────────

struct Node {
    bool leaf = true;
    std::size_t label = 0;
    std::size_t feature = 0;
    double threshold = 0.0;
    std::unique_ptr<Node> lo, hi;
};

std::size_t majority(const std::vector<const Sample*>& rows) {
    std::vector<std::size_t> count(rsc::kBaseNumOps, 0);
    for (const auto* r : rows) ++count[r->teacher];
    return static_cast<std::size_t>(std::max_element(count.begin(), count.end()) - count.begin());
}

double gini(const std::vector<std::size_t>& count, std::size_t total) {
    if (!total) return 0.0;
    double g = 1.0;
    for (std::size_t c : count) {
        const double p = static_cast<double>(c) / static_cast<double>(total);
        g -= p * p;
    }
    return g;
}

double gini_of_rows(const std::vector<const Sample*>& rows) {
    std::vector<std::size_t> count(rsc::kBaseNumOps, 0);
    for (const auto* r : rows) ++count[r->teacher];
    return gini(count, rows.size());
}

std::unique_ptr<Node> fit(const std::vector<const Sample*>& rows, std::size_t depth,
                          std::size_t max_depth, std::size_t min_rows) {
    auto node = std::make_unique<Node>();
    node->label = rows.empty() ? 0 : majority(rows);
    if (depth >= max_depth || rows.size() < min_rows) return node;

    // Pure already?
    bool pure = true;
    for (const auto* r : rows)
        if (r->teacher != rows.front()->teacher) { pure = false; break; }
    if (pure) return node;

    const std::size_t dim = rows.front()->x.size();
    double best_score = gini_of_rows(rows);
    std::size_t best_feature = 0;
    double best_threshold = 0.0;
    bool found = false;

    for (std::size_t f = 0; f < dim; ++f) {
        // Candidate thresholds at deciles rather than at every value: the
        // split quality is flat between adjacent samples, and this keeps
        // the fit from being quadratic in the sample count.
        std::vector<double> vals;
        vals.reserve(rows.size());
        for (const auto* r : rows) vals.push_back(r->x[f]);
        std::sort(vals.begin(), vals.end());
        for (int q = 1; q < 10; ++q) {
            const double t = vals[vals.size() * static_cast<std::size_t>(q) / 10];
            std::vector<std::size_t> lo_count(rsc::kBaseNumOps, 0), hi_count(rsc::kBaseNumOps, 0);
            std::size_t lo_n = 0, hi_n = 0;
            for (const auto* r : rows) {
                if (r->x[f] <= t) { ++lo_count[r->teacher]; ++lo_n; }
                else              { ++hi_count[r->teacher]; ++hi_n; }
            }
            if (!lo_n || !hi_n) continue;
            const double score =
                (static_cast<double>(lo_n) * gini(lo_count, lo_n) +
                 static_cast<double>(hi_n) * gini(hi_count, hi_n)) / static_cast<double>(rows.size());
            if (score < best_score - 1e-12) {
                best_score = score;
                best_feature = f;
                best_threshold = t;
                found = true;
            }
        }
    }
    if (!found) return node;

    std::vector<const Sample*> lo, hi;
    for (const auto* r : rows) (r->x[best_feature] <= best_threshold ? lo : hi).push_back(r);

    node->leaf = false;
    node->feature = best_feature;
    node->threshold = best_threshold;
    node->lo = fit(lo, depth + 1, max_depth, min_rows);
    node->hi = fit(hi, depth + 1, max_depth, min_rows);
    return node;
}

std::size_t predict(const Node& n, const std::vector<double>& x) {
    if (n.leaf) return n.label;
    return predict(x[n.feature] <= n.threshold ? *n.lo : *n.hi, x);
}

std::size_t count_leaves(const Node& n) {
    return n.leaf ? 1 : count_leaves(*n.lo) + count_leaves(*n.hi);
}

}  // namespace

int main() {
    constexpr std::size_t kHidden = 64;
    constexpr std::size_t kUpdates = 8000;
    constexpr std::size_t kBatch = 32;
    constexpr std::size_t kFitRows = 40000;
    constexpr std::size_t kEvalRows = 20000;

    rsc::Dispatcher model(rsc::kBaseFeatureDim, kHidden, rsc::kBaseNumOps, /*seed=*/0);

    // Same recipe as tests/test_rsc_base.cpp, including the per-domain
    // baselines -- a pooled baseline sits near the weighted average reward
    // and miscalibrates whichever domain's difficulty is far from it.
    rsc::BaseTaskGenerator train_gen(1);
    std::mt19937_64 rng(42);
    std::array<double, 4> baselines{0.5, 0.5, 0.5, 0.5};
    for (std::size_t u = 0; u < kUpdates; ++u) {
        rsc::Gradients sum;
        for (std::size_t b = 0; b < kBatch; ++b) {
            auto task = train_gen.sample();
            auto domain = rsc::domain_of(task);
            bool ok = false;
            auto g = rsc::reinforce_gradient(model, task, rsc::kBaseFeatureDim, 0, rng,
                                             baselines[static_cast<std::size_t>(domain)], ok,
                                             0.01, /*entropy_beta=*/0.20);
            rsc::accumulate_gradients(sum, g);
        }
        rsc::scale_gradients(sum, 1.0 / static_cast<double>(kBatch));
        model.apply_gradients(sum, 0.05);
    }

    // The network's own per-domain accuracy on the eval set -- the number
    // the tree has to be compared against, not an assumed one. Taken from
    // the raw-feature pass, since the network's inputs never change.
    const auto ref_rows = collect(model, /*seed=*/9876, kEvalRows, /*augmented=*/false);
    std::array<std::size_t, 4> net_ok{}, total{};
    for (const auto& s : ref_rows) {
        ++total[s.domain];
        if (s.teacher == s.truth) ++net_ok[s.domain];
    }

    std::println("Distilling a hidden-{} dispatcher into a decision tree", kHidden);
    std::println("{} rows to fit, {} to evaluate, {} ops, {} features",
                 kFitRows, kEvalRows, rsc::kBaseNumOps, rsc::kBaseFeatureDim);
    std::println("");
    std::print("  network accuracy per domain:");
    for (std::size_t d = 0; d < 4; ++d)
        std::print(" {}={:.3f}", d, static_cast<double>(net_ok[d]) / static_cast<double>(total[d]));
    std::println("");
    std::println("");
    for (bool augmented : {false, true}) {
    const auto fit_rows = collect(model, /*seed=*/1, kFitRows, augmented);
    const auto eval_rows = collect(model, /*seed=*/9876, kEvalRows, augmented);
    std::vector<const Sample*> ptrs;
    ptrs.reserve(fit_rows.size());
    for (const auto& s : fit_rows) ptrs.push_back(&s);

    std::println("  features: {}", augmented ? "raw + two Tier-1 residuals" : "raw");
    std::println("  {:>5} | {:>6} | {:>8} | {:>29} | {:>8}",
                 "depth", "leaves", "fidelity", "tree accuracy per domain", "worst gap");

    for (std::size_t depth : {2u, 4u, 6u, 8u, 10u, 12u}) {
        auto tree = fit(ptrs, 0, depth, /*min_rows=*/20);

        std::size_t agree = 0;
        std::array<std::size_t, 4> tree_ok{};
        for (const auto& s : eval_rows) {
            const std::size_t p = predict(*tree, s.x);
            if (p == s.teacher) ++agree;
            if (p == s.truth) ++tree_ok[s.domain];
        }

        double worst = 0.0;
        std::string acc;
        for (std::size_t d = 0; d < 4; ++d) {
            const double n = static_cast<double>(net_ok[d]) / static_cast<double>(total[d]);
            const double t = static_cast<double>(tree_ok[d]) / static_cast<double>(total[d]);
            worst = std::max(worst, n - t);
            acc += std::format(" {:.3f}", t);
        }

        std::println("  {:>5} | {:>6} | {:>7.3f} | {:>29} | {:>+7.3f}",
                     depth, count_leaves(*tree),
                     static_cast<double>(agree) / static_cast<double>(eval_rows.size()),
                     acc, -worst);
    }
    std::println("");
    }

    std::println("");
    std::println("Fidelity is agreement with the network; accuracy is agreement with the");
    std::println("registry. They fail differently, so neither stands in for the other.");
    std::println("The gap column is the worst per-domain accuracy lost, which is the");
    std::println("number that decides whether the tree can ship.");
    return 0;
}
