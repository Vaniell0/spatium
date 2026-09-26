// RSC labelling a continuous-collision choice: the one-sided rule has to
// reach the labels, or a model is trained to prefer answers that are late.
#include <catch2/catch_test_macros.hpp>

#include <comparison_task.hpp>

#include <cmath>
#include <random>

namespace {

// Three candidates for a time of impact, cheapest first: one that answers a
// hair late, one that answers far too early to be of use, one that answers
// a hair early. Only the last is admissible.
struct Toi { double toi; };

rsc::ComparisonTaskGenerator<double, Toi> make_generator(bool one_sided) {
    std::vector<rsc::Candidate<double, Toi>> c{
        {"late by 1e-6", [](const double& t) { return Toi{t + 1e-6}; }},
        {"early by 0.05", [](const double& t) { return Toi{t - 0.05}; }},
        {"early by 1e-5", [](const double& t) { return Toi{t - 1e-5}; }},
    };
    rsc::ComparisonTaskGenerator<double, Toi> g(
        std::move(c), [](std::mt19937_64& rng) { return std::uniform_real_distribution<double>(0.2, 0.8)(rng); },
        [](const double& t) { return std::vector<double>{t}; }, [](const double& t) { return Toi{t}; },
        [](const Toi& a, const Toi& b) { return std::abs(a.toi - b.toi); }, 1e-4, 3);
    if (one_sided)
        g.with_gate([](const Toi& cand, const Toi& ref) {
            // Never late; early, but by no more than it takes to stay useful.
            return cand.toi <= ref.toi && cand.toi >= ref.toi - 1e-3;
        });
    return g;
}

}  // namespace

// The symmetric gate took the cheapest candidate within 1e-4 of the
// reference either way, and the cheapest is late: every label it made
// taught a model to prefer the one answer continuous collision must never
// give.
TEST_CASE("A continuous-collision label is never a late answer", "[rsc][ccd]") {
    auto g = make_generator(/*one_sided=*/true);
    for (int i = 0; i < 50; ++i) {
        const auto task = g.sample();
        INFO("sample " << i << " labelled " << g.candidates()[task.op_index].name);
        CHECK(task.op_index == 2);
    }
}

// And without a gate the template behaves exactly as the seven domains
// built on it rely on.
TEST_CASE("Without a gate the comparison is the symmetric one it always was", "[rsc][ccd]") {
    auto g = make_generator(/*one_sided=*/false);
    CHECK(g.sample().op_index == 0);
}
