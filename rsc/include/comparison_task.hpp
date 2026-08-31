#pragma once

// Generalizes the "compare N candidates against a reference, dispatch to
// the cheapest one accurate enough" pattern precision_task.hpp (domain 1)
// and geodesic_task.hpp (domain 2) independently arrived at -- same shape,
// different specifics, hand-duplicated twice before this existed. Adding a
// comparison-based domain now reduces to: supply candidates (name,
// compute), a problem sampler, a feature extractor, a reference, and a
// distance metric -- sample() itself is written once, here, not per domain.
//
// Candidates are checked cheapest-first (index 0 = cheapest); the first
// one within `tolerance` of the reference wins -- the AlphaDev-style
// reward shape already decided: correctness gate first, cost second.
//
// Expensive per-problem work (geodesic's Heat-method Cholesky
// factorization, e.g.) does NOT have to happen inside sample(): Problem-
// Sampler/Candidate::compute/Reference are plain std::function closures,
// so a domain that needs to precompute something expensive once does it
// while building those closures, before ever calling sample() -- this
// template just calls whatever it's given, cheap or not, and doesn't know
// or care whether the closures are backed by fresh computation or a
// precomputed cache (see geodesic_task.hpp for the caching case).

#include <task.hpp>
#include <cstdint>
#include <functional>
#include <random>
#include <string>
#include <vector>

namespace rsc {

template<typename Problem, typename Output>
struct Candidate {
    std::string name;
    std::function<Output(const Problem&)> compute;
};

template<typename Problem, typename Output>
class ComparisonTaskGenerator {
public:
    using ProblemSampler = std::function<Problem(std::mt19937_64&)>;
    using FeatureExtractor = std::function<std::vector<double>(const Problem&)>;
    using ReferenceFn = std::function<Output(const Problem&)>;
    using DistanceFn = std::function<double(const Output&, const Output&)>;

    ComparisonTaskGenerator(std::vector<Candidate<Problem, Output>> candidates,
                             ProblemSampler sample_problem, FeatureExtractor extract_features,
                             ReferenceFn reference, DistanceFn distance, double tolerance,
                             std::uint64_t seed = 0)
        : candidates_(std::move(candidates)), sample_problem_(std::move(sample_problem)),
          extract_features_(std::move(extract_features)), reference_(std::move(reference)),
          distance_(std::move(distance)), tolerance_(tolerance), rng_(seed) {}

    Task sample() {
        Problem problem = sample_problem_(rng_);
        Output ref = reference_(problem);

        // Each candidate is judged against the *last* (presumably most
        // capable) candidate's own error, not a fixed absolute constant --
        // this is what makes the same criterion correct for two genuinely
        // different real cases at once, not by coincidence:
        //   - precision_task.hpp: the last candidate (solve_cubic_real50)
        //     *is* the reference, so its distance to the reference is 0,
        //     and this reduces to "cheap candidate within tolerance of the
        //     reference" -- a gate.
        //   - geodesic_task.hpp: the last candidate (Heat) is itself only
        //     an approximation of the exact reference, with its own
        //     nonzero error. This reduces to "is the cheaper candidate at
        //     least as close to exact as Heat is" -- a closest-wins
        //     comparison, not a fixed-threshold gate. Using a fixed
        //     absolute tolerance here (tried first, this domain's original
        //     shape before the refactor) silently changed geodesic's
        //     actual ground truth and broke test_rsc_geodesic.cpp -- caught
        //     by the regression test, not assumed correct.
        double last_distance = distance_(candidates_.back().compute(problem), ref);

        std::size_t chosen = candidates_.size() - 1; // fallback: the last one
        for (std::size_t i = 0; i < candidates_.size(); ++i) {
            if (distance_(candidates_[i].compute(problem), ref) <= last_distance + tolerance_) {
                chosen = i;
                break;
            }
        }
        return Task{chosen, extract_features_(problem), {}};
    }

    const std::vector<Candidate<Problem, Output>>& candidates() const { return candidates_; }

private:
    std::vector<Candidate<Problem, Output>> candidates_;
    ProblemSampler sample_problem_;
    FeatureExtractor extract_features_;
    ReferenceFn reference_;
    DistanceFn distance_;
    double tolerance_;
    std::mt19937_64 rng_;
};

} // namespace rsc
