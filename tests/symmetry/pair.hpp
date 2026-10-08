#pragma once

// A pair of paths to one answer, checked on inputs drawn by measure, against a
// tolerance that comes from what each path can promise rather than from the eye.
//
// The library's deepest bugs were found by two paths that disagreed, not by a
// test anybody wrote, so the pair is the unit. A pair has a reference (a path the
// library already defines), a route (another way to the same answer), a sampler
// (input i, from a low-discrepancy sequence, so the same inputs every run), a
// distance between two answers, and a bound: the sum of what the two paths are
// each allowed to be off by at that input. A disagreement is an error only when
// it is larger than both bounds together.
//
// The report says how much was checked, not that it was checked: "agreed on N
// inputs", and the measure of the set where it might not, by the rule of three
// (no failure in N draws puts the failing measure below 3/N at 95 %).
//
// Kind says what the reference is. A Function is another path computing the
// same value; a Relation is an axiom the route must satisfy (exp after log).
// They are not steps of one scale: a relation holds by definition and is checked
// on samples, a function agrees with a path that can itself be wrong.

#include <spatium/algebra/monte_carlo.hpp>

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace symmetry {

enum class Kind { Function, Relation };

struct Report {
    std::string name;
    Kind kind = Kind::Function;
    std::size_t samples = 0;
    std::size_t failures = 0;
    std::size_t first_failure = 0;     // the index of the first input that failed, 1-based; 0 if none
    double worst_ratio = 0.0;          // largest distance / bound over the inputs

    // The measure of the set of inputs on which the pair may disagree, with 95 % confidence
    // when nothing failed, and the observed share when something did.
    double unmatched_measure() const {
        if (samples == 0) return 1.0;
        return failures == 0 ? 3.0 / static_cast<double>(samples)
                             : static_cast<double>(failures) / static_cast<double>(samples);
    }
};

// Input k of the Halton sequence, in [0, 1): deterministic, and well spread in the few
// dimensions a pair needs. The library's own sequence generator is private to its engine.
inline double halton(std::size_t i, std::size_t k) {
    namespace d = spatium::monte_carlo_detail;
    return d::radical_inverse(static_cast<std::uint64_t>(i), d::kPrimes[k]);
}

// sample(i): the i-th input (i = 1..n). ref(s) and route(s): the two answers.
// distance(a, b): how far apart they are, as a number. bound(s): what the two paths
// together are allowed to differ by at s. A distance that is not finite is a failure.
template<class Sample, class Sampler, class Ref, class Route, class Distance, class Bound>
Report check(std::string name, Kind kind, std::size_t n, Sampler sample, Ref ref, Route route,
             Distance distance, Bound bound) {
    Report r;
    r.name = std::move(name);
    r.kind = kind;
    r.samples = n;
    for (std::size_t i = 1; i <= n; ++i) {
        const Sample s = sample(i);
        const double dist = distance(ref(s), route(s));
        const double tol = bound(s);
        const bool ok = std::isfinite(dist) && dist <= tol;
        if (!ok) {
            ++r.failures;
            if (r.first_failure == 0) r.first_failure = i;
        }
        const double ratio = std::isfinite(dist) ? dist / tol : std::numeric_limits<double>::infinity();
        if (ratio > r.worst_ratio) r.worst_ratio = ratio;
    }
    return r;
}

struct Entry {
    std::string name;
    std::size_t min_samples = 0;       // fewer than this and the check says too little to count
    std::function<Report()> run;
};

// The registry. A pair is a file in tests/symmetry/pairs/, and the file registers itself: a new
// pair is a new file, nothing shared is edited, so any number of them can be written at once and
// each compiles on its own. `build(n)` makes the entry for n inputs (the sweep asks for 300 000).
struct Registered {
    std::string name;
    std::size_t default_samples = 0;
    std::function<Entry(std::size_t)> build;
};

inline std::vector<Registered>& registry() {
    static std::vector<Registered> r;
    return r;
}

struct Registrar {
    Registrar(std::string name, std::size_t default_samples, std::function<Entry(std::size_t)> build) {
        registry().push_back({std::move(name), default_samples, std::move(build)});
    }
};

}  // namespace symmetry
