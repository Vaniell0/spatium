// Find a sequence of registered operations that reaches a target, under a
// restriction on which operations may be used.
//
// The skeleton of the synthesis kind from the catalogue in rsc/README.md,
// and deliberately with no model in it. That is not a simplification to
// be corrected later -- it is what the measurement says. A restriction
// like "using only addition" is a *restricted action space*, not a hint,
// so the space it produces is small, and on small spaces breadth-first
// beats a learned heuristic on every axis at once: measured 2026-09-18, a
// policy-ordered search expanded 1.6-1.7x more nodes than uninformed
// breadth-first while solving less often and returning a shortest chain
// half as often (rsc/tools/search_heuristic.cpp). A model belongs here
// only once a measurement shows search running out of budget, which is
// where the restriction is loose and the space is large.
//
// The state is the multiset of values available so far, not the last value
// computed. That matters: a chain over a single accumulator collapses,
// because different sequences reach the same number and the search folds
// them together -- 9 765 625 paths at depth 10 became 13 836 distinct
// states on the toy that measured it. A multiset keeps intermediate
// results usable and distinguishes sequences that a scalar accumulator
// would merge.
//
// The oracle is free, which is the part worth noticing. Nothing here is
// labelled by hand: a candidate is right when running it produces the
// target, and the registry's own operations do the running. The library
// is its own specification.
//
// A refusal is a result. "Unreachable within N steps under this
// restriction" answers the question that was asked, and it is the answer a
// search can give that a lookup cannot.
//
// Build: part of the rsc tools target.

#include <registry.hpp>
#include <tier1_ops.hpp>

#include <algorithm>
#include <cmath>
#include <deque>
#include <format>
#include <optional>
#include <print>
#include <set>
#include <string>
#include <vector>

namespace {

// One applied operation, kept so a found chain can be printed as the
// derivation it is rather than as a bare answer.
struct Step {
    std::size_t op;
    double a, b, result;
};

struct State {
    std::vector<double> values;   // sorted, so equivalent multisets compare equal
    std::vector<Step> history;
};

std::vector<double> sorted_with(const std::vector<double>& v, std::size_t drop_i,
                                std::size_t drop_j, double added) {
    std::vector<double> out;
    out.reserve(v.size() - 1);
    for (std::size_t k = 0; k < v.size(); ++k)
        if (k != drop_i && k != drop_j) out.push_back(v[k]);
    out.push_back(added);
    std::sort(out.begin(), out.end());
    return out;
}

// Binary operations only -- `in_size == 2 && out_size == 1`, the same
// arity rule Registry::chainable_ops() already uses to decide what can
// feed its own output back in.
std::vector<std::size_t> binary_ops(const rsc::Registry& reg,
                                    const std::vector<std::string>& allowed) {
    std::vector<std::size_t> out;
    for (std::size_t i = 0; i < reg.size(); ++i) {
        const auto& sig = reg[i].signature();
        if (sig.in_size != 2 || sig.out_size != 1) continue;
        if (!allowed.empty() &&
            std::find(allowed.begin(), allowed.end(), sig.name) == allowed.end())
            continue;
        out.push_back(i);
    }
    return out;
}

struct Result {
    bool found = false;
    std::vector<Step> chain;
    std::size_t expanded = 0;
};

Result synthesize(const rsc::Registry& reg, const std::vector<std::size_t>& ops,
                  const std::vector<double>& inputs, double target, std::size_t max_steps,
                  double tolerance = 1e-9) {
    Result res;
    std::vector<double> start = inputs;
    std::sort(start.begin(), start.end());

    std::deque<State> queue{State{start, {}}};
    std::set<std::vector<double>> seen{start};

    while (!queue.empty()) {
        const State cur = queue.front();
        queue.pop_front();
        ++res.expanded;
        if (cur.history.size() >= max_steps) continue;

        for (std::size_t i = 0; i < cur.values.size(); ++i) {
            for (std::size_t j = 0; j < cur.values.size(); ++j) {
                if (i == j) continue;
                for (std::size_t op : ops) {
                    const double in[2] = {cur.values[i], cur.values[j]};
                    double out = 0.0;
                    reg[op](std::span<const double>(in, 2), std::span<double>(&out, 1));
                    if (!std::isfinite(out)) continue;

                    State next{sorted_with(cur.values, i, j, out), cur.history};
                    next.history.push_back(Step{op, in[0], in[1], out});

                    if (std::abs(out - target) <= tolerance) {
                        res.found = true;
                        res.chain = next.history;
                        return res;
                    }
                    if (seen.insert(next.values).second) queue.push_back(std::move(next));
                }
            }
        }
    }
    return res;
}

std::string render(const rsc::Registry& reg, const std::vector<Step>& chain) {
    std::string out;
    for (const auto& s : chain)
        out += std::format("{}{}({}, {}) = {}", out.empty() ? "" : " -> ",
                           reg[s.op].signature().name, s.a, s.b, s.result);
    return out;
}

void run(const rsc::Registry& reg, const char* label, const std::vector<double>& inputs,
         double target, const std::vector<std::string>& allowed, std::size_t max_steps) {
    const auto ops = binary_ops(reg, allowed);
    std::string names;
    for (std::size_t op : ops) names += (names.empty() ? "" : ", ") + reg[op].signature().name;

    std::println("{}", label);
    std::println("  allowed: {{{}}}   max steps: {}", names, max_steps);

    const Result r = synthesize(reg, ops, inputs, target, max_steps);
    if (r.found)
        std::println("  found in {} steps, {} states expanded:\n    {}",
                     r.chain.size(), r.expanded, render(reg, r.chain));
    else
        std::println("  unreachable within {} steps under this restriction "
                     "({} states expanded)", max_steps, r.expanded);
    std::println("");
}

}  // namespace

int main() {
    const rsc::Registry reg = rsc::build_tier1_registry();

    std::println("Synthesis over the registry: find a sequence, not an answer");
    std::println("No model. The restriction is the action space, and on a small action");
    std::println("space breadth-first is the measured best option.");
    std::println("");

    // Reachable with both operations available.
    run(reg, "1. reach 21 from {1, 2, 3, 4}", {1, 2, 3, 4}, 21.0, {}, 3);

    // The same target, addition only. The answer changes, and so may the
    // existence of one -- which is the whole point of restricting.
    run(reg, "2. the same, using only addition", {1, 2, 3, 4}, 21.0, {"add"}, 3);

    // A refusal, stated as a result. 100 is not reachable from these four
    // values by three additions, and saying so is an answer.
    run(reg, "3. reach 100 from {1, 2, 3, 4}, addition only", {1, 2, 3, 4}, 100.0, {"add"}, 3);

    // Multiplication alone can reach it, which is the comparison the
    // restriction exists to make legible.
    run(reg, "4. reach 24 from {1, 2, 3, 4}, multiplication only", {1, 2, 3, 4}, 24.0,
        {"multiply"}, 3);

    std::println("A refusal is a result: \"unreachable in N steps under this restriction\"");
    std::println("is what was asked, and is the kind of answer a search can give and a");
    std::println("table of results cannot.");
    return 0;
}
