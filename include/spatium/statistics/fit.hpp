#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <algorithm>
#  include <concepts>
#  include <cstddef>
#  include <numeric>
#  include <random>
#  include <span>
#  include <type_traits>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::statistics {

// ── Fitting parameters to a sample ────────────────────────────────
//
// One pass over a sample in an order drawn by `rng`: for each item the caller's objective says
// how wrong the parameters are on it (`loss`) and which way to move them (`gradient`), and the
// step rule (algebra/analysis/optimize.hpp) moves them. What the parameters belong to -- a
// softmax, a net -- and what an item is are the caller's; the pass knows neither. It returns
// the mean loss before the steps, the number a training loop watches.

template <class Params> struct Fit {
    double loss = 0;
    Params gradient{};
};

template <class Params, class Rule, class Item, class Objective>
    requires requires(Rule r, Params& p, const Params& g, double rate) { r.step(p, g, rate); } &&
             std::convertible_to<std::invoke_result_t<Objective&, const Item&>, Fit<Params>>
double fit_pass(Params& params, Rule& rule, std::span<const Item> sample, double rate,
                std::mt19937_64& rng, Objective&& objective) {
    std::vector<std::size_t> order(sample.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::ranges::shuffle(order, rng);
    double total = 0;
    for (std::size_t i : order) {
        const Fit<Params> f = objective(sample[i]);
        total += f.loss;
        rule.step(params, f.gradient, rate);
    }
    return sample.empty() ? 0.0 : total / static_cast<double>(sample.size());
}

} // namespace spatium::statistics
