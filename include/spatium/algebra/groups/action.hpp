#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/concepts.hpp>
#  include <spatium/core/verify.hpp>
#  include <algorithm>
#  include <concepts>
#  include <cstddef>
#  include <ranges>
#  include <span>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::inline algebra {

// ── Group actions ─────────────────────────────────────────────────
//
// A group G acts on a set X when each element g moves the points of X,
// act(g, x), so that the identity moves nothing and acting by compose(g, h)
// is acting by h and then by g (compose applies its second argument first,
// as everywhere in algebra/). The same group can act on several sets: the
// signed permutations act on integer vectors and on the cells of a cube.
//
// The orbit of x is every point the group carries it to, its stabiliser the
// elements that leave it in place. For a finite group the two are tied by
// |G| = |orbit| * |stabiliser|, which `verify_action` checks along with the
// two axioms: a cheap pair, since the orbit and the stabiliser are found by
// separate loops.

// The two axioms are checked by verify_action, not at compile time.
template <class A, class G, class X>
concept GroupAction = Group<G>
    && requires(const A& a, const typename G::ElementType& g, const X& x) {
        { a.act(g, x) } -> std::convertible_to<X>;
    };

// A group that lists its elements.
template <class G>
concept FiniteGroup = Group<G>
    && requires(const G& g) {
        { g.elements() } -> std::ranges::forward_range;
    };

// Every point the group carries x to, sorted, each once.
template <FiniteGroup G, std::totally_ordered X, GroupAction<G, X> A>
std::vector<X> orbit(const G& group, const A& action, const X& x) {
    std::vector<X> out;
    out.reserve(static_cast<std::size_t>(std::ranges::distance(group.elements())));
    for (const auto& g : group.elements())
        out.push_back(action.act(g, x));
    std::ranges::sort(out);
    out.erase(std::ranges::unique(out).begin(), out.end());
    return out;
}

// The elements that leave x where it is, in the order elements() lists them.
template <FiniteGroup G, std::equality_comparable X, GroupAction<G, X> A>
std::vector<typename G::ElementType> stabilizer(const G& group, const A& action, const X& x) {
    std::vector<typename G::ElementType> out;
    for (const auto& g : group.elements())
        if (action.act(g, x) == x)
            out.push_back(g);
    return out;
}

// The orbits of a set of points: each point's orbit once, in the order of the
// smallest point of each. Points of `points` the group carries outside it are
// still listed in the orbits, so pass a set the group maps onto itself.
template <FiniteGroup G, std::totally_ordered X, GroupAction<G, X> A>
std::vector<std::vector<X>> orbits(const G& group, const A& action, std::span<const X> points) {
    std::vector<std::vector<X>> out;
    for (const X& x : points) {
        auto o = orbit(group, action, x);
        if (std::ranges::none_of(out, [&](const auto& seen) { return seen.front() == o.front(); }))
            out.push_back(std::move(o));
    }
    std::ranges::sort(out, {}, [](const auto& o) { return o.front(); });
    return out;
}

// The identity moves nothing, acting by compose(g, h) is acting by h then g,
// and |G| = |orbit| * |stabiliser| at every sample.
template <FiniteGroup G, std::totally_ordered X, GroupAction<G, X> A>
VerifyResult verify_action(const G& group, const A& action, std::span<const X> samples) {
    const auto order = static_cast<std::size_t>(std::ranges::distance(group.elements()));
    for (const X& x : samples) {
        if (!(action.act(group.identity(), x) == x))
            return VerifyResult::fail("act(identity, x) != x");
        for (const auto& g : group.elements())
            for (const auto& h : group.elements())
                if (!(action.act(group.compose(g, h), x) == action.act(g, action.act(h, x))))
                    return VerifyResult::fail("act(compose(g, h), x) != act(g, act(h, x))");
        if (orbit(group, action, x).size() * stabilizer(group, action, x).size() != order)
            return VerifyResult::fail("|orbit| * |stabilizer| != |G|");
    }
    return VerifyResult::ok();
}

} // namespace spatium::algebra
