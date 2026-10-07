#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/access.hpp>
#  include <spatium/core/concepts.hpp>
#  include <cmath>
#  include <initializer_list>
#  include <span>
#  include <string>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium {

// Runtime axiom verification for mathematical spaces.
// Generate random sample points, then call verify() to check
// that the space implementation satisfies the required axioms.
//
// These don't prove correctness, but catch most implementation bugs.

struct VerifyResult {
    bool passed = true;
    std::string failure;

    // Alias for discoverability
    const std::string& message() const { return failure; }

    explicit operator bool() const { return passed; }

    static VerifyResult ok() { return {}; }
    static VerifyResult fail(std::string msg) { return {false, std::move(msg)}; }
};

// ── Metric axioms ──────────────────────────────────────────────
// 1. d(x, y) >= 0                     (non-negativity)
// 2. d(x, y) == 0  iff  x == y        (identity of indiscernibles)
// 3. d(x, y) == d(y, x)               (symmetry)
// 4. d(x, z) <= d(x, y) + d(y, z)     (triangle inequality)

// Through core/access.hpp, so a space whose distance is a free function or
// derived from its log and metric is checked like one with a member.
template<class S>
    requires spaces::Distanced<S>
VerifyResult verify_metric(const S& space,
                           std::span<const spaces::point_t<S>> samples,
                           spaces::scalar_t<S> tolerance = spaces::scalar_t<S>{1e-8}) {
    using T = spaces::scalar_t<S>;
    using std::abs;
    const auto dist = [&](const auto& a, const auto& b) { return T(spaces::distance(space, a, b)); };

    for (std::size_t i = 0; i < samples.size(); ++i) {
        // Non-negativity
        auto d_ii = dist(samples[i], samples[i]);
        if (d_ii < T{0} || abs(d_ii) > tolerance)
            return VerifyResult::fail("d(x,x) != 0");

        for (std::size_t j = i + 1; j < samples.size(); ++j) {
            auto d_ij = dist(samples[i], samples[j]);
            auto d_ji = dist(samples[j], samples[i]);

            // Non-negativity
            if (d_ij < -tolerance)
                return VerifyResult::fail("d(x,y) < 0");

            // Symmetry
            if (abs(d_ij - d_ji) > tolerance)
                return VerifyResult::fail("d(x,y) != d(y,x)");

            // Triangle inequality
            for (std::size_t k = 0; k < samples.size(); ++k) {
                auto d_ik = dist(samples[i], samples[k]);
                auto d_kj = dist(samples[k], samples[j]);
                if (d_ij > d_ik + d_kj + tolerance)
                    return VerifyResult::fail("triangle inequality violated");
            }
        }
    }
    return VerifyResult::ok();
}

// ── Inner product axioms ───────────────────────────────────────
// 1. <u, v> == <v, u>                  (symmetry)
// 2. <u, u> >= 0                       (positive-definiteness)
// 3. <u, u> == 0  iff  u == 0         (definiteness)
// 4. <au, v> == a * <u, v>            (linearity in first arg)
// 5. <u+w, v> == <u,v> + <w,v>        (additivity)

template<InnerProductSpace S>
VerifyResult verify_inner_product(const S& space,
                                  std::span<const typename S::VectorType> samples,
                                  typename S::ScalarType tolerance = typename S::ScalarType{1e-8}) {
    using T = typename S::ScalarType;
    using std::abs;

    for (std::size_t i = 0; i < samples.size(); ++i) {
        // Positive-definiteness
        auto uu = space.inner(samples[i], samples[i]);
        if (uu < -tolerance)
            return VerifyResult::fail("<u,u> < 0");

        for (std::size_t j = 0; j < samples.size(); ++j) {
            // Symmetry
            auto uv = space.inner(samples[i], samples[j]);
            auto vu = space.inner(samples[j], samples[i]);
            if (abs(uv - vu) > tolerance)
                return VerifyResult::fail("<u,v> != <v,u>");

            // Linearity: <2u, v> == 2 * <u, v>
            auto scaled = samples[i] * T{2};
            auto two_uv = space.inner(scaled, samples[j]);
            if (abs(two_uv - T{2} * uv) > tolerance)
                return VerifyResult::fail("linearity violated");
        }
    }
    return VerifyResult::ok();
}

// ── Manifold axioms (exp/log roundtrip) ────────────────────────
// exp(p, log(p, q), 1) ≈ q

// Through core/access.hpp, like verify_metric.
template<class S>
    requires spaces::Exponential<S> && spaces::Logarithmic<S>
VerifyResult verify_exp_log(const S& space,
                            std::span<const spaces::point_t<S>> samples,
                            spaces::scalar_t<S> tolerance = spaces::scalar_t<S>{1e-6}) {
    using std::abs;
    using P = spaces::point_t<S>;
    using V = spaces::tangent_t<S>;

    for (std::size_t i = 0; i < samples.size(); ++i) {
        for (std::size_t j = 0; j < samples.size(); ++j) {
            if (i == j) continue;
            auto v = V(spaces::log_map(space, samples[i], samples[j]));
            auto recovered = P(spaces::exp_map(space, samples[i], v, spaces::scalar_t<S>{1}));

            if constexpr (spaces::Distanced<S>) {
                auto err = spaces::distance(space, samples[j], recovered);
                if (err > tolerance)
                    return VerifyResult::fail("exp(p, log(p,q), 1) != q");
            }
        }
    }
    return VerifyResult::ok();
}

// ── Norm consistency ───────────────────────────────────────────
// norm(v) == sqrt(inner(v, v))  (for InnerProductSpace)

template<InnerProductSpace S>
VerifyResult verify_norm_consistency(const S& space,
                                    std::span<const typename S::VectorType> samples,
                                    typename S::ScalarType tolerance = typename S::ScalarType{1e-8}) {
    using std::abs; using std::sqrt;

    for (const auto& v : samples) {
        auto n = space.norm(v);
        auto from_inner = sqrt(space.inner(v, v));
        if (abs(n - from_inner) > tolerance)
            return VerifyResult::fail("norm(v) != sqrt(<v,v>)");
    }
    return VerifyResult::ok();
}

// ── Convenience overloads (initializer_list) ──────────────────

template<class S>
    requires spaces::Distanced<S>
VerifyResult verify_metric(const S& space,
                           std::initializer_list<spaces::point_t<S>> samples,
                           spaces::scalar_t<S> tolerance = spaces::scalar_t<S>{1e-8}) {
    std::vector<spaces::point_t<S>> v(samples);
    return verify_metric(space, std::span<const spaces::point_t<S>>{v}, tolerance);
}

template<InnerProductSpace S>
VerifyResult verify_inner_product(const S& space,
                                  std::initializer_list<typename S::VectorType> samples,
                                  typename S::ScalarType tolerance = typename S::ScalarType{1e-8}) {
    std::vector<typename S::VectorType> v(samples);
    return verify_inner_product(space, std::span{v}, tolerance);
}

template<class S>
    requires spaces::Exponential<S> && spaces::Logarithmic<S>
VerifyResult verify_exp_log(const S& space,
                            std::initializer_list<spaces::point_t<S>> samples,
                            spaces::scalar_t<S> tolerance = spaces::scalar_t<S>{1e-6}) {
    std::vector<spaces::point_t<S>> v(samples);
    return verify_exp_log(space, std::span<const spaces::point_t<S>>{v}, tolerance);
}

template<InnerProductSpace S>
VerifyResult verify_norm_consistency(const S& space,
                                    std::initializer_list<typename S::VectorType> samples,
                                    typename S::ScalarType tolerance = typename S::ScalarType{1e-8}) {
    std::vector<typename S::VectorType> v(samples);
    return verify_norm_consistency(space, std::span{v}, tolerance);
}


// ── Lorentzian axioms ──────────────────────────────────────────
// A Lorentzian space has no distance, so `verify_metric` does not apply; these
// are the axioms it has instead.
// 1. interval(x, x) == 0, interval(x, y) == interval(y, x)
// 2. causal(x, y) agrees with the sign of interval(x, y)
// 3. the reverse triangle inequality: for x << y << z (each in the timelike
//    future of the one before), proper_time(x, z) >= proper_time(x, y) +
//    proper_time(y, z) -- the longest time between two events is the straight one
// Which is the one that tells a Lorentzian form from one with the wrong
// signature: with two time directions it fails.
template<class S>
    requires LorentzianManifold<S>
VerifyResult verify_lorentzian(const S& space,
                               std::span<const typename S::PointType> samples,
                               typename S::ScalarType tolerance = typename S::ScalarType{1e-9}) {
    using T = typename S::ScalarType;
    using std::abs;
    for (std::size_t i = 0; i < samples.size(); ++i) {
        if (abs(T(space.interval(samples[i], samples[i]))) > tolerance)
            return VerifyResult::fail("interval(x, x) != 0");
        for (std::size_t j = 0; j < samples.size(); ++j) {
            const T sij = space.interval(samples[i], samples[j]);
            if (abs(T(sij - space.interval(samples[j], samples[i]))) > tolerance)
                return VerifyResult::fail("interval is not symmetric");
            const Causal c = space.causal(samples[i], samples[j]);
            if (sij < T{-1e-6} && c != Causal::Timelike) return VerifyResult::fail("negative interval not classed timelike");
            if (sij > T{1e-6} && c != Causal::Spacelike) return VerifyResult::fail("positive interval not classed spacelike");
        }
    }
    for (std::size_t i = 0; i < samples.size(); ++i)
        for (std::size_t j = 0; j < samples.size(); ++j) {
            if (!(space.causal(samples[i], samples[j]) == Causal::Timelike && space.precedes(samples[i], samples[j]))) continue;
            for (std::size_t k = 0; k < samples.size(); ++k) {
                if (!(space.causal(samples[j], samples[k]) == Causal::Timelike && space.precedes(samples[j], samples[k]))) continue;
                const T direct = space.proper_time(samples[i], samples[k]);
                const T via = T(space.proper_time(samples[i], samples[j]) + space.proper_time(samples[j], samples[k]));
                if (direct + tolerance < via)
                    return VerifyResult::fail("reverse triangle inequality fails: a detour through a timelike event is longer");
            }
        }
    return VerifyResult::ok();
}

template<class S>
    requires LorentzianManifold<S>
VerifyResult verify_lorentzian(const S& space,
                               std::initializer_list<typename S::PointType> samples,
                               typename S::ScalarType tolerance = typename S::ScalarType{1e-9}) {
    std::vector<typename S::PointType> v(samples);
    return verify_lorentzian(space, std::span<const typename S::PointType>{v}, tolerance);
}

} // namespace spatium
