#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
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

template<MetricSpace S>
VerifyResult verify_metric(const S& space,
                           std::span<const typename S::PointType> samples,
                           typename S::ScalarType tolerance = typename S::ScalarType{1e-8}) {
    using T = typename S::ScalarType;
    using std::abs;

    for (std::size_t i = 0; i < samples.size(); ++i) {
        // Non-negativity
        auto d_ii = space.distance(samples[i], samples[i]);
        if (d_ii < T{0} || abs(d_ii) > tolerance)
            return VerifyResult::fail("d(x,x) != 0");

        for (std::size_t j = i + 1; j < samples.size(); ++j) {
            auto d_ij = space.distance(samples[i], samples[j]);
            auto d_ji = space.distance(samples[j], samples[i]);

            // Non-negativity
            if (d_ij < -tolerance)
                return VerifyResult::fail("d(x,y) < 0");

            // Symmetry
            if (abs(d_ij - d_ji) > tolerance)
                return VerifyResult::fail("d(x,y) != d(y,x)");

            // Triangle inequality
            for (std::size_t k = 0; k < samples.size(); ++k) {
                auto d_ik = space.distance(samples[i], samples[k]);
                auto d_kj = space.distance(samples[k], samples[j]);
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

template<Manifold S>
VerifyResult verify_exp_log(const S& space,
                            std::span<const typename S::PointType> samples,
                            typename S::ScalarType tolerance = typename S::ScalarType{1e-6}) {
    using std::abs;

    for (std::size_t i = 0; i < samples.size(); ++i) {
        for (std::size_t j = 0; j < samples.size(); ++j) {
            if (i == j) continue;
            auto v = space.log_map(samples[i], samples[j]);
            auto recovered = space.exp_map(samples[i], v, typename S::ScalarType{1});

            if constexpr (MetricSpace<S>) {
                auto err = space.distance(samples[j], recovered);
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

template<MetricSpace S>
VerifyResult verify_metric(const S& space,
                           std::initializer_list<typename S::PointType> samples,
                           typename S::ScalarType tolerance = typename S::ScalarType{1e-8}) {
    std::vector<typename S::PointType> v(samples);
    return verify_metric(space, std::span{v}, tolerance);
}

template<InnerProductSpace S>
VerifyResult verify_inner_product(const S& space,
                                  std::initializer_list<typename S::VectorType> samples,
                                  typename S::ScalarType tolerance = typename S::ScalarType{1e-8}) {
    std::vector<typename S::VectorType> v(samples);
    return verify_inner_product(space, std::span{v}, tolerance);
}

template<Manifold S>
VerifyResult verify_exp_log(const S& space,
                            std::initializer_list<typename S::PointType> samples,
                            typename S::ScalarType tolerance = typename S::ScalarType{1e-6}) {
    std::vector<typename S::PointType> v(samples);
    return verify_exp_log(space, std::span{v}, tolerance);
}

template<InnerProductSpace S>
VerifyResult verify_norm_consistency(const S& space,
                                    std::initializer_list<typename S::VectorType> samples,
                                    typename S::ScalarType tolerance = typename S::ScalarType{1e-8}) {
    std::vector<typename S::VectorType> v(samples);
    return verify_norm_consistency(space, std::span{v}, tolerance);
}

} // namespace spatium
