#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <array>
#  include <cmath>
#  include <cstddef>
#  include <cstdint>
#  include <random>
#endif

SPATIUM_EXPORT namespace spatium::inline algebra {

// ── A net of one hidden tanh layer, to a number in (-1, 1) ────────
//
// F inputs, H tanh units, one tanh output: value(x) = tanh(c + sum_j v_j tanh(b_j + w_j . x)).
// F H + 2 H + 1 weights, all in one array in the order w (unit by unit), b, v, c, so a step rule
// of algebra/analysis/optimize.hpp moves them as one vector. The loss is the squared error to a
// target, and `gradient` is its derivative in the weights by the backward pass, written out.
//
// Templated on the scalar like the rest of algebra/: over Dual<T> the same `loss` carries the
// derivative in one weight, which is the second route the backward pass is held to.

template <std::size_t F, std::size_t H, Scalar T = double> class TanhNet {
  public:
    using Input = std::array<T, F>;
    static constexpr std::size_t weight_count = F * H + 2 * H + 1;
    using Weights = std::array<T, weight_count>;

    TanhNet() = default;
    explicit TanhNet(const Weights& w) : w_(w) {}

    // Small random weights: with all equal the units would stay equal, and with all zero the
    // output would not depend on the input. The output layer starts near zero.
    static TanhNet random(std::uint64_t seed) {
        std::mt19937_64 rng(seed);
        std::normal_distribution<double> n(0, 1);
        TanhNet net;
        for (std::size_t j = 0; j < H; ++j) {
            for (std::size_t k = 0; k < F; ++k)
                net.w_[in(j, k)] = T(n(rng) / std::sqrt(static_cast<double>(F)));
            net.w_[bias(j)] = T{0};
            net.w_[out(j)] = T(0.1 * n(rng) / std::sqrt(static_cast<double>(H)));
        }
        net.w_[last()] = T{0};
        return net;
    }

    T value(const Input& x) const {
        std::array<T, H> h;
        return forward(x, h);
    }

    // (value(x) - target)^2
    T loss(const Input& x, T target) const {
        const T e = value(x) - target;
        return e * e;
    }

    // The derivative of `loss` in every weight, by the backward pass.
    Weights gradient(const Input& x, T target) const {
        std::array<T, H> h;
        const T y = forward(x, h);
        const T d_sum = T{2} * (y - target) * (T{1} - y * y);
        Weights g{};
        for (std::size_t j = 0; j < H; ++j) {
            const T d_unit = d_sum * w_[out(j)] * (T{1} - h[j] * h[j]);
            for (std::size_t k = 0; k < F; ++k)
                g[in(j, k)] = d_unit * x[k];
            g[bias(j)] = d_unit;
            g[out(j)] = d_sum * h[j];
        }
        g[last()] = d_sum;
        return g;
    }

    const Weights& weights() const { return w_; }
    Weights& weights() { return w_; }

  private:
    static constexpr std::size_t in(std::size_t j, std::size_t k) { return j * F + k; }
    static constexpr std::size_t bias(std::size_t j) { return F * H + j; }
    static constexpr std::size_t out(std::size_t j) { return F * H + H + j; }
    static constexpr std::size_t last() { return F * H + 2 * H; }

    T forward(const Input& x, std::array<T, H>& h) const {
        using std::tanh;
        T sum = w_[last()];
        for (std::size_t j = 0; j < H; ++j) {
            T pre = w_[bias(j)];
            for (std::size_t k = 0; k < F; ++k)
                pre = pre + w_[in(j, k)] * x[k];
            h[j] = tanh(pre);
            sum = sum + w_[out(j)] * h[j];
        }
        return tanh(sum);
    }

    Weights w_{};
};

} // namespace spatium::algebra
