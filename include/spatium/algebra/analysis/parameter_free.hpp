#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <algorithm>
#  include <cmath>
#  include <cstddef>
#  include <span>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::inline algebra {

// ── Step rules with no rate: the distance travelled over the gradients seen ─
//
// DoG (Ivgi, Hinder and Carmon 2023): the step at t is the farthest the parameters have gone from
// where they started, over the root of the sum of the squared gradient norms so far,
//
//   eta_t = max_{i <= t} |x_i - x_0| / sqrt(sum_{i <= t} |g_i|^2),
//
// so it starts tiny (the first distance is r_eps = alpha (1 + |x_0|)) and grows while the parameters
// keep going one way, settling near the best fixed rate without being told one. `LayerDoG` is their
// per-layer form, L-DoG: each block of the parameters (a matrix, a vector) has its own distance and
// its own sum, with 1e-8 added under the root. The rate the caller passes multiplies the step: 1 is
// the paper's choice, about 0.5 to 1.5 works, above that it diverges.
//
// `MuonDoG` joins it with Muon, which the paper leaves open (DoG with momentum): each matrix block
// moves along its orthogonalised momentum O (Newton-Schulz, as `Muon`), by the block's distance over
// the root of the sum of |O|^2 seen, so its step needs no rate either and grows with the distance
// travelled; the other blocks move by L-DoG on their gradients. This combination is ours, untested by
// the paper.
//
// Both take their shape when they run, not when they compile: the parameters and the gradient are
// spans, the blocks are offsets into them, and the state is sized on the first step. So one rule moves
// a fixed `std::array` net and a model whose sizes are read from a file alike. The orthogonal factor is
// the same five Newton-Schulz steps as `Muon<N, R, C>::orthogonalized`, at a shape given at run time.
//
// `PolyAverage` is the running weighted mean the paper evaluates through (polynomial-decay averaging,
// gamma 8): x_bar_t = (1 - w) x_bar_{t-1} + w x_t with w = (gamma + 1) / (t + gamma).

// A block of the parameters: `rows` by `cols` from `offset`; a vector is one row.
struct ParamBlock {
    std::size_t offset = 0, rows = 1, cols = 1;
    std::size_t size() const { return rows * cols; }
    bool matrix() const { return rows > 1 && cols > 1; }
};

namespace parameter_free_detail {

template <class T> T norm2(std::span<const T> x) {
    T s{0};
    for (const T& v : x)
        s = s + v * v;
    return s;
}

// The orthogonal factor of the R by C matrix g (row-major), approached by five Newton-Schulz steps
// with the coefficients `Muon` uses; worked as the wide one, so the products are min(R, C) square.
template <class T> std::vector<T> orthogonalized(std::span<const T> g, std::size_t R, std::size_t C, T epsilon) {
    using std::sqrt;
    const bool tall = R > C;
    const std::size_t m = tall ? C : R, n = tall ? R : C;
    std::vector<T> x(m * n);
    T norm{0};
    for (std::size_t i = 0; i < R; ++i)
        for (std::size_t j = 0; j < C; ++j) {
            const T v = g[i * C + j];
            x[tall ? j * n + i : i * n + j] = v;
            norm = norm + v * v;
        }
    norm = sqrt(norm) + epsilon;
    for (T& v : x)
        v = v / norm;
    const T a{3.4445}, b{-4.7750}, c{2.0315};
    std::vector<T> A(m * m), B(m * m), next(m * n);
    for (int k = 0; k < 5; ++k) {
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) {
                T s{0};
                for (std::size_t t = 0; t < n; ++t)
                    s = s + x[i * n + t] * x[j * n + t];
                A[i * m + j] = s;
            }
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) {
                T aa{0};
                for (std::size_t t = 0; t < m; ++t)
                    aa = aa + A[i * m + t] * A[t * m + j];
                B[i * m + j] = b * A[i * m + j] + c * aa;
            }
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                T bx{0};
                for (std::size_t t = 0; t < m; ++t)
                    bx = bx + B[i * m + t] * x[t * n + j];
                next[i * n + j] = a * x[i * n + j] + bx;
            }
        x.swap(next);
    }
    std::vector<T> out(R * C);
    for (std::size_t i = 0; i < R; ++i)
        for (std::size_t j = 0; j < C; ++j)
            out[i * C + j] = x[tall ? j * n + i : i * n + j];
    return out;
}

// The blocks as given, or the whole array as one when none are.
inline std::vector<ParamBlock> blocks_or_all(const std::vector<ParamBlock>& b, std::size_t n) {
    return b.empty() ? std::vector<ParamBlock>{{0, 1, n}} : b;
}

} // namespace parameter_free_detail

template <Scalar T = double> struct LayerDoG {
    std::vector<ParamBlock> blocks; // empty: one block, the whole array
    T alpha{1e-8}, epsilon{1e-8};
    std::vector<T> start;
    std::vector<T> reach, seen; // per block: the farthest distance from the start, the sum of |g|^2
    bool started = false;

    void step(std::span<T> params, std::span<const T> gradient, T scale = T{1}) {
        using std::sqrt;
        using namespace parameter_free_detail;
        const auto bs = blocks_or_all(blocks, params.size());
        if (!started) {
            start.assign(params.begin(), params.end()), started = true;
            reach.assign(bs.size(), T{0}), seen.assign(bs.size(), T{0});
            for (std::size_t k = 0; k < bs.size(); ++k)
                reach[k] = alpha * (T{1} + sqrt(norm2<T>(std::span<const T>(params.data() + bs[k].offset, bs[k].size()))));
        }
        for (std::size_t k = 0; k < bs.size(); ++k) {
            const ParamBlock& b = bs[k];
            T moved{0};
            for (std::size_t i = b.offset; i < b.offset + b.size(); ++i)
                moved = moved + (params[i] - start[i]) * (params[i] - start[i]);
            reach[k] = std::max(reach[k], sqrt(moved));
            seen[k] = seen[k] + norm2<T>(std::span<const T>(gradient.data() + b.offset, b.size()));
            const T eta = scale * reach[k] / sqrt(seen[k] + epsilon);
            for (std::size_t i = b.offset; i < b.offset + b.size(); ++i)
                params[i] = params[i] - eta * gradient[i];
        }
    }
};

template <Scalar T = double> struct MuonDoG {
    std::vector<ParamBlock> blocks; // the matrices among them move by Muon, the rest by L-DoG
    T momentum{0.95}, alpha{1e-8}, epsilon{1e-8};
    std::vector<T> start, buffer;
    std::vector<T> reach, seen;
    bool started = false;

    void step(std::span<T> params, std::span<const T> gradient, T scale = T{1}) {
        using std::sqrt;
        using namespace parameter_free_detail;
        const auto bs = blocks_or_all(blocks, params.size());
        if (!started) {
            start.assign(params.begin(), params.end()), started = true;
            buffer.assign(params.size(), T{0});
            reach.assign(bs.size(), T{0}), seen.assign(bs.size(), T{0});
            for (std::size_t k = 0; k < bs.size(); ++k)
                reach[k] = alpha * (T{1} + sqrt(norm2<T>(std::span<const T>(params.data() + bs[k].offset, bs[k].size()))));
        }
        for (std::size_t k = 0; k < bs.size(); ++k) {
            const ParamBlock& b = bs[k];
            T moved{0};
            for (std::size_t i = b.offset; i < b.offset + b.size(); ++i)
                moved = moved + (params[i] - start[i]) * (params[i] - start[i]);
            reach[k] = std::max(reach[k], sqrt(moved));
            if (!b.matrix()) {
                seen[k] = seen[k] + norm2<T>(std::span<const T>(gradient.data() + b.offset, b.size()));
                const T eta = scale * reach[k] / sqrt(seen[k] + epsilon);
                for (std::size_t i = b.offset; i < b.offset + b.size(); ++i)
                    params[i] = params[i] - eta * gradient[i];
                continue;
            }
            std::vector<T> u(b.size());
            for (std::size_t i = 0; i < b.size(); ++i) {
                T& m = buffer[b.offset + i];
                const T g = gradient[b.offset + i];
                m = momentum * m + (T{1} - momentum) * g;
                u[i] = (T{1} - momentum) * g + momentum * m;
            }
            const auto o = orthogonalized<T>(std::span<const T>(u), b.rows, b.cols, epsilon);
            const T shape = b.rows > b.cols ? sqrt(T(double(b.rows) / double(b.cols))) : T{1};
            seen[k] = seen[k] + shape * shape * norm2<T>(std::span<const T>(o));
            const T eta = scale * reach[k] / sqrt(seen[k] + epsilon);
            for (std::size_t i = 0; i < b.size(); ++i)
                params[b.offset + i] = params[b.offset + i] - eta * shape * o[i];
        }
    }
};

template <Scalar T = double> struct PolyAverage {
    T gamma{8};
    std::vector<T> mean;
    long count = 0;

    void add(std::span<const T> x) {
        if (mean.size() != x.size())
            mean.assign(x.size(), T{0}), count = 0;
        ++count;
        const T w = (gamma + T{1}) / (T(double(count)) + gamma);
        for (std::size_t i = 0; i < x.size(); ++i)
            mean[i] = (T{1} - w) * mean[i] + w * x[i];
    }
};

} // namespace spatium::algebra
