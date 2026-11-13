#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <array>
#  include <cmath>
#  include <cstddef>
#  include <cstdint>
#endif

SPATIUM_EXPORT namespace spatium::inline algebra {

// ── Step rules ────────────────────────────────────────────────────
//
// A step rule moves parameters given a gradient. It is deterministic in the gradient: where the
// gradient comes from -- one sample, a batch, the whole data -- is the caller's. The models of
// algebra/models/ give gradients; the rules here take them.

// AdaGrad (Duchi, Hazan and Singer 2011): each parameter moves by rate * g / sqrt(G), with G the
// sum of its squared gradients so far, so a parameter that has been informative settles and a
// rarely seen one can still move. The first step moves every parameter with a nonzero gradient
// by exactly the rate, whatever the gradient's scale; scaling every gradient by c > 0 changes
// no step at all.
template <std::size_t N, Scalar T = double> struct AdaGrad {
    std::array<T, N> squared{};

    void step(std::array<T, N>& params, const std::array<T, N>& gradient, T rate) {
        using std::sqrt;
        for (std::size_t j = 0; j < N; ++j) {
            squared[j] = squared[j] + gradient[j] * gradient[j];
            if (squared[j] > T{0})
                params[j] = params[j] - rate * gradient[j] / sqrt(squared[j]);
        }
    }
};

// Adam (Kingma and Ba 2015): a running mean m of the gradient and v of its square, each
// corrected for starting at zero, and a step of rate * m / (sqrt(v) + epsilon). The first step
// is rate * g / (|g| + epsilon): the rate, signed, for any gradient well above epsilon; scaling
// every gradient by c changes the steps only through epsilon.
template <std::size_t N, Scalar T = double> struct Adam {
    T beta1{0.9}, beta2{0.999}, epsilon{1e-8};
    std::array<T, N> mean{}, square{};
    std::uint64_t steps = 0;

    void step(std::array<T, N>& params, const std::array<T, N>& gradient, T rate) {
        using std::pow;
        using std::sqrt;
        ++steps;
        const T k = T(static_cast<double>(steps));
        const T c1 = T{1} - pow(beta1, k), c2 = T{1} - pow(beta2, k);
        for (std::size_t j = 0; j < N; ++j) {
            mean[j] = beta1 * mean[j] + (T{1} - beta1) * gradient[j];
            square[j] = beta2 * square[j] + (T{1} - beta2) * gradient[j] * gradient[j];
            params[j] = params[j] - rate * (mean[j] / c1) / (sqrt(square[j] / c2) + epsilon);
        }
    }
};

// Muon (Keller Jordan 2024, github.com/KellerJordan/Muon, `muon.py`, which this follows): for a
// matrix of parameters, the step is the gradient's momentum made orthogonal -- its singular values
// moved to about 1, its directions kept -- so a weight matrix moves as much along a rarely seen
// direction as along a dominant one. The first R * C parameters are the matrix, R rows of C in order
// (a `TanhNet`'s hidden layer, H by F); the rest, a net's biases and output, move by Adam, as the
// reference's auxiliary Adam does.
//
// As in the reference: the momentum is a running mean, buf <- beta buf + (1 - beta) g with beta =
// 0.95, and Nesterov's update (1 - beta) g + beta buf; the orthogonal factor is approached by five
// steps of the quintic Newton-Schulz iteration X <- a X + (b A + c A^2) X, A = X X^T, (a, b, c) =
// (3.4445, -4.7750, 2.0315), from the update over its Frobenius norm plus 1e-7, worked on the wide
// orientation; it leaves the singular values spread over about (0.5, 1.5), which the reference says
// does not hurt; and the step is the rate times that, times max(1, R / C)^(1/2). The rate is in units
// of the spectral norm per step (the reference's default 0.02); the rest moves by Adam at
// `rest_rate` times it, with the reference's betas (0.9, 0.95) and epsilon 1e-10 (its defaults: 3e-4
// against 0.02). Scaling every gradient by c > 0 changes no step of the matrix but through 1e-7.
template <std::size_t N, std::size_t R, std::size_t C, Scalar T = double> struct Muon {
    static_assert(R * C <= N && R > 0 && C > 0, "the matrix is the first R * C of the parameters");
    T momentum{0.95}, epsilon{1e-7}, rest_rate{0.015};
    std::array<T, R * C> buffer{};
    Adam<N - R * C, T> rest{T{0.9}, T{0.95}, T{1e-10}};

    void step(std::array<T, N>& params, const std::array<T, N>& gradient, T rate) {
        using std::sqrt;
        std::array<T, R * C> u;
        for (std::size_t i = 0; i < R * C; ++i) {
            buffer[i] = momentum * buffer[i] + (T{1} - momentum) * gradient[i];
            u[i] = (T{1} - momentum) * gradient[i] + momentum * buffer[i];
        }
        const std::array<T, R * C> o = orthogonalized(u, epsilon);
        const T scale = R > C ? sqrt(T(static_cast<double>(R) / static_cast<double>(C))) : T{1};
        for (std::size_t i = 0; i < R * C; ++i) params[i] = params[i] - rate * scale * o[i];
        if constexpr (N > R * C) {
            std::array<T, N - R * C> p, d;
            for (std::size_t i = 0; i < N - R * C; ++i) p[i] = params[R * C + i], d[i] = gradient[R * C + i];
            rest.step(p, d, rate * rest_rate);
            for (std::size_t i = 0; i < N - R * C; ++i) params[R * C + i] = p[i];
        }
    }

    // The orthogonal factor of the R by C matrix g, approached; worked as the wide m by n one, m <= n,
    // so the products are m by m.
    static std::array<T, R * C> orthogonalized(const std::array<T, R * C>& g, T epsilon) {
        using std::sqrt;
        constexpr bool tall = R > C;
        constexpr std::size_t m = tall ? C : R, n = tall ? R : C;
        std::array<T, m * n> x;
        T norm{0};
        for (std::size_t i = 0; i < R; ++i)
            for (std::size_t j = 0; j < C; ++j) {
                const T v = g[i * C + j];
                x[tall ? j * n + i : i * n + j] = v;
                norm = norm + v * v;
            }
        norm = sqrt(norm) + epsilon;
        for (T& v : x) v = v / norm;
        const T a{3.4445}, b{-4.7750}, c{2.0315};
        for (int k = 0; k < 5; ++k) {
            std::array<T, m * m> A{}, B{};
            for (std::size_t i = 0; i < m; ++i)
                for (std::size_t j = 0; j < m; ++j)
                    for (std::size_t t = 0; t < n; ++t) A[i * m + j] = A[i * m + j] + x[i * n + t] * x[j * n + t];
            for (std::size_t i = 0; i < m; ++i)
                for (std::size_t j = 0; j < m; ++j) {
                    T aa{0};
                    for (std::size_t t = 0; t < m; ++t) aa = aa + A[i * m + t] * A[t * m + j];
                    B[i * m + j] = b * A[i * m + j] + c * aa;
                }
            std::array<T, m * n> next;
            for (std::size_t i = 0; i < m; ++i)
                for (std::size_t j = 0; j < n; ++j) {
                    T bx{0};
                    for (std::size_t t = 0; t < m; ++t) bx = bx + B[i * m + t] * x[t * n + j];
                    next[i * n + j] = a * x[i * n + j] + bx;
                }
            x = next;
        }
        std::array<T, R * C> out;
        for (std::size_t i = 0; i < R; ++i)
            for (std::size_t j = 0; j < C; ++j) out[i * C + j] = x[tall ? j * n + i : i * n + j];
        return out;
    }
};

} // namespace spatium::algebra
