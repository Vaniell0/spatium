#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <array>
#  include <cmath>
#  include <cstddef>
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

} // namespace spatium::algebra
