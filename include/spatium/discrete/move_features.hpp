#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#include <spatium/discrete/lattice.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#endif

// ── A move as a row of numbers ────────────────────────────────────
//
// A game hands an agent its legal moves; an agent that scores them needs each one as a row of
// numbers. Nothing here knows a game. A row has two blocks. The first is geometric, made here
// from nothing but the piece's vectors and the two cells of the move: how far along how many
// axes it went (the axes in decreasing order of distance, so the row is the same whichever way
// the board is turned or mirrored), how many steps of the piece's vector that took, whether it
// slid, and how far the destination is from the centre of the box. The second block is the
// game's own: what a capture is worth, whether the move gives check, whatever the person
// writing the game wants weighed. `join` puts them side by side. The model that scores rows is
// algebra/models/softmax.hpp.

SPATIUM_EXPORT namespace spatium::discrete {

// ── The geometric block ───────────────────────────────────────────

// How many steps of one of the piece's vectors carry `from` to `to`: 1 for a leap or a slide
// of one step, k for a slide of k, 0 if no step of the piece goes there.
template <std::size_t N>
int steps_taken(const Piece<N>& piece, const Cell<N>& from, const Cell<N>& to) {
    int best = 0;
    for (const Step<N>& step : piece.steps()) {
        int k = 0;
        bool fits = true;
        for (std::size_t i = 0; i < N && fits; ++i) {
            const int d = to[i] - from[i], v = step.vector[i];
            if (v == 0) {
                fits = d == 0;
            } else if (d % v != 0 || d / v <= 0 || (k != 0 && d / v != k)) {
                fits = false;
            } else {
                k = d / v;
            }
        }
        if (!fits || k == 0 || (step.reach == Reach::Leap && k != 1))
            continue;
        best = best == 0 ? k : std::min(best, k);
    }
    return best;
}

template <std::size_t N> inline constexpr std::size_t geometric_feature_count = 2 * N + 3;

// The geometric row of a move, every entry in [0, 1] except the gain in centrality, which is
// in [-1, 1]. Entries 0..N-1 are the distances along each axis, largest first, over the length
// of the axis; N is the number of steps over the longest axis; N+1 is 1 for a slide of more
// than one step; N+2..2N+1 are the destination's distance from the centre along each axis,
// largest first; the last is how much nearer the centre the move went, averaged over the axes.
template <std::size_t N>
std::array<double, geometric_feature_count<N>> geometric_features(
    const Piece<N>& piece, const Shape<N>& shape, const Cell<N>& from, const Cell<N>& to
) {
    std::array<double, geometric_feature_count<N>> out{};
    std::array<double, N> moved{}, centre{};
    double gain = 0, longest = 1;
    for (std::size_t i = 0; i < N; ++i) {
        const double reach = std::max(shape.extent[i] - 1, 1);
        longest = std::max(longest, reach);
        moved[i] = std::abs(to[i] - from[i]) / reach;
        const double c_to = std::abs(2.0 * to[i] - (shape.extent[i] - 1)) / reach;
        const double c_from = std::abs(2.0 * from[i] - (shape.extent[i] - 1)) / reach;
        centre[i] = c_to;
        gain += c_from - c_to;
    }
    std::ranges::sort(moved, std::greater<>{});
    std::ranges::sort(centre, std::greater<>{});
    const int k = steps_taken(piece, from, to);
    for (std::size_t i = 0; i < N; ++i) {
        out[i] = moved[i];
        out[N + 2 + i] = centre[i];
    }
    out[N] = k / longest;
    out[N + 1] = k > 1 ? 1.0 : 0.0;
    out[2 * N + 2] = gain / static_cast<double>(N);
    return out;
}

// Two blocks of a row, side by side.
template <std::size_t A, std::size_t B>
constexpr std::array<double, A + B> join(
    const std::array<double, A>& a, const std::array<double, B>& b
) {
    std::array<double, A + B> out{};
    std::ranges::copy(a, out.begin());
    std::ranges::copy(b, out.begin() + A);
    return out;
}

} // namespace spatium::discrete
