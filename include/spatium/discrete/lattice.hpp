#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#include <spatium/algebra/groups/action.hpp>
#include <spatium/algebra/groups/hyperoctahedral.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <vector>
#endif

// ── A board of any number of axes, and pieces as vectors ──────────
//
// `Bitboard` is the 8x8 board in one word. This is the same question without
// the 8x8: a box of cells with N axes of any length, and a piece that is a set
// of vectors and nothing else.
//
// A piece is described by the shape of its step, not by a list of squares. The
// step (1, 2) is the knight because its orbit under the symmetries of the
// board -- every signed permutation of the axes, the hyperoctahedral group,
// N! 2^N elements -- is the eight knight moves. The same description in four
// axes (file, rank, turn, timeline) gives the rook along the axes, the bishop
// along the diagonals of any two, the unicorn along three and the dragon along
// four, with no case written for any of them:
//
//   orbit of (1, 0, 0, 0)   8 vectors   rook
//   orbit of (1, 1, 0, 0)  24           bishop
//   orbit of (1, 1, 1, 0)  32           unicorn
//   orbit of (1, 1, 1, 1)  16           dragon
//   orbit of (1, 2, 0, 0)  48           knight
//
// A piece reaches by leaping (the vector once) or by sliding (the vector over
// and over, to the edge of the box or the first piece in the way). The slide
// is a sum over k = 1, 2, ... of from + k v, cut where the box or a blocker
// ends it; the box is the cutoff, so nothing infinite is ever added.
//
// A piece with a direction -- the pawn -- is not carried to itself by the
// group, so it is not an orbit; it is built from `Step`s by hand and is the
// case the symmetry tests must reject.
//
// This is the definition: it walks the board. The word-parallel tables of
// attacks.hpp are the fast route to the same sets on 8x8, and are held to
// this one.

SPATIUM_EXPORT namespace spatium::discrete {

template <std::size_t N> using Cell = algebra::IntVec<N>;

template <std::size_t N> constexpr Cell<N> operator+(const Cell<N>& a, const Cell<N>& b) {
    Cell<N> out{};
    for (std::size_t i = 0; i < N; ++i)
        out[i] = a[i] + b[i];
    return out;
}

// ── The symmetries of a box ───────────────────────────────────────

// The signed permutations of the axes (algebra/groups/hyperoctahedral.hpp) act on the steps
// of a piece and on the cells of a cube; orbits and stabilisers are algebra::orbit and
// algebra::stabilizer over those two actions.
template <std::size_t N> using Symmetry = algebra::SignedPermutation<N>;

// ── The box ───────────────────────────────────────────────────────

template <std::size_t N> struct Shape {
    std::array<int, N> extent{};

    constexpr bool contains(const Cell<N>& c) const {
        for (std::size_t i = 0; i < N; ++i)
            if (c[i] < 0 || c[i] >= extent[i])
                return false;
        return true;
    }

    constexpr std::size_t size() const {
        std::size_t n = 1;
        for (int e : extent)
            n *= static_cast<std::size_t>(e);
        return n;
    }

    constexpr bool is_cube() const {
        for (std::size_t i = 1; i < N; ++i)
            if (extent[i] != extent[0])
                return false;
        return true;
    }

    // Position of a cell in a dense array, the first axis varying fastest.
    constexpr std::size_t index(const Cell<N>& c) const {
        std::size_t at = 0, stride = 1;
        for (std::size_t i = 0; i < N; ++i) {
            at += static_cast<std::size_t>(c[i]) * stride;
            stride *= static_cast<std::size_t>(extent[i]);
        }
        return at;
    }

    constexpr Cell<N> cell_at(std::size_t at) const {
        Cell<N> c{};
        for (std::size_t i = 0; i < N; ++i) {
            c[i] = static_cast<int>(at % static_cast<std::size_t>(extent[i]));
            at /= static_cast<std::size_t>(extent[i]);
        }
        return c;
    }

    // Where the symmetry carries a cell. Only a box with equal axes is carried onto itself by
    // every signed permutation, so this is for a cube.
    constexpr Cell<N> transformed(const Symmetry<N>& s, const Cell<N>& c) const {
        return algebra::CubeCellAction<N>{extent[0]}.act(s, c);
    }
};

// ── Pieces ────────────────────────────────────────────────────────

enum class Reach : std::uint8_t { Leap, Slide };

template <std::size_t N> struct Step {
    Cell<N> vector{};
    Reach reach = Reach::Leap;

    bool operator==(const Step&) const = default;
};

template <std::size_t N> class Piece {
  public:
    Piece() = default;

    // Every image of `v` under the symmetries of the cube, as steps of the given reach.
    static Piece orbit(const Cell<N>& v, Reach reach) {
        Piece p;
        for (const Cell<N>& w :
             algebra::orbit(algebra::HyperoctahedralGroup<N>{}, algebra::SignedPermutationAction<N>{}, v))
            p.steps_.push_back({w, reach});
        p.normalise();
        return p;
    }

    // A piece from steps given by hand, for the ones with a direction.
    static Piece of(std::vector<Step<N>> steps) {
        Piece p;
        p.steps_ = std::move(steps);
        p.normalise();
        return p;
    }

    // Both pieces' moves. A vector in both stays a slide if either slides.
    Piece operator+(const Piece& other) const {
        Piece p = *this;
        p.steps_.insert(p.steps_.end(), other.steps_.begin(), other.steps_.end());
        p.normalise();
        return p;
    }

    const std::vector<Step<N>>& steps() const { return steps_; }

    // The number of distinct vectors: how many directions or jumps the piece has.
    std::size_t size() const { return steps_.size(); }

  private:
    // Sorted by vector, one step per vector, a slide winning over a leap.
    void normalise() {
        std::ranges::sort(steps_, [](const Step<N>& a, const Step<N>& b) {
            return a.vector != b.vector ? a.vector < b.vector : a.reach > b.reach;
        });
        const auto dup = std::ranges::unique(steps_, [](const Step<N>& a, const Step<N>& b) {
            return a.vector == b.vector;
        });
        steps_.erase(dup.begin(), dup.end());
    }

    std::vector<Step<N>> steps_;
};

// The step with k ones on the first k axes: (1,0,..) rook, (1,1,0,..) bishop, ...
template <std::size_t N> Piece<N> diagonal_mover(std::size_t k, Reach reach) {
    Cell<N> v{};
    for (std::size_t i = 0; i < k && i < N; ++i)
        v[i] = 1;
    return Piece<N>::orbit(v, reach);
}

template <std::size_t N> Piece<N> rook() {
    return diagonal_mover<N>(1, Reach::Slide);
}
template <std::size_t N> Piece<N> bishop() {
    static_assert(N >= 2, "a bishop moves along two axes at once");
    return diagonal_mover<N>(2, Reach::Slide);
}
template <std::size_t N> Piece<N> unicorn() {
    static_assert(N >= 3, "a unicorn moves along three axes at once");
    return diagonal_mover<N>(3, Reach::Slide);
}
template <std::size_t N> Piece<N> dragon() {
    static_assert(N >= 4, "a dragon moves along four axes at once");
    return diagonal_mover<N>(4, Reach::Slide);
}

template <std::size_t N> Piece<N> knight() {
    static_assert(N >= 2, "a knight needs two axes");
    Cell<N> v{};
    v[0] = 1;
    v[1] = 2;
    return Piece<N>::orbit(v, Reach::Leap);
}

// One step along any set of axes, or a slide along any set of axes: every other cell of the
// 3^N neighbourhood of a cell, or the lines through it.
template <std::size_t N> Piece<N> king() {
    Piece<N> p;
    for (std::size_t k = 1; k <= N; ++k)
        p = p + diagonal_mover<N>(k, Reach::Leap);
    return p;
}
template <std::size_t N> Piece<N> queen() {
    Piece<N> p;
    for (std::size_t k = 1; k <= N; ++k)
        p = p + diagonal_mover<N>(k, Reach::Slide);
    return p;
}

// ── Where a piece reaches ─────────────────────────────────────────

// Calls `visit` on each cell the piece reaches from `from`, a blocker included. `occupied`
// says whether a cell holds a piece.
template <std::size_t N, class Occupied, class Visit>
void for_each_attack(
    const Piece<N>& piece, const Shape<N>& shape, const Cell<N>& from, Occupied&& occupied,
    Visit&& visit
) {
    for (const Step<N>& step : piece.steps())
        for (Cell<N> c = from + step.vector; shape.contains(c); c = c + step.vector) {
            visit(c);
            if (step.reach == Reach::Leap || occupied(c))
                break;
        }
}

// The same, as a sorted list of cells.
template <std::size_t N, class Occupied>
std::vector<Cell<N>> attacks(
    const Piece<N>& piece, const Shape<N>& shape, const Cell<N>& from, Occupied&& occupied
) {
    std::vector<Cell<N>> out;
    std::size_t bound = 0; // a leap reaches one cell, a slide no more than the longest axis
    for (const Step<N>& step : piece.steps())
        bound += step.reach == Reach::Leap ? 1 : static_cast<std::size_t>(*std::ranges::max_element(shape.extent));
    out.reserve(bound);
    for_each_attack(piece, shape, from, occupied, [&](const Cell<N>& c) { out.push_back(c); });
    std::ranges::sort(out);
    out.erase(std::ranges::unique(out).begin(), out.end());
    return out;
}

} // namespace spatium::discrete
