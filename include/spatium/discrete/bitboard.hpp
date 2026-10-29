#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#include <spatium/algebra/concepts.hpp>

#include <array>
#include <bit>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <iterator>
#endif

// ── The 8x8 board as a discrete space ─────────────────────────────
//
// `Square` is a point of the board, `Bitboard` a set of points held in one
// 64-bit word, and `BoardSymmetryGroup` the eight symmetries of the square
// that carry the board onto itself (the dihedral group D4).
//
// Layout: a1 = 0, b1 = 1, ... h1 = 7, a2 = 8, ... h8 = 63 -- `index = rank * 8
// + file`, so a step north is +8 and a step east is +1.
//
// A set is a word because every question asked of a board (which squares does
// this piece reach, which are occupied, which are attacked twice) is a union,
// an intersection or a complement, and on a word those are one instruction.
// `FiniteSet` answers the same questions on a sorted vector with an allocation
// per operation; the cost of that choice is measured in
// tests/test_bitboard.cpp's header comment.
//
// The symmetry group is exact: a square is addressed by doubled, centred
// coordinates (u, v) = (2 file - 7, 2 rank - 7), both odd integers in [-7, 7],
// so every symmetry is a signed permutation matrix acting on integers with no
// rounding and no case analysis. Composition is the matrix product, the
// inverse the transpose.

SPATIUM_EXPORT namespace spatium::discrete {

inline constexpr int board_size = 8;
inline constexpr int square_count = board_size * board_size;

// ── Square ────────────────────────────────────────────────────────

struct Square {
    std::uint8_t index = 0;

    constexpr Square() = default;
    constexpr explicit Square(int i) : index(static_cast<std::uint8_t>(i)) {}

    // Unchecked: the caller guarantees 0 <= file, rank < board_size.
    static constexpr Square at(int file, int rank) { return Square(rank * board_size + file); }

    constexpr int file() const { return index % board_size; }
    constexpr int rank() const { return index / board_size; }

    constexpr bool operator==(const Square&) const = default;
    constexpr auto operator<=>(const Square&) const = default;
};

constexpr bool on_board(int file, int rank) {
    return 0 <= file && file < board_size && 0 <= rank && rank < board_size;
}

// ── BoardSymmetry: one element of D4 ──────────────────────────────

// The matrix [[a, b], [c, d]] acting on the doubled centred coordinates (u, v).
// Only the eight signed permutation matrices are valid elements; the named
// constructors below are the way to get them.
struct BoardSymmetry {
    std::int8_t a = 1, b = 0, c = 0, d = 1;

    constexpr bool operator==(const BoardSymmetry&) const = default;

    // The image of a square.
    constexpr Square apply(Square s) const {
        const int u = 2 * s.file() - (board_size - 1);
        const int v = 2 * s.rank() - (board_size - 1);
        const int nu = a * u + b * v;
        const int nv = c * u + d * v;
        return Square::at((nu + (board_size - 1)) / 2, (nv + (board_size - 1)) / 2);
    }
};

inline constexpr BoardSymmetry identity_symmetry{1, 0, 0, 1};
inline constexpr BoardSymmetry rotate_90{0, -1, 1, 0}; // quarter turn, counter-clockwise
inline constexpr BoardSymmetry rotate_180{-1, 0, 0, -1};
inline constexpr BoardSymmetry rotate_270{0, 1, -1, 0};
inline constexpr BoardSymmetry mirror_files{-1, 0, 0, 1};    // a <-> h: left-right mirror
inline constexpr BoardSymmetry flip_ranks{1, 0, 0, -1};      // 1 <-> 8: top-bottom flip
inline constexpr BoardSymmetry transpose{0, 1, 1, 0};        // about the a1-h8 diagonal
inline constexpr BoardSymmetry anti_transpose{0, -1, -1, 0}; // about the a8-h1 diagonal

inline constexpr std::array<BoardSymmetry, 8> all_symmetries{identity_symmetry, rotate_90,
                                                             rotate_180,        rotate_270,
                                                             mirror_files,      flip_ranks,
                                                             transpose,         anti_transpose};

// The group D4 in the shape `algebra::power` and the other generic group
// algorithms expect: `compose(a, b)` applies b first, then a.
struct BoardSymmetryGroup {
    using ElementType = BoardSymmetry;

    constexpr ElementType identity() const { return identity_symmetry; }

    constexpr ElementType compose(const ElementType& x, const ElementType& y) const {
        return {
            static_cast<std::int8_t>(x.a * y.a + x.b * y.c),
            static_cast<std::int8_t>(x.a * y.b + x.b * y.d),
            static_cast<std::int8_t>(x.c * y.a + x.d * y.c),
            static_cast<std::int8_t>(x.c * y.b + x.d * y.d)
        };
    }

    // Orthogonal matrices invert by transposition.
    constexpr ElementType inverse(const ElementType& x) const { return {x.a, x.c, x.b, x.d}; }
};
static_assert(algebra::Group<BoardSymmetryGroup>);

// ── Bitboard ──────────────────────────────────────────────────────

class Bitboard {
  public:
    constexpr Bitboard() = default;
    constexpr explicit Bitboard(std::uint64_t bits) : bits_(bits) {}

    static constexpr Bitboard of(Square s) { return Bitboard(std::uint64_t{1} << s.index); }

    constexpr std::uint64_t bits() const { return bits_; }
    constexpr bool empty() const { return bits_ == 0; }
    constexpr int count() const { return std::popcount(bits_); }
    constexpr bool contains(Square s) const { return (bits_ >> s.index) & 1; }

    constexpr Bitboard with(Square s) const { return Bitboard(bits_ | of(s).bits_); }
    constexpr Bitboard without(Square s) const { return Bitboard(bits_ & ~of(s).bits_); }

    // The square with the smallest / largest index. Precondition: not empty.
    constexpr Square lowest() const { return Square(std::countr_zero(bits_)); }
    constexpr Square highest() const { return Square(63 - std::countl_zero(bits_)); }

    constexpr Bitboard operator&(Bitboard o) const { return Bitboard(bits_ & o.bits_); }
    constexpr Bitboard operator|(Bitboard o) const { return Bitboard(bits_ | o.bits_); }
    constexpr Bitboard operator^(Bitboard o) const { return Bitboard(bits_ ^ o.bits_); }
    constexpr Bitboard operator~() const { return Bitboard(~bits_); }
    constexpr Bitboard& operator&=(Bitboard o) {
        bits_ &= o.bits_;
        return *this;
    }
    constexpr Bitboard& operator|=(Bitboard o) {
        bits_ |= o.bits_;
        return *this;
    }
    constexpr Bitboard& operator^=(Bitboard o) {
        bits_ ^= o.bits_;
        return *this;
    }
    constexpr bool operator==(const Bitboard&) const = default;

    // Squares from the lowest index to the highest, without copying anything:
    // each step clears the lowest set bit of a local word.
    class Iterator {
      public:
        using value_type = Square;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::input_iterator_tag;

        constexpr Iterator() = default;
        constexpr explicit Iterator(std::uint64_t bits) : bits_(bits) {}
        constexpr Square operator*() const { return Square(std::countr_zero(bits_)); }
        constexpr Iterator& operator++() {
            bits_ &= bits_ - 1;
            return *this;
        }
        constexpr Iterator operator++(int) {
            Iterator t = *this;
            ++*this;
            return t;
        }
        constexpr bool operator==(const Iterator&) const = default;

      private:
        std::uint64_t bits_ = 0;
    };

    constexpr Iterator begin() const { return Iterator(bits_); }
    constexpr Iterator end() const { return Iterator(0); }

    // ── Symmetries of the board acting on a set of squares ────────

    // Top-bottom flip: rank r goes to rank 7 - r. One byte swap.
    constexpr Bitboard flipped_ranks() const { return Bitboard(std::byteswap(bits_)); }

    // Left-right mirror: file f goes to file 7 - f.
    constexpr Bitboard mirrored_files() const {
        constexpr std::uint64_t k1 = 0x5555555555555555ull, k2 = 0x3333333333333333ull,
                                k4 = 0x0f0f0f0f0f0f0f0full;
        std::uint64_t x = bits_;
        x = ((x >> 1) & k1) | ((x & k1) << 1);
        x = ((x >> 2) & k2) | ((x & k2) << 2);
        x = ((x >> 4) & k4) | ((x & k4) << 4);
        return Bitboard(x);
    }

    // Transpose about the a1-h8 diagonal: (file, rank) goes to (rank, file).
    constexpr Bitboard transposed() const {
        constexpr std::uint64_t k1 = 0x5500550055005500ull, k2 = 0x3333000033330000ull,
                                k4 = 0x0f0f0f0f00000000ull;
        std::uint64_t x = bits_;
        std::uint64_t t = k4 & (x ^ (x << 28));
        x ^= t ^ (t >> 28);
        t = k2 & (x ^ (x << 14));
        x ^= t ^ (t >> 14);
        t = k1 & (x ^ (x << 7));
        x ^= t ^ (t >> 7);
        return Bitboard(x);
    }

    // Any element of D4, square by square. The fast forms above must agree
    // with it; it is the definition they are tested against.
    constexpr Bitboard transformed(const BoardSymmetry& g) const {
        Bitboard out;
        for (Square s : *this)
            out = out.with(g.apply(s));
        return out;
    }

  private:
    std::uint64_t bits_ = 0;
};

}  // namespace spatium::discrete
