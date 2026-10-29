#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#include <spatium/discrete/bitboard.hpp>

#include <array>
#include <cstdint>
#include <utility>
#endif

// ── Which squares a piece reaches ─────────────────────────────────
//
// A piece's one-move set on the empty board is a pattern of vectors; what makes
// it depend on the position is only the edge of the board and, for a piece that
// slides, the pieces in the way. So the work splits the same way:
//
//   leapers   (knight, king, ...)  the set is a table, one word per square,
//             built at compile time from the orbit of a vector (a, b) under D4
//             -- the board edge is cut once, when the table is built;
//   sliders   (rook, bishop, queen) the set is a ray per direction, cut at the
//             first piece in the way: ray(from) ^ ray(first blocker).
//
// `attacks` are the squares a piece reaches, a blocker included. Removing the
// piece's own side from the result is the caller's mask, not part of the
// geometry.

SPATIUM_EXPORT namespace spatium::discrete {

// ── Leapers ───────────────────────────────────────────────────────

namespace detail {

constexpr std::array<Bitboard, square_count> make_leaper_table(int a, int b) {
    std::array<Bitboard, square_count> table{};
    for (int s = 0; s < square_count; ++s) {
        const Square from(s);
        std::uint64_t bits = 0;
        for (int image = 0; image < 8; ++image) { // signs and the swap of the two axes
            int u = (image & 1) ? -a : a;
            int v = (image & 2) ? -b : b;
            if (image & 4)
                std::swap(u, v);
            const int file = from.file() + u, rank = from.rank() + v;
            if (on_board(file, rank))
                bits |= Bitboard::of(Square::at(file, rank)).bits();
        }
        table[s] = Bitboard(bits);
    }
    return table;
}

template <int A, int B>
inline constexpr std::array<Bitboard, square_count> leaper_table = make_leaper_table(A, B);

} // namespace detail

// Every image of the vector (A, B) under the symmetries of the square:
// (1,2) is the knight, (1,0) the one-step orthogonal, (1,1) the one-step diagonal,
// (1,3) the camel. The pieces below are named instances of this.
template <int A, int B> constexpr Bitboard leaper_attacks(Square from) {
    static_assert(A != 0 || B != 0, "the zero vector is not a move");
    return detail::leaper_table<A, B>[from.index];
}

constexpr Bitboard knight_attacks(Square from) {
    return leaper_attacks<1, 2>(from);
}
constexpr Bitboard king_attacks(Square from) {
    return leaper_attacks<1, 0>(from) | leaper_attacks<1, 1>(from);
}

// ── Sliders ───────────────────────────────────────────────────────

enum class Direction : std::uint8_t {
    North,
    NorthEast,
    East,
    SouthEast,
    South,
    SouthWest,
    West,
    NorthWest
};

namespace detail {

inline constexpr std::array<std::pair<int, int>, 8> direction_step{
    {{0, 1}, {1, 1}, {1, 0}, {1, -1}, {0, -1}, {-1, -1}, {-1, 0}, {-1, 1}}
};

// A step in these directions raises the square index, so the nearest blocker
// is the lowest set bit; in the others it is the highest.
constexpr bool raises_index(Direction d) {
    return d == Direction::North || d == Direction::NorthEast || d == Direction::East ||
           d == Direction::NorthWest;
}

constexpr std::array<std::array<Bitboard, square_count>, 8> make_ray_table() {
    std::array<std::array<Bitboard, square_count>, 8> table{};
    for (int d = 0; d < 8; ++d)
        for (int s = 0; s < square_count; ++s) {
            const auto [df, dr] = direction_step[d];
            std::uint64_t bits = 0;
            for (int file = Square(s).file() + df, rank = Square(s).rank() + dr;
                 on_board(file, rank); file += df, rank += dr)
                bits |= Bitboard::of(Square::at(file, rank)).bits();
            table[d][s] = Bitboard(bits);
        }
    return table;
}

inline constexpr std::array<std::array<Bitboard, square_count>, 8> ray_table = make_ray_table();

} // namespace detail

// The squares beyond `from` in direction `d`, to the edge, ignoring everything on the board.
constexpr Bitboard ray(Direction d, Square from) {
    return detail::ray_table[static_cast<std::size_t>(d)][from.index];
}

// The same ray cut at the first occupied square, which it includes.
constexpr Bitboard ray_attacks(Direction d, Square from, Bitboard occupied) {
    const Bitboard full = ray(d, from);
    const Bitboard blockers = full & occupied;
    if (blockers.empty())
        return full;
    const Square first = detail::raises_index(d) ? blockers.lowest() : blockers.highest();
    return full ^ ray(d, first);
}

constexpr Bitboard rook_attacks(Square from, Bitboard occupied) {
    return ray_attacks(Direction::North, from, occupied) |
           ray_attacks(Direction::East, from, occupied) |
           ray_attacks(Direction::South, from, occupied) |
           ray_attacks(Direction::West, from, occupied);
}

constexpr Bitboard bishop_attacks(Square from, Bitboard occupied) {
    return ray_attacks(Direction::NorthEast, from, occupied) |
           ray_attacks(Direction::SouthEast, from, occupied) |
           ray_attacks(Direction::SouthWest, from, occupied) |
           ray_attacks(Direction::NorthWest, from, occupied);
}

constexpr Bitboard queen_attacks(Square from, Bitboard occupied) {
    return rook_attacks(from, occupied) | bishop_attacks(from, occupied);
}

// ── The pawn: the one piece whose pattern depends on the side ─────
//
// A pawn attacks the two squares diagonally ahead of it, and "ahead" is north
// for White and south for Black. So it is not carried to itself by every
// symmetry of the board: the left-right mirror keeps it, and the top-bottom
// flip turns White's pattern into Black's -- the colour symmetry the tests hold
// it to.

enum class Color : std::uint8_t { White, Black };

constexpr Color opposite(Color c) {
    return c == Color::White ? Color::Black : Color::White;
}

namespace detail {

constexpr std::array<std::array<Bitboard, square_count>, 2> make_pawn_table() {
    std::array<std::array<Bitboard, square_count>, 2> table{};
    for (int c = 0; c < 2; ++c)
        for (int s = 0; s < square_count; ++s) {
            const int rank = Square(s).rank() + (c == 0 ? 1 : -1);
            std::uint64_t bits = 0;
            for (int df : {-1, 1})
                if (on_board(Square(s).file() + df, rank))
                    bits |= Bitboard::of(Square::at(Square(s).file() + df, rank)).bits();
            table[c][s] = Bitboard(bits);
        }
    return table;
}

inline constexpr std::array<std::array<Bitboard, square_count>, 2> pawn_table =
    make_pawn_table();

} // namespace detail

// The squares a pawn of colour `c` on `from` attacks. Pushes are not attacks and are not here.
constexpr Bitboard pawn_attacks(Color c, Square from) {
    return detail::pawn_table[static_cast<std::size_t>(c)][from.index];
}

}  // namespace spatium::discrete
