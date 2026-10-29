// The 8x8 board as a discrete space: squares, sets of squares in one word, the symmetries
// of the board, and the squares a piece reaches.
//
// Every answer is held against a second route to it: the fast word operations against the
// square-by-square definition, the table-driven attacks against a plain walk of the board,
// and the attack sets against the symmetry group (a piece that looks the same from every
// side must be carried to itself by each of the eight symmetries). Each verifier gets a
// case it must fail.
//
// The word is chosen over FiniteSet<int> for speed, not for looks. Measured once, clang 19
// -O2, intersecting two sets of 22 and 13 squares: 126 ns for FiniteSet (a sorted vector,
// one allocation per operation) against 2.6 ns for the word, about 50 times; the word side
// was a volatile-guarded AND plus popcount, so the ratio is rough.

#include <spatium/algebra/functions.hpp>
#include <spatium/discrete/attacks.hpp>
#include <spatium/discrete/bitboard.hpp>

#include <algorithm>
#include <bit>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <initializer_list>
#include <ranges>
#include <utility>

using namespace spatium::discrete;

namespace {

struct SplitMix64 {
    std::uint64_t state;
    std::uint64_t next() {
        std::uint64_t z = (state += 0x9e3779b97f4a7c15ull);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
        return z ^ (z >> 31);
    }
};

// Sparse and dense boards both: a quarter of the squares, or half.
Bitboard random_board(SplitMix64& rng, int i) {
    return Bitboard(i % 2 ? rng.next() & rng.next() : rng.next());
}

Square random_square(SplitMix64& rng) {
    return Square(static_cast<int>(rng.next() % 64));
}

// The reference routes: a plain walk of the board by coordinates, sharing no table with the code
// under test.
Bitboard walk(Square from, std::initializer_list<std::pair<int, int>> steps, Bitboard occupied) {
    Bitboard out;
    for (auto [df, dr] : steps)
        for (int f = from.file() + df, r = from.rank() + dr; on_board(f, r); f += df, r += dr) {
            out = out.with(Square::at(f, r));
            if (occupied.contains(Square::at(f, r)))
                break;
        }
    return out;
}

Bitboard jump(Square from, std::initializer_list<std::pair<int, int>> deltas) {
    Bitboard out;
    for (auto [df, dr] : deltas)
        if (on_board(from.file() + df, from.rank() + dr))
            out = out.with(Square::at(from.file() + df, from.rank() + dr));
    return out;
}

constexpr std::initializer_list<std::pair<int, int>> orthogonal{{0, 1}, {1, 0}, {0, -1}, {-1, 0}};
constexpr std::initializer_list<std::pair<int, int>> diagonal{{1, 1}, {1, -1}, {-1, -1}, {-1, 1}};
constexpr std::initializer_list<std::pair<int, int>> knight_jumps{{1, 2},  {2, 1},   {2, -1},
                                                                  {1, -2}, {-1, -2}, {-2, -1},
                                                                  {-2, 1}, {-1, 2}};
constexpr std::initializer_list<std::pair<int, int>> king_steps{{0, 1},  {1, 1},  {1, 0},
                                                                {1, -1}, {0, -1}, {-1, -1},
                                                                {-1, 0}, {-1, 1}};

// Failures of: a piece's attack set, moved by a symmetry, equals the attack set of the moved
// position.
template <class Attacks> int equivariance_failures(Attacks attacks, int samples) {
    SplitMix64 rng{0xC0FFEE};
    int failures = 0;
    for (const BoardSymmetry& g : all_symmetries)
        for (int i = 0; i < samples; ++i) {
            const Square s = random_square(rng);
            const Bitboard occupied = random_board(rng, i);
            failures +=
                attacks(g.apply(s), occupied.transformed(g)) != attacks(s, occupied).transformed(g);
        }
    return failures;
}

// Failures of: composing two symmetries and applying the result acts as applying them one after the
// other.
template <class Group> int homomorphism_failures(const Group& group) {
    int failures = 0;
    for (const auto& g : all_symmetries)
        for (const auto& h : all_symmetries)
            for (int s = 0; s < square_count; ++s)
                failures += group.compose(g, h).apply(Square(s)) != g.apply(h.apply(Square(s)));
    return failures;
}

// A composition with its arguments the wrong way round: not a homomorphism, and the check must say
// so.
struct ReversedGroup : BoardSymmetryGroup {
    constexpr ElementType compose(const ElementType& x, const ElementType& y) const {
        return BoardSymmetryGroup::compose(y, x);
    }
};

} // namespace

TEST_CASE("a square's file and rank give back its index", "[bitboard]") {
    for (int s = 0; s < square_count; ++s) {
        const Square sq(s);
        CHECK(Square::at(sq.file(), sq.rank()) == sq);
        CHECK(sq.index == sq.rank() * 8 + sq.file());
    }
    CHECK(Square::at(0, 0).index == 0);  // a1
    CHECK(Square::at(7, 7).index == 63); // h8
    CHECK_FALSE(on_board(-1, 0));
    CHECK_FALSE(on_board(0, 8));
}

TEST_CASE("a bitboard is a set of squares with iteration in index order", "[bitboard]") {
    static_assert(std::ranges::input_range<Bitboard>);

    SplitMix64 rng{1};
    for (int i = 0; i < 200; ++i) {
        const Bitboard board = random_board(rng, i);
        int n = 0, previous = -1;
        for (Square s : board) {
            CHECK(board.contains(s));
            CHECK(s.index > previous);
            previous = s.index;
            ++n;
        }
        CHECK(n == board.count());
        CHECK(board.empty() == (n == 0));
        if (!board.empty()) {
            CHECK(board.lowest().index == std::countr_zero(board.bits()));
            CHECK(board.highest().index == 63 - std::countl_zero(board.bits()));
        }
    }

    const Square e4 = Square::at(4, 3);
    const Bitboard one = Bitboard::of(e4);
    CHECK(one.count() == 1);
    CHECK(one.with(e4) == one);
    CHECK(one.without(e4).empty());
    CHECK((one | Bitboard::of(Square(0))).count() == 2);
    CHECK((one & ~one).empty());
    CHECK((one ^ one).empty());
}

TEST_CASE(
    "the word forms of the board symmetries agree with the square-by-square definition",
    "[bitboard][symmetry]"
) {
    SplitMix64 rng{2};
    for (int i = 0; i < 500; ++i) {
        const Bitboard board = random_board(rng, i);
        CHECK(board.flipped_ranks() == board.transformed(flip_ranks));
        CHECK(board.mirrored_files() == board.transformed(mirror_files));
        CHECK(board.transposed() == board.transformed(transpose));
    }
}

TEST_CASE(
    "a word form that is not a symmetry is reported against the definition", "[bitboard][symmetry]"
) {
    // Shifting every rank up by one (a rotation of the word) is not the top-bottom flip.
    SplitMix64 rng{3};
    int disagreements = 0;
    for (int i = 0; i < 100; ++i) {
        const Bitboard board = random_board(rng, i);
        disagreements += Bitboard(std::rotl(board.bits(), 8)) != board.transformed(flip_ranks);
    }
    CHECK(disagreements > 90);
}

TEST_CASE("the eight symmetries of the board form the group D4", "[bitboard][symmetry]") {
    const BoardSymmetryGroup group;

    // Eight distinct elements, closed under composition, with the identity and inverses.
    for (std::size_t i = 0; i < all_symmetries.size(); ++i)
        for (std::size_t j = i + 1; j < all_symmetries.size(); ++j)
            CHECK(all_symmetries[i] != all_symmetries[j]);

    for (const auto& g : all_symmetries) {
        CHECK(group.compose(g, group.identity()) == g);
        CHECK(group.compose(group.identity(), g) == g);
        CHECK(group.compose(g, group.inverse(g)) == group.identity());
        for (const auto& h : all_symmetries) {
            const auto gh = group.compose(g, h);
            CHECK(std::ranges::find(all_symmetries, gh) != all_symmetries.end());
            for (const auto& k : all_symmetries)
                CHECK(
                    group.compose(group.compose(g, h), k) == group.compose(g, group.compose(h, k))
                );
        }
    }

    // It is not commutative: a quarter turn and a mirror do not commute.
    CHECK(group.compose(rotate_90, mirror_files) != group.compose(mirror_files, rotate_90));

    // Orders, through the library's own generic group algorithm.
    CHECK(spatium::algebra::power(group, rotate_90, 4) == identity_symmetry);
    CHECK(spatium::algebra::power(group, rotate_90, 2) == rotate_180);
    CHECK(spatium::algebra::power(group, mirror_files, 2) == identity_symmetry);
    CHECK(spatium::algebra::power(group, transpose, 2) == identity_symmetry);
}

TEST_CASE("the symmetries act on squares as the group says", "[bitboard][symmetry]") {
    const BoardSymmetryGroup group;
    CHECK(homomorphism_failures(group) == 0);

    // Each is a permutation of the sixty-four squares.
    for (const auto& g : all_symmetries) {
        Bitboard image;
        for (int s = 0; s < square_count; ++s)
            image = image.with(g.apply(Square(s)));
        CHECK(image.count() == square_count);
    }

    // Corners: a quarter turn counter-clockwise carries a1 to h1.
    const Square a1 = Square::at(0, 0), h1 = Square::at(7, 0), a8 = Square::at(0, 7),
                 h8 = Square::at(7, 7);
    CHECK(rotate_90.apply(a1) == h1);
    CHECK(rotate_180.apply(a1) == h8);
    CHECK(mirror_files.apply(a1) == h1);
    CHECK(flip_ranks.apply(a1) == a8);
    CHECK(transpose.apply(h1) == a8);
    CHECK(anti_transpose.apply(a1) == h8);
}

TEST_CASE(
    "a composition in the wrong order is caught by the action check", "[bitboard][symmetry]"
) {
    CHECK(homomorphism_failures(ReversedGroup{}) > 0);
}

TEST_CASE(
    "the table attacks of the knight and the king match a plain walk of the board",
    "[bitboard][attacks]"
) {
    for (int s = 0; s < square_count; ++s) {
        const Square from(s);
        CHECK(knight_attacks(from) == jump(from, knight_jumps));
        CHECK(king_attacks(from) == jump(from, king_steps));
    }
    // The knight is the orbit of (1, 2); the camel, the orbit of (1, 3), is a leaper of the same
    // kind.
    for (int s = 0; s < square_count; ++s) {
        const Square from(s);
        CHECK(leaper_attacks<1, 2>(from) == knight_attacks(from));
        CHECK(
            leaper_attacks<1, 3>(from) ==
            jump(from, {{1, 3}, {3, 1}, {3, -1}, {1, -3}, {-1, -3}, {-3, -1}, {-3, 1}, {-1, 3}})
        );
    }
}

TEST_CASE("the ray attacks of the sliders match a plain walk of the board", "[bitboard][attacks]") {
    SplitMix64 rng{4};
    for (int s = 0; s < square_count; ++s)
        for (int i = 0; i < 300; ++i) {
            const Square from(s);
            const Bitboard occupied = random_board(rng, i);
            CHECK(rook_attacks(from, occupied) == walk(from, orthogonal, occupied));
            CHECK(bishop_attacks(from, occupied) == walk(from, diagonal, occupied));
            CHECK(
                queen_attacks(from, occupied) ==
                (walk(from, orthogonal, occupied) | walk(from, diagonal, occupied))
            );
        }
}

TEST_CASE("attack sets have the counts every board shows", "[bitboard][attacks]") {
    const Bitboard empty;
    for (int s = 0; s < square_count; ++s)
        CHECK(rook_attacks(Square(s), empty).count() == 14);

    CHECK(bishop_attacks(Square::at(4, 3), empty).count() == 13); // e4
    CHECK(queen_attacks(Square::at(3, 3), empty).count() == 27);  // d4
    CHECK(knight_attacks(Square::at(0, 0)).count() == 2);         // a1
    CHECK(knight_attacks(Square::at(3, 3)).count() == 8);         // d4
    CHECK(king_attacks(Square::at(0, 0)).count() == 3);           // a1
    CHECK(king_attacks(Square::at(4, 3)).count() == 8);           // e4

    // A blocker is attacked and hides what is behind it.
    const Square a1 = Square::at(0, 0), a4 = Square::at(0, 3), d1 = Square::at(3, 0);
    const Bitboard blocked = Bitboard::of(a4).with(d1);
    const Bitboard reached = rook_attacks(a1, blocked);
    CHECK(reached.count() == 6); // a2 a3 a4 and b1 c1 d1
    CHECK(reached.contains(a4));
    CHECK(reached.contains(d1));
    CHECK_FALSE(reached.contains(Square::at(0, 4))); // a5, behind a4
    CHECK_FALSE(reached.contains(Square::at(4, 0))); // e1, behind d1
}

TEST_CASE(
    "pieces that look the same from every side are carried to themselves by every symmetry",
    "[bitboard][attacks][symmetry]"
) {
    CHECK(equivariance_failures([](Square s, Bitboard o) { return rook_attacks(s, o); }, 200) == 0);
    CHECK(
        equivariance_failures([](Square s, Bitboard o) { return bishop_attacks(s, o); }, 200) == 0
    );
    CHECK(
        equivariance_failures([](Square s, Bitboard o) { return queen_attacks(s, o); }, 200) == 0
    );
    CHECK(equivariance_failures([](Square s, Bitboard) { return knight_attacks(s); }, 200) == 0);
    CHECK(equivariance_failures([](Square s, Bitboard) { return king_attacks(s); }, 200) == 0);
    CHECK(
        equivariance_failures([](Square s, Bitboard) { return leaper_attacks<1, 3>(s); }, 200) == 0
    );
}

TEST_CASE(
    "a piece with a preferred direction is reported by the symmetry check",
    "[bitboard][attacks][symmetry]"
) {
    // A ray that only goes north looks different from each side, so some symmetry must expose it.
    CHECK(
        equivariance_failures(
            [](Square s, Bitboard o) { return ray_attacks(Direction::North, s, o); }, 200
        ) > 0
    );
    // And a leaper with an uneven pattern: one knight jump only.
    CHECK(equivariance_failures([](Square s, Bitboard) { return jump(s, {{1, 2}}); }, 200) > 0);
}

TEST_CASE(
    "a pawn attacks the two squares diagonally ahead, the way its side faces", "[bitboard][attacks]"
) {
    const Square e4 = Square::at(4, 3), a2 = Square::at(0, 1), h7 = Square::at(7, 6);
    CHECK(pawn_attacks(Color::White, e4) == Bitboard::of(Square::at(3, 4)).with(Square::at(5, 4)));
    CHECK(pawn_attacks(Color::Black, e4) == Bitboard::of(Square::at(3, 2)).with(Square::at(5, 2)));
    // On the edge it has one diagonal; on the last rank none.
    CHECK(pawn_attacks(Color::White, a2) == Bitboard::of(Square::at(1, 2)));
    CHECK(pawn_attacks(Color::Black, h7) == Bitboard::of(Square::at(6, 5)));
    for (int file = 0; file < 8; ++file) {
        CHECK(pawn_attacks(Color::White, Square::at(file, 7)).empty());
        CHECK(pawn_attacks(Color::Black, Square::at(file, 0)).empty());
    }
    CHECK(opposite(Color::White) == Color::Black);
}

TEST_CASE(
    "the pawn is carried to itself by the mirror and to the other colour by the flip",
    "[bitboard][attacks][symmetry]"
) {
    for (int s = 0; s < square_count; ++s) {
        const Square from(s);
        for (Color c : {Color::White, Color::Black}) {
            // Left-right: the same pattern, mirrored.
            CHECK(
                pawn_attacks(c, mirror_files.apply(from)) ==
                pawn_attacks(c, from).transformed(mirror_files)
            );
            // Top-bottom: White's pattern becomes Black's.
            CHECK(
                pawn_attacks(opposite(c), flip_ranks.apply(from)) ==
                pawn_attacks(c, from).transformed(flip_ranks)
            );
        }
    }
}

TEST_CASE(
    "a pawn is not carried to itself by every symmetry, and the check says so",
    "[bitboard][attacks][symmetry]"
) {
    // A quarter turn sends the diagonals ahead to the diagonals beside.
    int failures = 0;
    for (int s = 0; s < square_count; ++s)
        failures += pawn_attacks(Color::White, rotate_90.apply(Square(s))) !=
                    pawn_attacks(Color::White, Square(s)).transformed(rotate_90);
    CHECK(failures > 0);
}
