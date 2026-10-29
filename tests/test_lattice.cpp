// A board of any number of axes, and pieces that are vectors.
//
// Every claim is held against a second route to it: the size of a piece's orbit against the
// count by choosing axes and signs, the reach of a piece on an empty board against a count
// from the shape of the displacement, the walking definition against the word-parallel tables
// of attacks.hpp on 8x8, and the reach of every piece against the symmetries of the box (a
// piece that looks the same from every side is carried to itself by each of them). Each check
// gets a case it must fail: a knight with the wrong vector, a slider that ignores blockers,
// and a pawn, whose step has a direction.

#include <spatium/discrete/attacks.hpp>
#include <spatium/discrete/bitboard.hpp>
#include <spatium/discrete/lattice.hpp>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <random>
#include <vector>

using namespace spatium::discrete;

namespace {

std::size_t choose(std::size_t n, std::size_t k) {
    std::size_t c = 1;
    for (std::size_t i = 1; i <= k; ++i)
        c = c * (n - k + i) / i;
    return c;
}

std::size_t power(std::size_t base, std::size_t exp) {
    std::size_t p = 1;
    for (std::size_t i = 0; i < exp; ++i)
        p *= base;
    return p;
}

template <std::size_t N> auto nobody() {
    return [](const Cell<N>&) { return false; };
}

template <std::size_t N> Shape<N> cube(int side) {
    Shape<N> s;
    s.extent.fill(side);
    return s;
}

template <std::size_t N>
bool carried_to_itself(
    const Piece<N>& piece, const Shape<N>& shape, std::mt19937_64& rng, int occupancies
) {
    std::bernoulli_distribution coin(0.3);
    for (int round = 0; round < occupancies; ++round) {
        std::vector<bool> occ(shape.size());
        for (std::size_t i = 0; i < occ.size(); ++i)
            occ[i] = coin(rng);
        for (const Symmetry<N>& g : spatium::algebra::hyperoctahedral_elements<N>) {
            std::vector<bool> moved(occ.size());
            for (std::size_t i = 0; i < occ.size(); ++i)
                moved[shape.index(shape.transformed(g, shape.cell_at(i)))] = occ[i];
            for (std::size_t i = 0; i < shape.size(); ++i) {
                const Cell<N> from = shape.cell_at(i);
                auto reach = attacks(piece, shape, from, [&](const Cell<N>& c) {
                    return occ[shape.index(c)];
                });
                for (Cell<N>& c : reach)
                    c = shape.transformed(g, c);
                std::ranges::sort(reach);
                const auto image =
                    attacks(piece, shape, shape.transformed(g, from), [&](const Cell<N>& c) {
                        return moved[shape.index(c)];
                    });
                if (reach != image)
                    return false;
            }
        }
    }
    return true;
}

// The pawn of White: one step north, no symmetry to carry it to itself.
Piece<2> pawn_pushes() {
    return Piece<2>::of({Step<2>{{0, 1}, Reach::Leap}});
}

} // namespace

TEST_CASE("the orbit of a step has the size that choosing axes and signs gives", "[lattice]") {
    auto check = []<std::size_t N>() {
        // rook: one axis and a sign; the piece along k axes: which axes, which signs
        CHECK(rook<N>().size() == 2 * N);
        CHECK(diagonal_mover<N>(N, Reach::Slide).size() == power(2, N));
        for (std::size_t k = 1; k <= N; ++k)
            CHECK(diagonal_mover<N>(k, Reach::Slide).size() == choose(N, k) * power(2, k));
        // king and queen: every other cell of the 3^N neighbourhood
        CHECK(king<N>().size() == power(3, N) - 1);
        CHECK(queen<N>().size() == power(3, N) - 1);
        // knight: which axis jumps two, which jumps one, and two signs
        CHECK(knight<N>().size() == 4 * N * (N - 1));
    };
    check.template operator()<2>();
    check.template operator()<3>();
    check.template operator()<4>();
    check.template operator()<5>();

    // The numbers the pieces of four axes are known by.
    CHECK(rook<4>().size() == 8);
    CHECK(bishop<4>().size() == 24);
    CHECK(unicorn<4>().size() == 32);
    CHECK(dragon<4>().size() == 16);
    CHECK(knight<4>().size() == 48);
    CHECK(queen<4>().size() == 80);
    // and on the plane, the pieces of the board
    CHECK(rook<2>().size() == 4);
    CHECK(bishop<2>().size() == 4);
    CHECK(knight<2>().size() == 8);
    CHECK(king<2>().size() == 8);
}

TEST_CASE("a slide wins over a leap along the same vector", "[lattice]") {
    const Piece<2> both = Piece<2>::orbit({1, 0}, Reach::Leap) + rook<2>();
    CHECK(both.size() == 4);
    for (const Step<2>& s : both.steps())
        CHECK(s.reach == Reach::Slide);
    // a king and a rook share the four orthogonal vectors: the king's step is the slide's first
    const Piece<2> kr = king<2>() + rook<2>();
    CHECK(kr.size() == 8);
    CHECK(std::ranges::count(kr.steps(), Reach::Slide, &Step<2>::reach) == 4);
}

TEST_CASE(
    "the symmetries of the square are the same group as the board's",
    "[lattice][symmetry]"
) {
    SECTION("on two axes it is the dihedral group of the board") {
        const Shape<2> board = cube<2>(board_size);
        for (const Symmetry<2>& s : spatium::algebra::hyperoctahedral_elements<2>) {
            // the matrix of s: row i has the entry sign[i] in column source[i]
            const int a = s.source[0] == 0 ? s.sign[0] : 0, b = s.source[0] == 1 ? s.sign[0] : 0;
            const int c = s.source[1] == 0 ? s.sign[1] : 0, d = s.source[1] == 1 ? s.sign[1] : 0;
            const auto match = std::ranges::find_if(all_symmetries, [&](const BoardSymmetry& t) {
                return t.a == a && t.b == b && t.c == c && t.d == d;
            });
            REQUIRE(match != all_symmetries.end());
            for (int sq = 0; sq < square_count; ++sq) {
                const Square from(sq);
                const Cell<2> image = board.transformed(s, {from.file(), from.rank()});
                const Square moved = match->apply(from);
                REQUIRE(image == Cell<2>{moved.file(), moved.rank()});
            }
        }
    }
}

TEST_CASE("the reach on an empty box is counted from the shape of the displacement", "[lattice]") {
    // From a corner of a box with side L, the queen of N axes reaches the displacements whose
    // non-zero entries are all the same positive number: for each count k of axes, C(N,k)
    // choices of axes and L-1 sizes.
    SECTION("queen from a corner of eight cells to the side on four axes") {
        const Shape<4> box = cube<4>(8);
        std::size_t expected = 0;
        for (std::size_t k = 1; k <= 4; ++k)
            expected += choose(4, k) * 7;
        CHECK(expected == 105);
        CHECK(attacks(queen<4>(), box, {0, 0, 0, 0}, nobody<4>()).size() == expected);
        CHECK(attacks(king<4>(), box, {0, 0, 0, 0}, nobody<4>()).size() == 15);
    }
    // From the centre of a cube of side 5, every direction runs two cells, and a knight's
    // jump stays inside.
    SECTION("queen from the centre of a cube of five cells") {
        CHECK(attacks(queen<3>(), cube<3>(5), {2, 2, 2}, nobody<3>()).size() == 2 * 26);
        CHECK(attacks(queen<4>(), cube<4>(5), {2, 2, 2, 2}, nobody<4>()).size() == 2 * 80);
        CHECK(attacks(knight<4>(), cube<4>(5), {2, 2, 2, 2}, nobody<4>()).size() == 48);
    }
    SECTION("the visitor sees each cell the list holds, without keeping them") {
        const Shape<3> box = cube<3>(4);
        std::size_t seen = 0;
        for_each_attack(queen<3>(), box, {1, 1, 1}, nobody<3>(), [&](const Cell<3>&) { ++seen; });
        CHECK(seen == attacks(queen<3>(), box, {1, 1, 1}, nobody<3>()).size());
    }
    SECTION("a blocker is reached and nothing past it") {
        const Shape<2> board = cube<2>(8);
        auto at_d4 = [](const Cell<2>& c) { return c == Cell<2>{3, 3}; };
        const auto reach = attacks(rook<2>(), board, {3, 0}, at_d4);
        CHECK(std::ranges::count(reach, Cell<2>{3, 3}) == 1);
        CHECK(std::ranges::count(reach, Cell<2>{3, 4}) == 0);
        CHECK(reach.size() == 3 + 7);
    }
}

TEST_CASE("the walking definition agrees with the tables of the 8x8 board", "[lattice][bitboard]") {
    const Shape<2> board = cube<2>(board_size);
    std::mt19937_64 rng(2026);

    struct Row {
        Piece<2> piece;
        Bitboard (*fast)(Square, Bitboard);
    };
    const Row rows[] = {
        {knight<2>(), [](Square s, Bitboard) { return knight_attacks(s); }},
        {king<2>(), [](Square s, Bitboard) { return king_attacks(s); }},
        {rook<2>(), [](Square s, Bitboard o) { return rook_attacks(s, o); }},
        {bishop<2>(), [](Square s, Bitboard o) { return bishop_attacks(s, o); }},
        {queen<2>(), [](Square s, Bitboard o) { return queen_attacks(s, o); }},
    };
    auto as_set = [&](const std::vector<Cell<2>>& cells) {
        Bitboard b;
        for (const Cell<2>& c : cells)
            b = b.with(Square::at(c[0], c[1]));
        return b;
    };

    for (const Row& row : rows)
        for (int round = 0; round < 200; ++round) {
            const Bitboard occ(rng());
            auto occupied = [&](const Cell<2>& c) { return occ.contains(Square::at(c[0], c[1])); };
            for (int sq = 0; sq < square_count; ++sq) {
                const Square from(sq);
                REQUIRE(
                    as_set(attacks(row.piece, board, {from.file(), from.rank()}, occupied)) ==
                    row.fast(from, occ)
                );
            }
        }

    SECTION("a wrong knight vector and a slider that ignores blockers both disagree") {
        const Piece<2> camel = Piece<2>::orbit({1, 3}, Reach::Leap);
        const Square e4 = Square::at(4, 3);
        CHECK(as_set(attacks(camel, board, {4, 3}, nobody<2>())) != knight_attacks(e4));

        const Bitboard occ(0x0000'0010'0000'1000ull);
        auto occupied = [&](const Cell<2>& c) { return occ.contains(Square::at(c[0], c[1])); };
        const auto blind = as_set(attacks(rook<2>(), board, {4, 3}, nobody<2>()));
        CHECK(blind != rook_attacks(e4, occ));
        CHECK(as_set(attacks(rook<2>(), board, {4, 3}, occupied)) == rook_attacks(e4, occ));
    }
}

TEST_CASE(
    "every piece that looks the same from every side is carried to itself by the cube's symmetries",
    "[lattice][symmetry]"
) {
    std::mt19937_64 rng(11);

    SECTION("three axes, sides of three and four") {
        for (int side : {3, 4}) {
            const Shape<3> box = cube<3>(side);
            CHECK(carried_to_itself<3>(rook<3>(), box, rng, 3));
            CHECK(carried_to_itself<3>(bishop<3>(), box, rng, 3));
            CHECK(carried_to_itself<3>(unicorn<3>(), box, rng, 3));
            CHECK(carried_to_itself<3>(queen<3>(), box, rng, 3));
            CHECK(carried_to_itself<3>(king<3>(), box, rng, 3));
            CHECK(carried_to_itself<3>(knight<3>(), box, rng, 3));
        }
    }
    SECTION("four axes, side three") {
        const Shape<4> box = cube<4>(3);
        CHECK(carried_to_itself<4>(rook<4>(), box, rng, 1));
        CHECK(carried_to_itself<4>(dragon<4>(), box, rng, 1));
        CHECK(carried_to_itself<4>(knight<4>(), box, rng, 1));
    }
    SECTION("the pawn is not: a step with a direction is not an orbit") {
        CHECK_FALSE(carried_to_itself<2>(pawn_pushes(), cube<2>(board_size), rng, 1));
        // and the orbit of the same vector is: the pawn's four-fold image is the king's orthogonal
        // step
        CHECK(
            carried_to_itself<2>(Piece<2>::orbit({0, 1}, Reach::Leap), cube<2>(board_size), rng, 1)
        );
    }
}
