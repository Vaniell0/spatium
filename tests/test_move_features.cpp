// A move as a row of numbers.
//
// Each claim is held against a second route: `steps_taken` against the cells `attacks` walks
// to on an empty box, and the geometric row against every symmetry of the box (on two, three
// and four axes, a move and its image under a turn or a mirror give the same row). Each check
// gets a case it must fail: a piece with the wrong vector, and a row that keeps the axes in the
// order given.

#include <spatium/discrete/lattice.hpp>
#include <spatium/discrete/move_features.hpp>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <random>
#include <vector>

using namespace spatium::discrete;
using Catch::Matchers::WithinAbs;

namespace {

template <std::size_t N> std::vector<Cell<N>> all_cells(const Shape<N>& shape) {
    std::vector<Cell<N>> out;
    for (std::size_t i = 0; i < shape.size(); ++i)
        out.push_back(shape.cell_at(i));
    return out;
}

template <std::size_t N> Piece<N> piece_named(int which) {
    switch (which) {
    case 0:
        return rook<N>();
    case 1:
        return knight<N>();
    case 2:
        return queen<N>();
    default:
        return king<N>();
    }
}

} // namespace

TEST_CASE("move features: steps_taken agrees with where the piece reaches", "[move_features]") {
    const Shape<3> shape{{5, 5, 5}};
    const auto cells = all_cells(shape);
    for (int which = 0; which < 4; ++which) {
        const Piece<3> piece = piece_named<3>(which);
        for (const Cell<3>& from : {Cell<3>{0, 0, 0}, Cell<3>{2, 2, 2}, Cell<3>{4, 1, 3}}) {
            const auto reached = attacks(piece, shape, from, [](const Cell<3>&) { return false; });
            for (const Cell<3>& to : cells) {
                const bool listed = std::ranges::binary_search(reached, to);
                REQUIRE((steps_taken(piece, from, to) > 0) == listed);
            }
        }
    }
    SECTION("the count is the number of steps") {
        const Piece<3> r = rook<3>();
        REQUIRE(steps_taken(r, {1, 1, 1}, {4, 1, 1}) == 3);
        REQUIRE(steps_taken(r, {4, 1, 1}, {4, 1, 0}) == 1);
        REQUIRE(steps_taken(r, {1, 1, 1}, {2, 2, 1}) == 0);
        REQUIRE(steps_taken(knight<3>(), {2, 2, 2}, {3, 4, 2}) == 1);
        REQUIRE(steps_taken(knight<3>(), {2, 2, 2}, {4, 6, 2}) == 0);
    }
    SECTION("a piece with the wrong vector is caught") {
        // a rook is not a knight: the same pair of cells must give different answers
        REQUIRE(steps_taken(rook<3>(), {2, 2, 2}, {3, 4, 2}) == 0);
        REQUIRE(steps_taken(knight<3>(), {2, 2, 2}, {3, 4, 2}) == 1);
    }
}

namespace {

template <std::size_t N> void geometric_row_is_symmetric(int side) {
    const Shape<N> shape{[&] {
        std::array<int, N> e;
        e.fill(side);
        return e;
    }()};
    const auto cells = all_cells(shape);
    std::mt19937 rng(7 + N);
    for (int which = 0; which < 4; ++which) {
        const Piece<N> piece = piece_named<N>(which);
        for (int round = 0; round < 60; ++round) {
            const Cell<N> from = cells[rng() % cells.size()], to = cells[rng() % cells.size()];
            const auto row = geometric_features(piece, shape, from, to);
            for (const Symmetry<N>& s : spatium::algebra::hyperoctahedral_elements<N>) {
                const auto moved = geometric_features(
                    piece, shape, shape.transformed(s, from), shape.transformed(s, to)
                );
                for (std::size_t i = 0; i < row.size(); ++i)
                    REQUIRE_THAT(moved[i], WithinAbs(row[i], 1e-12));
            }
        }
    }
}

} // namespace

TEST_CASE(
    "move features: the geometric row is the same under every symmetry of the box", "[move_features]"
) {
    geometric_row_is_symmetric<2>(8);
    geometric_row_is_symmetric<3>(5);
    geometric_row_is_symmetric<4>(4);

    SECTION("a row that keeps the axes in the order given is not") {
        const Shape<3> shape{{5, 5, 5}};
        const Cell<3> from{0, 0, 0}, to{3, 1, 0};
        int differs = 0;
        for (const Symmetry<3>& s : spatium::algebra::hyperoctahedral_elements<3>) {
            const Cell<3> a = shape.transformed(s, from), b = shape.transformed(s, to);
            const std::array<int, 3> unsorted{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
            const std::array<int, 3> plain{to[0] - from[0], to[1] - from[1], to[2] - from[2]};
            differs += unsorted != plain;
        }
        REQUIRE(differs > 0);
    }
    SECTION("the entries are in range") {
        const Shape<2> shape{{8, 8}};
        const auto row = geometric_features(queen<2>(), shape, {0, 0}, {7, 7});
        REQUIRE_THAT(row[0], WithinAbs(1.0, 1e-12)); // seven along each axis of seven
        REQUIRE_THAT(row[2], WithinAbs(1.0, 1e-12)); // seven steps of seven
        REQUIRE_THAT(row[3], WithinAbs(1.0, 1e-12)); // a slide
        REQUIRE_THAT(row[4], WithinAbs(1.0, 1e-12)); // the far corner
        REQUIRE_THAT(row[6], WithinAbs(0.0, 1e-12)); // corner to corner: no nearer the centre
        const Shape<2> five{{5, 5}};
        const auto out = geometric_features(queen<2>(), five, {2, 2}, {4, 4});
        REQUIRE_THAT(out[6], WithinAbs(-1.0, 1e-12)); // from the centre to a corner
    }
}
