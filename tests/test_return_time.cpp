// How long a piece walking at random takes to come back.
//
// A piece on the empty 8x8 board moves to one of the cells it reaches, each with the same
// chance. The walk is a Markov chain on the 64 cells, and Kac's lemma gives the mean time to
// return to the cell v: 1 / pi(v) = 2|E| / deg(v), the edges counted in the part of the move
// graph that v is in. Four routes to that number, each built from a different part of the
// library:
//
//   the formula      the degrees, from lattice.hpp's attacks
//   a linear system  the mean hitting times, solved by algebra's solve_direct
//   the lumped chain the walk on the ten classes of cells the square's symmetries leave
//                    (algebra::orbits over CubeCellAction), solved the same way, times the
//                    size of v's class
//   the series       the sum of n f_n over the first-return probabilities f_n, walked out
//                    step by step until what is left cannot move the sum
//
// The lumped chain is the check on the symmetry: it exists only if every cell of a class
// moves into each other class with the same chance, which holds because each piece is an
// orbit of the group. Each route gets a case it must fail: the bishop, whose move graph is in
// two parts, breaks the formula counted over the whole board, and a piece with a direction
// cannot be lumped.

#include <spatium/algebra/groups/action.hpp>
#include <spatium/algebra/groups/hyperoctahedral.hpp>
#include <spatium/algebra/linear_solve.hpp>
#include <spatium/discrete/lattice.hpp>

#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

using namespace spatium;
using namespace spatium::discrete;

namespace {

constexpr std::size_t cells = 64;
constexpr std::size_t classes = 10; // the closed form m (m + 1) / 2 at m = 4, checked below

const Shape<2> board{{8, 8}};

using Chain = std::array<std::array<double, cells>, cells>;

struct Walk {
    Chain p{};                          // p[u][w], the chance to go from u to w
    std::array<std::size_t, cells> deg{};
};

Walk walk_of(const Piece<2>& piece) {
    Walk out;
    const auto empty = [](const Cell<2>&) { return false; };
    for (std::size_t u = 0; u < cells; ++u) {
        const auto reach = attacks(piece, board, board.cell_at(u), empty);
        out.deg[u] = reach.size();
        for (const Cell<2>& c : reach)
            out.p[u][board.index(c)] = 1.0 / static_cast<double>(reach.size());
    }
    return out;
}

// The cells the walk can reach from v.
std::array<bool, cells> part_of(const Walk& w, std::size_t v) {
    std::array<bool, cells> seen{};
    std::vector<std::size_t> todo{v};
    seen[v] = true;
    while (!todo.empty()) {
        const std::size_t u = todo.back();
        todo.pop_back();
        for (std::size_t x = 0; x < cells; ++x)
            if (w.p[u][x] > 0 && !seen[x]) {
                seen[x] = true;
                todo.push_back(x);
            }
    }
    return seen;
}

double kac_formula(const Walk& w, std::size_t v, const std::array<bool, cells>& counted) {
    std::size_t degree_sum = 0; // 2|E|
    for (std::size_t u = 0; u < cells; ++u)
        if (counted[u])
            degree_sum += w.deg[u];
    return static_cast<double>(degree_sum) / static_cast<double>(w.deg[v]);
}

// Mean return time to v from the mean hitting times h: h(v) = 0, h(u) = 1 + sum_x p(u, x) h(x).
// A cell the walk cannot reach gets h = 0 and takes no part.
template <std::size_t M>
double return_time_by_solve(const std::array<std::array<double, M>, M>& p, std::size_t v,
                            const std::array<bool, M>& part) {
    Matrix<double, M, M> a{};
    Vec<double, M> b{};
    for (std::size_t u = 0; u < M; ++u) {
        a(u, u) = 1;
        if (u == v || !part[u])
            continue;
        for (std::size_t x = 0; x < M; ++x)
            if (x != v)
                a(u, x) -= p[u][x];
        b[u] = 1;
    }
    const auto h = solve_direct(a, b);
    REQUIRE(h.has_value());
    double t = 1;
    for (std::size_t x = 0; x < M; ++x)
        t += p[v][x] * (*h)[x];
    return t;
}

struct Lumped {
    std::array<std::array<double, classes>, classes> q{};
    std::array<std::size_t, cells> class_of{};
    std::array<std::size_t, classes> size{};
};

// The walk on the classes of cells under the square's symmetries, or nothing if two cells of
// one class move into some class with different chances.
std::optional<Lumped> lumped(const Walk& w) {
    std::vector<Cell<2>> all;
    for (std::size_t u = 0; u < cells; ++u)
        all.push_back(board.cell_at(u));
    const auto classes_found =
        algebra::orbits(HyperoctahedralGroup<2>{}, CubeCellAction<2>{8}, std::span<const Cell<2>>(all));
    REQUIRE(classes_found.size() == classes);

    Lumped out;
    for (std::size_t k = 0; k < classes; ++k) {
        out.size[k] = classes_found[k].size();
        for (const Cell<2>& c : classes_found[k])
            out.class_of[board.index(c)] = k;
    }
    std::array<bool, classes> filled{};
    for (std::size_t u = 0; u < cells; ++u) {
        std::array<double, classes> row{};
        for (std::size_t x = 0; x < cells; ++x)
            row[out.class_of[x]] += w.p[u][x];
        const std::size_t k = out.class_of[u];
        if (!filled[k]) {
            out.q[k] = row;
            filled[k] = true;
            continue;
        }
        for (std::size_t j = 0; j < classes; ++j)
            if (std::abs(out.q[k][j] - row[j]) > 1e-12)
                return std::nullopt;
    }
    return out;
}

// sum n f_n, f_n the chance that the first return to v is at step n. `mass` is the chance of
// not having returned yet; the walk goes on until n * mass, a bound on the order of what is
// left of the sum when the mass falls geometrically, is below the tolerance.
double return_time_by_series(const Walk& w, std::size_t v) {
    std::array<double, cells> away = w.p[v]; // where the walk is, not yet back, after step 1
    away[v] = 0;
    double sum = w.p[v][v], mass = 1 - w.p[v][v];
    for (int n = 2; n < 100'000 && n * mass > 1e-11; ++n) {
        std::array<double, cells> next{};
        double back = 0;
        for (std::size_t u = 0; u < cells; ++u) {
            if (away[u] == 0)
                continue;
            back += away[u] * w.p[u][v];
            for (std::size_t x = 0; x < cells; ++x)
                if (x != v)
                    next[x] += away[u] * w.p[u][x];
        }
        sum += n * back;
        away = next;
        mass = 0; // summed afresh: subtracting each step's return would stall on rounding
        for (double a : away)
            mass += a;
    }
    return sum;
}

std::size_t at(int file, int rank) { return board.index({file, rank}); }

} // namespace

TEST_CASE("mean return times on the empty board: four routes agree", "[discrete][return-time]") {
    struct Case {
        const char* name;
        Piece<2> piece;
    };
    const std::array<Case, 4> pieces{{
        {"knight", knight<2>()}, {"king", king<2>()}, {"rook", rook<2>()}, {"queen", queen<2>()}
    }};
    for (const Case& c : pieces) {
        INFO(c.name);
        const Walk w = walk_of(c.piece);
        const auto lump = lumped(w);
        REQUIRE(lump.has_value());
        std::array<bool, classes> every_class{};
        every_class.fill(true);
        for (std::size_t v = 0; v < cells; ++v) {
            INFO("cell " << v);
            const auto part = part_of(w, v);
            const double formula = kac_formula(w, v, part);
            CHECK(return_time_by_solve(w.p, v, part) == Catch::Approx(formula).epsilon(1e-10));
            const std::size_t k = lump->class_of[v];
            const double by_class =
                return_time_by_solve(lump->q, k, every_class) * static_cast<double>(lump->size[k]);
            CHECK(by_class == Catch::Approx(formula).epsilon(1e-10));
        }
        // the series is slow, so at the corner and one inner cell only
        for (std::size_t v : {at(0, 0), at(3, 4)})
            CHECK(return_time_by_series(w, v) == Catch::Approx(kac_formula(w, v, part_of(w, v)))
                                                     .epsilon(1e-8));
    }

    // The numbers: 2|E| = 336 for the knight, 420 for the king, 64 * 14 for the rook.
    const auto corner = [](const Piece<2>& piece, std::size_t v) {
        const Walk w = walk_of(piece);
        return kac_formula(w, v, part_of(w, v));
    };
    CHECK(corner(knight<2>(), at(0, 0)) == 168);
    CHECK(corner(king<2>(), at(0, 0)) == 140);
    CHECK(corner(rook<2>(), at(3, 4)) == 64);
}

TEST_CASE("the routes fail where their assumptions do", "[discrete][return-time]") {
    SECTION("the bishop's graph is in two parts: counted over the whole board, twice too long") {
        const Walk w = walk_of(bishop<2>());
        std::array<bool, cells> whole{};
        whole.fill(true);
        for (std::size_t v : {at(0, 0), at(3, 4), at(7, 0)}) {
            const auto part = part_of(w, v);
            const double right = return_time_by_solve(w.p, v, part);
            CHECK(kac_formula(w, v, part) == Catch::Approx(right).epsilon(1e-10));
            CHECK(kac_formula(w, v, whole) == Catch::Approx(2 * right).epsilon(1e-10));
        }
    }

    SECTION("a piece with a direction cannot be lumped by the square's symmetries") {
        const Piece<2> lopsided = Piece<2>::of({
            {{1, 0}, Reach::Leap}, {{-1, 0}, Reach::Leap}, {{0, 1}, Reach::Leap},
            {{0, -1}, Reach::Leap}, {{1, 1}, Reach::Leap}
        });
        CHECK_FALSE(lumped(walk_of(lopsided)).has_value());
        CHECK(lumped(walk_of(king<2>())).has_value());
    }
}
