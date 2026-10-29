// The hyperoctahedral group and group actions.
//
// Each claim has a second route to it: the order against N! 2^N, acting by a composite against
// acting twice, the number of orbits found by walking them against Burnside's count (the mean
// number of points each element fixes), and that count against the closed form for the cells
// of a square. Each check gets a case it must fail: a composition that runs the wrong way, and
// an action that moves points under the identity.

#include <spatium/algebra/groups/action.hpp>
#include <spatium/algebra/groups/hyperoctahedral.hpp>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <vector>

using namespace spatium;

namespace {

template <std::size_t N> std::vector<IntVec<N>> cells(int side) {
    std::vector<IntVec<N>> out;
    std::size_t total = 1;
    for (std::size_t i = 0; i < N; ++i)
        total *= static_cast<std::size_t>(side);
    for (std::size_t at = 0; at < total; ++at) {
        IntVec<N> c{};
        std::size_t rest = at;
        for (std::size_t i = 0; i < N; ++i) {
            c[i] = static_cast<int>(rest % static_cast<std::size_t>(side));
            rest /= static_cast<std::size_t>(side);
        }
        out.push_back(c);
    }
    return out;
}

// Burnside: the number of orbits is the mean over the group of the number of points each
// element leaves in place. It never looks at an orbit.
template <class G, class A, class X>
std::size_t burnside_count(const G& group, const A& action, const std::vector<X>& points) {
    std::size_t fixed = 0;
    for (const auto& g : group.elements())
        for (const X& x : points)
            fixed += action.act(g, x) == x;
    return fixed / group.elements().size();
}

template <std::size_t N> struct Backwards : HyperoctahedralGroup<N> {
    constexpr SignedPermutation<N> compose(
        const SignedPermutation<N>& x, const SignedPermutation<N>& y
    ) const {
        return HyperoctahedralGroup<N>::compose(y, x);
    }
};

// Centres the cells without doubling them: the identity then halves every coordinate, so it
// moves every cell but the origin.
template <std::size_t N> struct UndoubledCellAction {
    int side = 0;
    constexpr IntVec<N> act(const SignedPermutation<N>& g, const IntVec<N>& c) const {
        IntVec<N> centred{};
        for (std::size_t i = 0; i < N; ++i)
            centred[i] = c[i] - (side - 1);
        const IntVec<N> moved = g.apply(centred);
        IntVec<N> out{};
        for (std::size_t i = 0; i < N; ++i)
            out[i] = (moved[i] + (side - 1)) / 2;
        return out;
    }
};

} // namespace

static_assert(FiniteGroup<HyperoctahedralGroup<3>>);
static_assert(GroupAction<SignedPermutationAction<3>, HyperoctahedralGroup<3>, IntVec<3>>);
static_assert(GroupAction<CubeCellAction<3>, HyperoctahedralGroup<3>, IntVec<3>>);

TEST_CASE("the hyperoctahedral group has N! 2^N elements and is a group", "[algebra][group]") {
    STATIC_CHECK(hyperoctahedral_order(1) == 2);
    STATIC_CHECK(hyperoctahedral_order(2) == 8);
    STATIC_CHECK(hyperoctahedral_order(3) == 48);
    STATIC_CHECK(hyperoctahedral_order(4) == 384);
    STATIC_CHECK(hyperoctahedral_order(5) == 3840);

    auto all = hyperoctahedral_elements<4>;
    std::ranges::sort(all);
    CHECK(std::ranges::adjacent_find(all) == all.end());

    const HyperoctahedralGroup<3> group;
    const auto& g = group.elements();
    const auto member = [&](const SignedPermutation<3>& s) {
        return std::ranges::find(g, s) != g.end();
    };
    for (const auto& x : g) {
        CHECK(group.compose(x, group.identity()) == x);
        CHECK(group.compose(group.identity(), x) == x);
        CHECK(group.compose(x, group.inverse(x)) == group.identity());
        CHECK(group.compose(group.inverse(x), x) == group.identity());
        for (const auto& y : g) {
            REQUIRE(member(group.compose(x, y)));
            for (const auto& z : g)
                REQUIRE(
                    group.compose(group.compose(x, y), z) == group.compose(x, group.compose(y, z))
                );
        }
    }
    const IntVec<4> v{3, -1, 4, 2};
    CHECK(HyperoctahedralGroup<4>{}.identity().apply(v) == v);
}

TEST_CASE("the signed permutations act on vectors and on the cells of a cube", "[algebra][group]") {
    const std::vector<IntVec<3>> vectors{{3, -1, 4}, {1, 2, 0}, {0, 0, 0}, {1, 1, 1}, {-2, 5, 5}};
    CHECK(verify_action(HyperoctahedralGroup<3>{}, SignedPermutationAction<3>{},
                        std::span<const IntVec<3>>(vectors)));

    for (int side : {1, 2, 3, 4, 8}) {
        const auto square = cells<2>(side);
        CHECK(verify_action(HyperoctahedralGroup<2>{}, CubeCellAction<2>{side},
                            std::span<const IntVec<2>>(square)));
    }
    const auto cube = cells<3>(4);
    CHECK(verify_action(HyperoctahedralGroup<3>{}, CubeCellAction<3>{4},
                        std::span<const IntVec<3>>(cube)));
    const auto tesseract = cells<4>(3);
    CHECK(verify_action(HyperoctahedralGroup<4>{}, CubeCellAction<4>{3},
                        std::span<const IntVec<4>>(tesseract)));
}

TEST_CASE("verify_action rejects a backwards composition and a moving identity", "[algebra][group]") {
    const std::vector<IntVec<3>> vectors{{3, -1, 4}, {1, 2, 0}};
    const auto backwards = verify_action(Backwards<3>{}, SignedPermutationAction<3>{},
                                         std::span<const IntVec<3>>(vectors));
    CHECK_FALSE(backwards);
    CHECK(backwards.failure == "act(compose(g, h), x) != act(g, act(h, x))");

    // (0, 0) is left in place by its identity, (3, 5) is not
    const std::vector<IntVec<2>> moved{{3, 5}, {0, 0}};
    const auto undoubled = verify_action(HyperoctahedralGroup<2>{}, UndoubledCellAction<2>{8},
                                         std::span<const IntVec<2>>(moved));
    CHECK_FALSE(undoubled);
    CHECK(undoubled.failure == "act(identity, x) != x");
}

TEST_CASE("orbits: walked, counted by Burnside, and by the closed form", "[algebra][group]") {
    SECTION("a step's orbit is every signed arrangement of its entries") {
        const HyperoctahedralGroup<4> group;
        const SignedPermutationAction<4> on_vectors;
        // entries with k nonzero, all equal: C(4, k) 2^k; the knight (1, 2, 0, 0): 4 * 3 * 4
        CHECK(orbit(group, on_vectors, IntVec<4>{1, 0, 0, 0}).size() == 8);
        CHECK(orbit(group, on_vectors, IntVec<4>{1, 1, 0, 0}).size() == 24);
        CHECK(orbit(group, on_vectors, IntVec<4>{1, 1, 1, 0}).size() == 32);
        CHECK(orbit(group, on_vectors, IntVec<4>{1, 1, 1, 1}).size() == 16);
        CHECK(orbit(group, on_vectors, IntVec<4>{1, 2, 0, 0}).size() == 48);
        CHECK(stabilizer(group, on_vectors, IntVec<4>{0, 0, 0, 0}).size() == 384);
    }

    SECTION("the cells of a square side n fall into floor((n+1)/2) (floor((n+1)/2) + 1) / 2 "
            "classes") {
        const HyperoctahedralGroup<2> group;
        for (int side = 1; side <= 9; ++side) {
            const auto square = cells<2>(side);
            const CubeCellAction<2> on_cells{side};
            const auto walked = orbits(group, on_cells, std::span<const IntVec<2>>(square));
            const std::size_t m = static_cast<std::size_t>((side + 1) / 2);
            CHECK(walked.size() == m * (m + 1) / 2);
            CHECK(burnside_count(group, on_cells, square) == walked.size());
        }
        // the board: ten classes, the corner's orbit the four corners
        const CubeCellAction<2> board{8};
        CHECK(orbit(group, board, IntVec<2>{0, 0}) ==
              std::vector<IntVec<2>>{{0, 0}, {0, 7}, {7, 0}, {7, 7}});
    }

    SECTION("on the cube and the tesseract Burnside agrees with the walk") {
        const auto cube = cells<3>(4);
        const auto tesseract = cells<4>(3);
        const CubeCellAction<3> on_cube{4};
        const CubeCellAction<4> on_tesseract{3};
        CHECK(orbits(HyperoctahedralGroup<3>{}, on_cube, std::span<const IntVec<3>>(cube)).size() ==
              burnside_count(HyperoctahedralGroup<3>{}, on_cube, cube));
        CHECK(orbits(HyperoctahedralGroup<4>{}, on_tesseract, std::span<const IntVec<4>>(tesseract))
                  .size() == burnside_count(HyperoctahedralGroup<4>{}, on_tesseract, tesseract));
    }
}
