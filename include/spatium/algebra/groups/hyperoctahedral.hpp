#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/concepts.hpp>
#  include <algorithm>
#  include <array>
#  include <cstddef>
#  include <cstdint>
#  include <numeric>
#endif

SPATIUM_EXPORT namespace spatium::inline algebra {

// ── The hyperoctahedral group B_N ─────────────────────────────────
//
// The symmetries of the N-cube: every signed permutation of N axes, N! 2^N
// elements. N = 2 is the dihedral group of the square (8 elements), N = 3 the
// full symmetry group of the cube (48).
//
// It acts on integer vectors (`SignedPermutationAction`) and on the cells of a
// cube of any side (`CubeCellAction`); both are `GroupAction`s in
// groups/action.hpp, so orbits, stabilisers and the action axioms come from
// there rather than from a loop written for this group.

template <std::size_t N> using IntVec = std::array<int, N>;

// The image of v has the entry sign[i] * v[source[i]] in place i.
template <std::size_t N> struct SignedPermutation {
    std::array<std::uint8_t, N> source{};
    std::array<std::int8_t, N> sign{};

    constexpr IntVec<N> apply(const IntVec<N>& v) const {
        IntVec<N> out{};
        for (std::size_t i = 0; i < N; ++i)
            out[i] = sign[i] * v[source[i]];
        return out;
    }

    constexpr bool operator==(const SignedPermutation&) const = default;
    constexpr auto operator<=>(const SignedPermutation&) const = default;
};

constexpr std::size_t hyperoctahedral_order(std::size_t n) {
    std::size_t order = std::size_t{1} << n;
    for (std::size_t i = 2; i <= n; ++i)
        order *= i;
    return order;
}

namespace hyperoctahedral_detail {

template <std::size_t N>
constexpr std::array<SignedPermutation<N>, hyperoctahedral_order(N)> make_hyperoctahedral() {
    std::array<SignedPermutation<N>, hyperoctahedral_order(N)> out{};
    std::array<std::uint8_t, N> perm{};
    std::iota(perm.begin(), perm.end(), std::uint8_t{0});
    std::size_t k = 0;
    do {
        for (std::size_t mask = 0; mask < (std::size_t{1} << N); ++mask) {
            SignedPermutation<N>& s = out[k++];
            s.source = perm;
            for (std::size_t i = 0; i < N; ++i)
                s.sign[i] = (mask >> i) & 1 ? std::int8_t{-1} : std::int8_t{1};
        }
    } while (std::next_permutation(perm.begin(), perm.end()));
    return out;
}

} // namespace hyperoctahedral_detail

// Every element once; the identity is the first.
template <std::size_t N>
inline constexpr std::array<SignedPermutation<N>, hyperoctahedral_order(N)>
    hyperoctahedral_elements = hyperoctahedral_detail::make_hyperoctahedral<N>();

// compose(x, y) applies y first.
template <std::size_t N> struct HyperoctahedralGroup {
    using ElementType = SignedPermutation<N>;

    constexpr ElementType identity() const { return hyperoctahedral_elements<N>[0]; }

    constexpr ElementType compose(const ElementType& x, const ElementType& y) const {
        ElementType out;
        for (std::size_t i = 0; i < N; ++i) {
            out.source[i] = y.source[x.source[i]];
            out.sign[i] = static_cast<std::int8_t>(x.sign[i] * y.sign[x.source[i]]);
        }
        return out;
    }

    constexpr ElementType inverse(const ElementType& x) const {
        ElementType out;
        for (std::size_t i = 0; i < N; ++i) {
            out.source[x.source[i]] = static_cast<std::uint8_t>(i);
            out.sign[x.source[i]] = x.sign[i];
        }
        return out;
    }

    constexpr const auto& elements() const { return hyperoctahedral_elements<N>; }
};

// On integer vectors: the step of a piece, a direction, a displacement.
template <std::size_t N> struct SignedPermutationAction {
    constexpr IntVec<N> act(const SignedPermutation<N>& g, const IntVec<N>& v) const {
        return g.apply(v);
    }
};

// On the cells 0..side-1 of each axis. The signed permutation acts on the doubled, centred
// coordinates 2c - (side - 1), which the cube's symmetries carry onto themselves.
template <std::size_t N> struct CubeCellAction {
    int side = 0;

    constexpr IntVec<N> act(const SignedPermutation<N>& g, const IntVec<N>& c) const {
        IntVec<N> centred{};
        for (std::size_t i = 0; i < N; ++i)
            centred[i] = 2 * c[i] - (side - 1);
        const IntVec<N> moved = g.apply(centred);
        IntVec<N> out{};
        for (std::size_t i = 0; i < N; ++i)
            out[i] = (moved[i] + (side - 1)) / 2;
        return out;
    }
};

static_assert(Group<HyperoctahedralGroup<4>>);

} // namespace spatium::algebra
