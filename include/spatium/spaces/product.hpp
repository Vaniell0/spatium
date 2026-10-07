#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/access.hpp>
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <cmath>
#  include <concepts>
#  include <cstddef>
#  include <tuple>
#  include <type_traits>
#endif

SPATIUM_EXPORT namespace spatium {

// Cartesian product of two spaces.
//
// A factor is anything core/access.hpp can reach -- by its members, by
// free functions found by ADL, or by derivation -- so a space someone else
// wrote is a factor as readily as one of ours. Every operation of the
// product is its factors' operation, through the customization points:
// distance sqrt(d1^2 + d2^2), exp and log componentwise, the metric the
// sum of the two, projection part by part.
//
// Points: when both factors' points and tangents are Vecs, a product point
// is one Vec of the combined size, [0, A1) the first factor, [A1, A1 + A2)
// the second -- what the mesh and viewer code reads. Any other point type
// (a matrix, for SPD's affine-invariant form) makes the product's points a
// Pair carrying the vector arithmetic generic algorithms use on tangents.
//
// It used to ask both factors for the member-based MetricSpace and for Vec
// points, so neither a matrix-valued space nor one reached by ADL could be
// a factor; the connectivity matrix showed both, as products of spaces
// green on their own and red together.

namespace product_detail {

// A Vec, exactly: a Matrix stores its entries in `data` as well, and
// asking for that member alone took a matrix for a vector.
template<class P> struct is_vec : std::false_type {};
template<class T, std::size_t N> struct is_vec<Vec<T, N>> : std::true_type {};
template<class P>
concept VecLike = is_vec<std::remove_cvref_t<P>>::value;

template<class P>
constexpr std::size_t size_of() {
    if constexpr (VecLike<P>) return std::tuple_size_v<decltype(P{}.data)>;
    else return 0;
}

// Two components with the vector arithmetic tangents need: sums, scaling,
// negation, comparison.
template<class A, class B>
struct Pair {
    A first{};
    B second{};
    friend constexpr Pair operator+(const Pair& x, const Pair& y) { return {A(x.first + y.first), B(x.second + y.second)}; }
    friend constexpr Pair operator-(const Pair& x, const Pair& y) { return {A(x.first - y.first), B(x.second - y.second)}; }
    friend constexpr Pair operator-(const Pair& x) { return {A(-x.first), B(-x.second)}; }
    template<class S> requires requires(const A& a, S s) { A(a * s); }
    friend constexpr Pair operator*(const Pair& x, S s) { return {A(x.first * s), B(x.second * s)}; }
    template<class S> requires requires(const A& a, S s) { A(a * s); }
    friend constexpr Pair operator*(S s, const Pair& x) { return {A(x.first * s), B(x.second * s)}; }
    friend constexpr bool operator==(const Pair& x, const Pair& y) { return x.first == y.first && x.second == y.second; }
};

// A factor's intrinsic dimension: its own, or, for a space that states
// none, the length of its tangent vectors.
template<class S>
constexpr Dimension dimension_of() {
    if constexpr (requires { S::dimension; }) return Dimension(S::dimension);
    else return Dimension(size_of<spaces::tangent_t<S>>());
}

template<class S>
constexpr bool complete_of() {
    if constexpr (requires { S::is_complete; }) return S::is_complete;
    else return false;
}

}  // namespace product_detail

template<typename S1, typename S2>
    requires std::same_as<spaces::scalar_t<S1>, spaces::scalar_t<S2>>
struct ProductSpace {
    using T = spaces::scalar_t<S1>;
    using ScalarType = T;
    using P1 = spaces::point_t<S1>;
    using P2 = spaces::point_t<S2>;
    using V1 = spaces::tangent_t<S1>;
    using V2 = spaces::tangent_t<S2>;

    // One Vec when every component is one and tangents are the size of
    // points; a Pair otherwise.
    static constexpr std::size_t A1 = product_detail::size_of<P1>();
    static constexpr std::size_t A2 = product_detail::size_of<P2>();
    static constexpr bool flat = product_detail::VecLike<P1> && product_detail::VecLike<P2> &&
                                 product_detail::VecLike<V1> && product_detail::VecLike<V2> &&
                                 product_detail::size_of<V1>() == A1 && product_detail::size_of<V2>() == A2;
    static constexpr std::size_t total_ambient = A1 + A2;

    using PointType = std::conditional_t<flat, Vec<T, A1 + A2>, product_detail::Pair<P1, P2>>;
    using TangentVector = std::conditional_t<flat, Vec<T, A1 + A2>, product_detail::Pair<V1, V2>>;

    // finite factors add; an infinite one makes the product infinite, a dynamic one dynamic
    static constexpr Dimension dimension =
        product_detail::dimension_of<S1>() + product_detail::dimension_of<S2>();
    static constexpr bool is_complete = product_detail::complete_of<S1>() && product_detail::complete_of<S2>();

    S1 space1;
    S2 space2;

    // Split and join. For a flat product the same three serve points and
    // tangents, which share a type there.
    template<class X, class A, class B>
    static X join_as(const A& a, const B& b) {
        if constexpr (flat) {
            X p;
            for (std::size_t i = 0; i < A1; ++i) p[i] = a[i];
            for (std::size_t i = 0; i < A2; ++i) p[A1 + i] = b[i];
            return p;
        } else {
            return X{a, b};
        }
    }
    template<class Part, class X>
    static Part first_of(const X& x) {
        if constexpr (flat) {
            Part r;
            for (std::size_t i = 0; i < A1; ++i) r[i] = x[i];
            return r;
        } else {
            return x.first;
        }
    }
    template<class Part, class X>
    static Part second_of(const X& x) {
        if constexpr (flat) {
            Part r;
            for (std::size_t i = 0; i < A2; ++i) r[i] = x[A1 + i];
            return r;
        } else {
            return x.second;
        }
    }

    P1 first(const PointType& p) const { return first_of<P1>(p); }
    P2 second(const PointType& p) const { return second_of<P2>(p); }
    PointType join(const P1& a, const P2& b) const { return join_as<PointType>(a, b); }

    // TopologicalSpace: in both factors; a factor with no test of its own
    // takes every point.
    bool contains(const PointType& p) const {
        const auto in = [](const auto& s, const auto& x) {
            if constexpr (requires { s.contains(x); }) return s.contains(x);
            else return true;
        };
        return in(space1, first(p)) && in(space2, second(p));
    }

    // MetricSpace: d = sqrt(d1^2 + d2^2), each factor's distance its own or
    // derived.
    ScalarType distance(const PointType& a, const PointType& b) const
        requires spaces::Distanced<S1> && spaces::Distanced<S2>
    {
        using std::sqrt;
        const T d1 = spaces::distance(space1, first(a), first(b));
        const T d2 = spaces::distance(space2, second(a), second(b));
        return sqrt(d1 * d1 + d2 * d2);
    }

    // Manifold: componentwise exp and log.
    PointType exp_map(const PointType& p, const TangentVector& v, ScalarType t) const
        requires spaces::Exponential<S1> && spaces::Exponential<S2>
    {
        return join(P1(spaces::exp_map(space1, first(p), first_of<V1>(v), t)),
                    P2(spaces::exp_map(space2, second(p), second_of<V2>(v), t)));
    }

    TangentVector log_map(const PointType& p, const PointType& q) const
        requires spaces::Logarithmic<S1> && spaces::Logarithmic<S2>
    {
        return join_as<TangentVector>(V1(spaces::log_map(space1, first(p), first(q))),
                                      V2(spaces::log_map(space2, second(p), second(q))));
    }

    // RiemannianManifold: the sum of the factors' metrics.
    ScalarType metric_at(const PointType& p, const TangentVector& u, const TangentVector& v) const
        requires spaces::Metrized<S1> && spaces::Metrized<S2>
    {
        return spaces::metric(space1, first(p), first_of<V1>(u), first_of<V1>(v)) +
               spaces::metric(space2, second(p), second_of<V2>(u), second_of<V2>(v));
    }

    // Componentwise projection: each part by its own space's project, a
    // part whose space has none -- one where every coordinate is a point,
    // like SPD's log-Euclidean chart -- left as it is.
    PointType project(const PointType& p) const
        requires requires(const S1& s, const P1& x) { spaces::project(s, x); } ||
                 requires(const S2& s, const P2& x) { spaces::project(s, x); }
    {
        const auto part = [](const auto& space, const auto& x) {
            using X = std::remove_cvref_t<decltype(x)>;
            if constexpr (requires { spaces::project(space, x); }) return X(spaces::project(space, x));
            else return x;
        };
        return join(part(space1, first(p)), part(space2, second(p)));
    }

    // Surface: componentwise normal, for a flat product of Surfaces.
    TangentVector normal(const PointType& p) const
        requires flat && Surface<S1> && Surface<S2>
    {
        return join_as<TangentVector>(space1.normal(first(p)), space2.normal(second(p)));
    }
};

} // namespace spatium
