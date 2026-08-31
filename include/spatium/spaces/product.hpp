#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <cmath>
#  include <tuple>
#endif

SPATIUM_EXPORT namespace spatium {

// Cartesian product of two spaces.
// Points are stored as a single Vec of combined dimension for mesh/viewer compat.
// Components split/join at indices [0, D1) and [D1, D1+D2).

template<typename S1, typename S2>
    requires MetricSpace<S1> && MetricSpace<S2>
          && std::same_as<typename S1::ScalarType, typename S2::ScalarType>
struct ProductSpace {
    using T = typename S1::ScalarType;
    using ScalarType = T;
    // Ambient dimensions from PointType array size (safe, no sizeof tricks)
    static constexpr std::size_t A1 = std::tuple_size_v<decltype(typename S1::PointType{}.data)>;
    static constexpr std::size_t A2 = std::tuple_size_v<decltype(typename S2::PointType{}.data)>;
    static constexpr std::size_t total_ambient = A1 + A2;

    using PointType = Vec<T, total_ambient>;
    using TangentVector = Vec<T, total_ambient>;

    static constexpr std::size_t dimension = S1::dimension + S2::dimension;
    static constexpr bool is_complete = S1::is_complete && S2::is_complete;

    S1 space1;
    S2 space2;

    // Split combined point into components
    typename S1::PointType first(const PointType& p) const {
        typename S1::PointType r;
        for (std::size_t i = 0; i < A1; ++i) r[i] = p[i];
        return r;
    }

    typename S2::PointType second(const PointType& p) const {
        typename S2::PointType r;
        for (std::size_t i = 0; i < A2; ++i) r[i] = p[A1 + i];
        return r;
    }

    // Join two component points into combined
    PointType join(const typename S1::PointType& a, const typename S2::PointType& b) const {
        PointType p;
        for (std::size_t i = 0; i < A1; ++i) p[i] = a[i];
        for (std::size_t i = 0; i < A2; ++i) p[A1 + i] = b[i];
        return p;
    }

    // TopologicalSpace
    bool contains(const PointType& p) const {
        return space1.contains(first(p)) && space2.contains(second(p));
    }

    // MetricSpace: product metric d = sqrt(d1^2 + d2^2)
    ScalarType distance(const PointType& a, const PointType& b) const {
        using std::sqrt;
        auto d1 = space1.distance(first(a), first(b));
        auto d2 = space2.distance(second(a), second(b));
        return sqrt(d1 * d1 + d2 * d2);
    }

    // Manifold: componentwise exp/log
    PointType exp_map(const PointType& p, const TangentVector& v, ScalarType t) const
        requires Manifold<S1> && Manifold<S2>
    {
        typename S1::TangentVector v1;
        for (std::size_t i = 0; i < A1; ++i) v1[i] = v[i];
        typename S2::TangentVector v2;
        for (std::size_t i = 0; i < A2; ++i) v2[i] = v[A1 + i];

        auto r1 = space1.exp_map(first(p), v1, t);
        auto r2 = space2.exp_map(second(p), v2, t);
        return join(r1, r2);
    }

    TangentVector log_map(const PointType& p, const PointType& q) const
        requires Manifold<S1> && Manifold<S2>
    {
        auto v1 = space1.log_map(first(p), first(q));
        auto v2 = space2.log_map(second(p), second(q));
        TangentVector v;
        for (std::size_t i = 0; i < A1; ++i) v[i] = v1[i];
        for (std::size_t i = 0; i < A2; ++i) v[A1 + i] = v2[i];
        return v;
    }

    // RiemannianManifold
    ScalarType metric_at(const PointType& p, const TangentVector& u, const TangentVector& v) const
        requires RiemannianManifold<S1> && RiemannianManifold<S2>
    {
        typename S1::TangentVector u1, v1;
        typename S2::TangentVector u2, v2;
        for (std::size_t i = 0; i < A1; ++i) { u1[i] = u[i]; v1[i] = v[i]; }
        for (std::size_t i = 0; i < A2; ++i) { u2[i] = u[A1 + i]; v2[i] = v[A1 + i]; }
        return space1.metric_at(first(p), u1, v1) + space2.metric_at(second(p), u2, v2);
    }

    // Surface: componentwise project/normal
    PointType project(const PointType& p) const
        requires Surface<S1> && Surface<S2>
    {
        return join(space1.project(first(p)), space2.project(second(p)));
    }

    TangentVector normal(const PointType& p) const
        requires Surface<S1> && Surface<S2>
    {
        auto n1 = space1.normal(first(p));
        auto n2 = space2.normal(second(p));
        TangentVector n;
        for (std::size_t i = 0; i < A1; ++i) n[i] = n1[i];
        for (std::size_t i = 0; i < A2; ++i) n[A1 + i] = n2[i];
        return n;
    }
};

} // namespace spatium
