#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium {

template<std::size_t N, Scalar T = double>
struct Euclidean {
    using ScalarType    = T;
    using PointType     = Vec<T, N>;
    using VectorType    = Vec<T, N>;
    using TangentVector = Vec<T, N>;

    static constexpr std::size_t dimension = N;
    static constexpr bool is_complete = true;

    // TopologicalSpace: R^N contains everything
    constexpr bool contains(const PointType&) const { return true; }

    // MetricSpace
    ScalarType distance(const PointType& a, const PointType& b) const {
        return (a - b).norm();
    }

    // NormedSpace
    ScalarType norm(const VectorType& v) const {
        return v.norm();
    }

    // InnerProductSpace
    constexpr ScalarType inner(const VectorType& u, const VectorType& v) const {
        return u.dot(v);
    }

    // Manifold: trivial (flat space)
    constexpr PointType exp_map(const PointType& p, const TangentVector& v, ScalarType t) const {
        return p + v * t;
    }

    constexpr TangentVector log_map(const PointType& p, const PointType& q) const {
        return q - p;
    }

    // RiemannianManifold: flat metric = standard inner product
    constexpr ScalarType metric_at(const PointType&,
                                   const TangentVector& u,
                                   const TangentVector& v) const {
        return u.dot(v);
    }

    // Surface: identity projection (flat space, all points are "on surface")
    constexpr PointType project(const PointType& p) const { return p; }

    // Surface: constant normal (last basis vector for N>=3, zero for N<3)
    constexpr TangentVector normal(const PointType&) const {
        TangentVector n{};
        if constexpr (N >= 3) n[2] = T{1};
        else if constexpr (N == 2) n[1] = T{1};
        else n[0] = T{1};
        return n;
    }
};

// Verify concept satisfaction
static_assert(Set<Euclidean<3>>);
static_assert(TopologicalSpace<Euclidean<3>>);
static_assert(MetricSpace<Euclidean<3>>);
static_assert(VectorSpace<Euclidean<3>>);
static_assert(NormedSpace<Euclidean<3>>);
static_assert(InnerProductSpace<Euclidean<3>>);
static_assert(Complete<Euclidean<3>>);
static_assert(BanachSpace<Euclidean<3>>);
static_assert(HilbertSpace<Euclidean<3>>);
static_assert(EuclideanSpace<Euclidean<3>>);
static_assert(Manifold<Euclidean<3>>);
static_assert(RiemannianManifold<Euclidean<3>>);
static_assert(Surface<Euclidean<3>>);

// Works for any dimension
static_assert(EuclideanSpace<Euclidean<1>>);
static_assert(EuclideanSpace<Euclidean<2>>);
static_assert(EuclideanSpace<Euclidean<100>>);

// Works with float
static_assert(EuclideanSpace<Euclidean<3, float>>);

// Common aliases
using E2 = Euclidean<2>;
using E3 = Euclidean<3>;
using E4 = Euclidean<4>;

} // namespace spatium
