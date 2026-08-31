#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <algorithm>
#  include <cmath>
#endif

SPATIUM_EXPORT namespace spatium {

// N-sphere embedded in R^{N+1}.
// Sphere<2> = standard 2-sphere (surface of a ball in 3D).
// Points are Vec<T, N+1> on the sphere surface.

template<std::size_t N, Scalar T = double>
struct Sphere {
    using ScalarType    = T;
    using PointType     = Vec<T, N + 1>;
    using TangentVector = Vec<T, N + 1>;

    static constexpr std::size_t dimension = N;
    static constexpr bool is_complete = true;

    T radius = T{1};

    bool contains(const PointType& p) const {
        using std::abs;
        auto diff = p.norm_squared() - radius * radius;
        return abs(diff) < epsilon<T>() * radius * radius;
    }

    ScalarType distance(const PointType& a, const PointType& b) const {
        using std::acos; using std::clamp;
        auto cos_angle = clamp(a.dot(b) / (radius * radius), T{-1}, T{1});
        return radius * acos(cos_angle);
    }

    PointType exp_map(const PointType& p, const TangentVector& v, ScalarType t) const {
        using std::cos; using std::sin;
        auto tv = v * t;
        auto theta = tv.norm() / radius;
        if (theta < epsilon<T>()) return p;
        auto dir = tv / (theta * radius);
        return p * cos(theta) + dir * (radius * sin(theta));
    }

    TangentVector log_map(const PointType& p, const PointType& q) const {
        using std::acos; using std::clamp;
        auto cos_angle = clamp(p.dot(q) / (radius * radius), T{-1}, T{1});
        auto theta = acos(cos_angle);
        if (theta < epsilon<T>()) return TangentVector{};
        auto proj = q - p * (p.dot(q) / p.dot(p));
        auto proj_norm = proj.norm();
        if (proj_norm < epsilon<T>()) {
            // Antipodal: pick canonical tangent direction perpendicular to p
            TangentVector canonical{};
            for (std::size_t i = 0; i <= N; ++i) {
                using std::abs;
                if (abs(p[i]) < T{0.9} * radius) {
                    canonical[i] = T{1};
                    break;
                }
            }
            // Gram-Schmidt: remove p component
            canonical = canonical - p * (p.dot(canonical) / p.dot(p));
            auto cn = canonical.norm();
            if (cn < epsilon<T>()) return TangentVector{};
            return canonical * (theta / cn);
        }
        return proj * (theta / proj_norm);
    }

    constexpr ScalarType metric_at(const PointType&,
                                   const TangentVector& u,
                                   const TangentVector& v) const {
        return u.dot(v);
    }

    PointType project(const PointType& p) const {
        auto n = p.norm();
        if (n < epsilon<T>()) {
            PointType result{};
            result[N] = radius;
            return result;
        }
        return p * (radius / n);
    }

    TangentVector normal(const PointType& p) const {
        return p / p.norm();
    }

    TangentVector project_tangent(const PointType& p, const TangentVector& v) const {
        auto n = normal(p);
        return v - n * n.dot(v);
    }
};

static_assert(RiemannianManifold<Sphere<2>>);
static_assert(Surface<Sphere<2>>);
static_assert(Complete<Sphere<2>>);
static_assert(RiemannianManifold<Sphere<1>>);
static_assert(Surface<Sphere<3>>);

using S1 = Sphere<1>;
using S2 = Sphere<2>;
using S3 = Sphere<3>;

} // namespace spatium
