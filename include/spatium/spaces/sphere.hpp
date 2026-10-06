#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/algebra/vector.hpp>
#  include <spatium/spaces/cos_sinc.hpp>
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

    // The angle is t |v| / r, signed and linear in t, and the direction
    // v / |v| does not depend on t: exp is smooth in t through t = 0, so a
    // Dual t carries d/dt exp_p(t v) = v there. Measuring |t v| instead (as
    // this did) is not differentiable at t = 0, and its small-angle branch,
    // which compares only a Dual's value, returned p with a zero derivative
    // -- found by the connectivity matrix's derivative cells.
    PointType exp_map(const PointType& p, const TangentVector& v, ScalarType t) const {
        // p cos(theta) + v t sinc(theta), theta^2 = t^2 |v|^2 / r^2: smooth
        // in v through v = 0 as well (see spaces/cos_sinc.hpp), where the
        // |v| form returned p with a zero derivative for a Dual tangent.
        const T theta2 = t * t * v.dot(v) / (radius * radius);
        const auto [c, sinc] = cos_sinc_of_square(theta2);
        return PointType{p * c + v * (t * sinc)};
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
            return canonical * (radius * theta / cn);
        }
        // The length is the distance along the sphere, radius * theta,
        // because that is what exp_map reads -- it divides by the radius.
        // It was theta alone, the two inverse only on the unit sphere.
        return proj * (radius * theta / proj_norm);
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
