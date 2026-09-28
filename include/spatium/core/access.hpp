#pragma once
// Customization points for spaces: the std::ranges::begin pattern, applied
// to what a space is rather than to what a range is.
//
// Every operation here is an object, not a function, and finds its
// implementation in a fixed order:
//
//   1. a member        space.exp_map(p, v, t)       -- every space in spaces/
//   2. a free function exp_map(space, p, v, t)      -- found by ADL, so a type
//                                                      the author cannot edit
//                                                      can still be a space
//   3. a derivation    from the operations that do exist
//
// Step 3 is the point. A space states a minimal basis -- exp, log and a
// metric -- and whatever follows from that basis is supplied, not
// re-implemented per space:
//
//   distance(p, q)      = |log_p q|_g
//   norm_at(p, v)       = sqrt(g_p(v, v))
//   geodesic(p, q, t)   = exp_p(t log_p q)
//   midpoint(p, q)      = geodesic(p, q, 1/2)
//
// A derivation never overrides an operation the space provides, so a closed
// form (Sphere's acos, Hyperbolic's acosh) is always the one used; the
// derived form is what a space gets when it has none, and what a test holds
// the closed form against.
//
// The member-based concepts in core/concepts.hpp are unchanged: every type
// that satisfied them satisfies the ones below. What changes is who else
// can -- see tests/test_space_access.cpp for a space with no members at all.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <cmath>
#  include <concepts>
#endif

SPATIUM_EXPORT namespace spatium::spaces {

// ── Associated types ───────────────────────────────────────────
// Read from the space's nested names by default. A type that cannot be
// edited specialises space_traits instead. The primary template is empty
// so that a type with neither is simply not a space -- a concept returns
// false instead of the compiler stopping inside a class template.

template<class S>
struct space_traits {};

template<class S>
    requires requires {
        typename S::PointType;
        typename S::TangentVector;
        typename S::ScalarType;
    }
struct space_traits<S> {
    using point_type   = typename S::PointType;
    using tangent_type = typename S::TangentVector;
    using scalar_type  = typename S::ScalarType;
};

template<class S> using point_t   = typename space_traits<S>::point_type;
template<class S> using tangent_t = typename space_traits<S>::tangent_type;
template<class S> using scalar_t  = typename space_traits<S>::scalar_type;

namespace access_detail {

// Poison pills. An unqualified call inside a customization point must find
// the user's overload by ADL and nothing else -- not the customization
// point object itself, and not an unrelated spatium:: function that happens
// to share the name (algebra's distance(Vec, Vec), for one).
void exp_map() = delete;
void log_map() = delete;
void metric_at() = delete;
void distance() = delete;
void project() = delete;
void normal() = delete;

struct exp_map_fn {
    template<class S, class P, class V, class T>
    constexpr auto operator()(const S& s, const P& p, const V& v, const T& t) const
        requires requires { s.exp_map(p, v, t); } || requires { exp_map(s, p, v, t); }
    {
        if constexpr (requires { s.exp_map(p, v, t); }) return s.exp_map(p, v, t);
        else                                            return exp_map(s, p, v, t);
    }
    // t = 1: the point a tangent vector reaches.
    template<class S, class P, class V>
    constexpr auto operator()(const S& s, const P& p, const V& v) const
        requires requires { s.exp_map(p, v, scalar_t<S>{1}); }
              || requires { exp_map(s, p, v, scalar_t<S>{1}); }
    {
        const scalar_t<S> one{1};
        if constexpr (requires { s.exp_map(p, v, one); }) return s.exp_map(p, v, one);
        else                                              return exp_map(s, p, v, one);
    }
};

struct log_map_fn {
    template<class S, class P>
    constexpr auto operator()(const S& s, const P& p, const P& q) const
        requires requires { s.log_map(p, q); } || requires { log_map(s, p, q); }
    {
        if constexpr (requires { s.log_map(p, q); }) return s.log_map(p, q);
        else                                         return log_map(s, p, q);
    }
};

struct metric_fn {
    template<class S, class P, class V>
    constexpr auto operator()(const S& s, const P& p, const V& u, const V& v) const
        requires requires { s.metric_at(p, u, v); } || requires { metric_at(s, p, u, v); }
    {
        if constexpr (requires { s.metric_at(p, u, v); }) return s.metric_at(p, u, v);
        else                                              return metric_at(s, p, u, v);
    }
};

struct project_fn {
    template<class S, class P>
    constexpr auto operator()(const S& s, const P& p) const
        requires requires { s.project(p); } || requires { project(s, p); }
    {
        if constexpr (requires { s.project(p); }) return s.project(p);
        else                                      return project(s, p);
    }
};

struct normal_fn {
    template<class S, class P>
    constexpr auto operator()(const S& s, const P& p) const
        requires requires { s.normal(p); } || requires { normal(s, p); }
    {
        if constexpr (requires { s.normal(p); }) return s.normal(p);
        else                                     return normal(s, p);
    }
};

} // namespace access_detail

// The objects. `inline namespace cpo` keeps them from colliding with a
// hidden friend of the same name in a user's type (the ranges idiom).
inline namespace cpo {
inline constexpr access_detail::exp_map_fn exp_map{};
inline constexpr access_detail::log_map_fn log_map{};
inline constexpr access_detail::metric_fn  metric{};
inline constexpr access_detail::project_fn project{};
inline constexpr access_detail::normal_fn  normal{};
} // namespace cpo

// ── Derived operations ─────────────────────────────────────────

namespace access_detail {

struct norm_at_fn {
    template<class S, class P, class V>
    constexpr auto operator()(const S& s, const P& p, const V& v) const
        requires requires { spaces::metric(s, p, v, v); }
    {
        using std::sqrt;
        return sqrt(spaces::metric(s, p, v, v));
    }
};

struct distance_fn {
    template<class S, class P>
    constexpr auto operator()(const S& s, const P& p, const P& q) const
        requires requires { s.distance(p, q); }
              || requires { distance(s, p, q); }
              || requires { spaces::metric(s, p, spaces::log_map(s, p, q), spaces::log_map(s, p, q)); }
    {
        if constexpr (requires { s.distance(p, q); })   return s.distance(p, q);
        else if constexpr (requires { distance(s, p, q); }) return distance(s, p, q);
        else {
            // |log_p q|_g: exact wherever log is -- inside the injectivity
            // radius, which is where a minimising geodesic is unique.
            using std::sqrt;
            const auto v = spaces::log_map(s, p, q);
            return sqrt(spaces::metric(s, p, v, v));
        }
    }
};

} // namespace access_detail

inline namespace cpo {
inline constexpr access_detail::norm_at_fn  norm_at{};
inline constexpr access_detail::distance_fn distance{};
} // namespace cpo

// ── Capabilities ───────────────────────────────────────────────
// What a space can be asked, whichever of the three routes answers. Named
// apart from core/concepts.hpp's member-based HasExpMap/MetricSpace/... so
// the two never read as the same requirement.

template<class S>
concept Exponential = requires(const S& s, const point_t<S>& p, const tangent_t<S>& v,
                               scalar_t<S> t) {
    { spaces::exp_map(s, p, v, t) } -> std::convertible_to<point_t<S>>;
};

template<class S>
concept Logarithmic = requires(const S& s, const point_t<S>& p, const point_t<S>& q) {
    { spaces::log_map(s, p, q) } -> std::convertible_to<tangent_t<S>>;
};

template<class S>
concept Metrized = requires(const S& s, const point_t<S>& p, const tangent_t<S>& u,
                            const tangent_t<S>& v) {
    { spaces::metric(s, p, u, v) } -> std::convertible_to<scalar_t<S>>;
};

template<class S>
concept Distanced = requires(const S& s, const point_t<S>& p, const point_t<S>& q) {
    { spaces::distance(s, p, q) } -> std::convertible_to<scalar_t<S>>;
};

// The minimal basis of a Riemannian manifold. distance follows from it, so
// it is not asked for.
template<class S>
concept Riemannian = Exponential<S> && Logarithmic<S> && Metrized<S>;

// ── Algorithms that need only the basis ────────────────────────

// The point a fraction t of the way along the geodesic from p to q.
template<class S>
    requires Exponential<S> && Logarithmic<S>
constexpr point_t<S> geodesic(const S& s, const point_t<S>& p, const point_t<S>& q,
                              scalar_t<S> t) {
    return spaces::exp_map(s, p, tangent_t<S>{spaces::log_map(s, p, q)}, t);
}

template<class S>
    requires Exponential<S> && Logarithmic<S>
constexpr point_t<S> midpoint(const S& s, const point_t<S>& p, const point_t<S>& q) {
    return geodesic(s, p, q, scalar_t<S>{1} / scalar_t<S>{2});
}

} // namespace spatium::spaces
