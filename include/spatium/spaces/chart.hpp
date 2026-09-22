#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/concepts.hpp>
#  include <spatium/spaces/parametric.hpp>
#  include <spatium/spaces/sphere.hpp>
#  include <cmath>
#  include <concepts>
#endif

// A chart is a (u,v) parametrization of a piece of R^3, and it is what
// every operation in the scene DSL actually consumes -- not `Surface`.
//
// That claim is checkable, and it is the reason this header exists.
// Reading the three operations:
//
//   offset_surface()        evaluate, normal_at, domain, periodic_u/v
//   sample_surface_uniform() area_element, domain
//   parametric_mesh()       evaluate, domain
//
// None of them calls project(), normal() or exp_map(). So `Surface` is
// the wrong name for the requirement: generalising the DSL to it would
// hand the operations methods nobody calls while losing the ones they
// need. The requirement they share had no name at all, and travelled
// through the signatures as the concrete type `ParametricSurface<T>`,
// which is why nothing else could enter the DSL. Naming it is the whole
// of this change.
//
// This is deliberately NOT part of the Set -> Manifold -> Surface
// hierarchy in core/concepts.hpp. That hierarchy says what a space *is*;
// a chart is a choice of coordinates *on* one, and a space can have many
// or none. Sphere<2, T> is a perfectly good RiemannianManifold and has no
// evaluate(u, v) anywhere -- chart_of() below is its first chart, not a
// generalisation of an existing one.
//
// The extension point is `chart_of(space)`, found by ADL: the same shape
// as `point_to(p, surf)` in the contact concepts, and for the same
// reason. A new space joins the DSL by having an overload written for it
// in its own namespace, with no edit to this header and no runtime
// registry -- the dispatch is a compile-time overload resolution, so a
// chart that a caller never uses costs nothing, and one it does use
// inlines.

SPATIUM_EXPORT namespace spatium {

// What the DSL operations require. `ParametricSurface<T>` satisfies it
// today; it is the erased form every chart is converted into, because a
// trace holds nodes of one type and cannot be templated on each node's
// space.
template<typename C>
concept Chart =
    requires { typename C::ScalarType; } &&
    Scalar<typename C::ScalarType> &&
    requires(const C& c, typename C::ScalarType s) {
        { c.evaluate(s, s) }     -> std::convertible_to<Vec<typename C::ScalarType, 3>>;
        { c.normal_at(s, s) }    -> std::convertible_to<Vec<typename C::ScalarType, 3>>;
        { c.area_element(s, s) } -> std::convertible_to<typename C::ScalarType>;
        { c.domain().u_min }     -> std::convertible_to<typename C::ScalarType>;
        { c.domain().u_max }     -> std::convertible_to<typename C::ScalarType>;
        { c.domain().v_min }     -> std::convertible_to<typename C::ScalarType>;
        { c.domain().v_max }     -> std::convertible_to<typename C::ScalarType>;
        { c.periodic_u() }       -> std::convertible_to<bool>;
        { c.periodic_v() }       -> std::convertible_to<bool>;
    };

// A ParametricSurface is already a chart; this is the identity that lets
// the DSL take one through the same door as everything else.
template<Scalar T>
ParametricSurface<T> chart_of(const ParametricSurface<T>& surface) {
    return surface;
}

// Sphere<2, T>'s first parametrization: the standard spherical chart,
// azimuth u in [0, 2pi) around the z axis, polar angle v in [0, pi] from
// the north pole.
//
// Note which sphere this is. Sphere<N, T> is the N-sphere embedded in
// R^{N+1}, so the ordinary sphere in R^3 -- the only one a (u,v) -> R^3
// chart can describe -- is Sphere<2, T>. Sphere<3, T> lives in R^4 and
// gets no overload here, which is not an omission: there is no chart of
// that shape to write.
//
// periodic_u is true (the seam at u = 0 is a seam, nothing more) and
// periodic_v is false, which is the honest answer and not the whole
// story: the v edges do not bound the surface, they collapse to the two
// poles. That is the middle row of the chart-edge table in ROADMAP --
// closed, but not by periodicity -- and it is exactly why is_closed()
// tests the geometry rather than reading these two flags. The poles are
// also where this chart degenerates: area_element goes to zero there (it
// is r^2 sin v), and the finite-difference normal_at with it. A pole is
// a property of the chart, not a defect of the sphere.
template<Scalar T>
ParametricSurface<T> chart_of(const Sphere<2, T>& sphere) {
    using std::sin, std::cos, std::acos;
    const T r  = sphere.radius;
    const T pi = acos(T{-1});
    // The closed forms are handed over because this chart degenerates at
    // its poles and the generic search cannot converge there -- see
    // `ParametricSurface::with_closed_forms`. A sphere's projection and
    // normal are one line each, so searching for them was always the
    // wrong trade even where the search happened to work.
    return ParametricSurface<T>(
               [r](T u, T v) -> Vec<T, 3> {
                   return {r * sin(v) * cos(u), r * sin(v) * sin(u), r * cos(v)};
               },
               {T{0}, T{2} * pi, T{0}, pi},
               /*periodic_u=*/true, /*periodic_v=*/false)
        .with_closed_forms(
            [r](const Vec<T, 3>& p) -> Vec<T, 3> {
                const T n = p.norm();
                // The centre projects nowhere in particular; the pole is
                // as good an answer as any and keeps this total.
                if (n <= epsilon<T>()) return {T{0}, T{0}, r};
                return Vec<T, 3>{p * (r / n)};
            },
            [](const Vec<T, 3>& p) -> Vec<T, 3> {
                const T n = p.norm();
                if (n <= epsilon<T>()) return {T{0}, T{0}, T{1}};
                return Vec<T, 3>{p * (T{1} / n)};
            });
}

// Anything the DSL can take: a chart, or a space that has one written
// for it. Unqualified so a user's overload in the user's own namespace
// is found by ADL at the point of instantiation.
template<typename S, typename T>
concept Chartable = requires(const S& s) {
    { chart_of(s) } -> std::convertible_to<ParametricSurface<T>>;
};

} // namespace spatium
