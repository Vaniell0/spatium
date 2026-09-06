#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/core/error.hpp>
#  include <spatium/geometry/concepts.hpp>
#  include <spatium/geometry/convex_hull.hpp>
#  include <spatium/geometry/make.hpp>
#  include <spatium/spaces/sphere.hpp>
#  include <algorithm>
#  include <cassert>
#  include <cmath>
#  include <format>
#  include <numbers>
#  include <type_traits>
#  include <utility>
#  include <vector>
#endif

// ── SphericalPolygon<T> ──────────────────────────────────────────
//
// A polygon on the 2-sphere (spatium::Sphere<2,T>): an ordered list of
// ambient vertices (Vec<T,3>, all at distance `radius` from the
// ambient origin -- exactly Sphere<2,T>'s own point convention, see
// spaces/sphere.hpp), with edges being MINOR-ARC great-circle
// segments between consecutive vertices -- the sphere's own
// geodesics. This is the curved analogue of geometry/polygon.hpp's
// Polygon and geometry/boolean.hpp's Sutherland-Hodgman machinery.
//
// Why the sphere first, and only the sphere: Spatium's differentiator
// is Surface as a first-class concept rather than a hardcoded type
// (see docs/architecture.md), and boolean set-operations on curved
// surfaces are the natural next boolean-ops target after the planar
// case in boolean.hpp -- Google's S2 library already does this for
// the sphere specifically. What makes the sphere tractable in closed
// form is that every geodesic is a great circle, i.e. the sphere's
// intersection with a plane through its own center -- so "clip by a
// half-space" is *exactly* the planar sign test (clip_half_space
// below), and two clip planes meet in a line through the center,
// giving two antipodal candidate crossing points instead of the
// planar case's unique one (detail::great_circle_edge_crossing). A
// general Surface has no such closed form: geodesic-geodesic
// intersection there needs numerical shooting/root-finding, which is
// deliberately out of scope for this file.
//
// RESTRICTION -- convex polygons only. Boolean ops on concave
// spherical polygons need real topological bookkeeping this file does
// not attempt: a concave clip region can require several half-spaces
// combined by OR (not the AND this file composes), can split one
// input polygon into several disjoint outputs, and can cross the same
// great circle more than once -- none of which "clip by each
// edge-plane of the other polygon, in sequence" can represent, convex
// or not. Rather than silently mishandling concave input, every
// operation below assumes and documents convexity, exactly mirroring
// Polygon's own operator&/-/+ restriction to convex input in
// boolean.hpp.

SPATIUM_EXPORT namespace spatium::geometry {

template<Scalar T = double>
struct SphericalPolygon {
    using ScalarType = T;
    using PointType = Vec<T, 3>;
    static constexpr std::size_t ambient_dimension = 3;

    std::vector<PointType> vertices; // ordered, CCW as seen from outside the sphere
    T radius = T{1};

    std::size_t size() const { return vertices.size(); }

    // The minor-arc great-circle edge between consecutive vertices i, i+1.
    struct Edge { PointType a, b; };
    Edge edge(std::size_t i) const {
        return {vertices[i], vertices[(i + 1) % vertices.size()]};
    }

    // Angular length of edge i, in radians (arc length = radius * this).
    T edge_angle(std::size_t i) const {
        using std::acos; using std::clamp;
        auto e = edge(i);
        auto cos_a = clamp(e.a.dot(e.b) / (radius * radius), T{-1}, T{1});
        return acos(cos_a);
    }

    T perimeter() const {
        T total{0};
        for (std::size_t i = 0; i < vertices.size(); ++i)
            total += radius * edge_angle(i);
        return total;
    }

    // Ambient centroid, re-projected onto the sphere. A cheap
    // representative point for the Shape concept, not a claim that
    // this is the spherical (area-weighted) centroid.
    PointType centroid() const {
        PointType sum{};
        for (auto& v : vertices) sum = Vec<T, 3>{sum + v};
        auto n = sum.norm();
        if (n < epsilon<T>()) return PointType{};
        return Vec<T, 3>{sum * (radius / n)};
    }

    // Spherical excess (Girard's theorem): area = (angle sum over
    // vertices - (n-2)*pi) * radius^2 -- the curved analogue of
    // Polygon::measure()'s planar shoelace formula. This is the one
    // real implementation; area() below is a pure alias, per
    // docs/conventions.md's measure()-alias rule.
    T measure() const {
        auto n = vertices.size();
        if (n < 3) return T{0};
        T angle_sum{0};
        for (std::size_t i = 0; i < n; ++i) {
            auto& prev = vertices[(i + n - 1) % n];
            auto& curr = vertices[i];
            auto& next = vertices[(i + 1) % n];
            angle_sum += interior_angle(prev, curr, next);
        }
        auto excess = angle_sum - T(n - 2) * std::numbers::pi_v<T>;
        return excess * radius * radius;
    }
    T area() const { return measure(); }

    // Convex-only point containment: `p` is inside iff it is on the
    // non-negative side of every edge's great-circle plane -- the
    // same "inside every clip half-space" test intersection uses.
    bool contains(const PointType& p) const {
        auto n = vertices.size();
        auto tol = epsilon<T>() * radius * radius * radius;
        for (std::size_t i = 0; i < n; ++i) {
            auto e = edge(i);
            auto normal = Vec<T, 3>{e.a.cross(e.b)};
            if (normal.dot(p) < -tol) return false;
        }
        return true;
    }

private:
    // Interior angle at `curr`, between the geodesics curr->prev and
    // curr->next: project each neighbor's chord onto the tangent
    // plane at `curr` (the exact same "strip off the component along
    // the point" projection Sphere<N,T>::log_map uses to turn an
    // ambient point into a tangent vector) and take the angle between
    // the two tangent directions.
    static T interior_angle(const PointType& prev, const PointType& curr, const PointType& next) {
        using std::acos; using std::clamp;
        auto to_tangent = [&](const PointType& other) {
            auto d = curr.dot(curr);
            return Vec<T, 3>{other - curr * (curr.dot(other) / d)};
        };
        auto t_prev = to_tangent(prev);
        auto t_next = to_tangent(next);
        auto denom = t_prev.norm() * t_next.norm();
        if (denom < epsilon<T>()) return T{0}; // degenerate vertex (antipodal neighbor)
        auto cos_angle = clamp(t_prev.dot(t_next) / denom, T{-1}, T{1});
        return acos(cos_angle);
    }
};

static_assert(std::is_same_v<typename SphericalPolygon<double>::PointType,
                              typename Sphere<2, double>::PointType>,
              "SphericalPolygon<T> vertices use Sphere<2,T>'s own ambient point convention");

static_assert(Shape<SphericalPolygon<>>);
static_assert(ClosedShape<SphericalPolygon<>>);
static_assert(Measurable<SphericalPolygon<>>);

// internal -- do not use, no API stability. Great-circle-plane
// intersection and arc-membership helpers shared by clip_half_space,
// spherical_intersection and spherical_difference, plus the gnomonic
// projection pair spherical_union uses to reuse the existing planar
// convex_hull().
namespace detail {

// Two planes through the sphere's center (normals n1, n2, assumed non-
// parallel) meet in a line through the center; that line meets the
// sphere at two antipodal points. Returns one of them (the caller
// picks the sign it needs, e.g. via on_minor_arc below). Parallel
// normals mean the "two" great circles are actually the same one (or
// disjoint copies, impossible for planes through a common center) --
// an internal invariant violation for every call site below, not a
// reachable boundary-input error, so this asserts rather than
// returning Result<T>.
template<Scalar T>
Vec<T, 3> great_circle_crossing_point(const Vec<T, 3>& n1, const Vec<T, 3>& n2, T radius) {
    auto line_dir = Vec<T, 3>{n1.cross(n2)};
    auto len = line_dir.norm();
    assert(len > epsilon<T>() * radius * radius && "great circles are parallel or coincident");
    return Vec<T, 3>{line_dir * (radius / len)};
}

// Is `x` -- already known to lie on the great circle through `a` and
// `b` -- on the minor arc from `a` to `b`? Same cross-product
// sidedness test as a planar point-on-segment check (detail::
// sutherland_hodgman's is_inside in boolean.hpp), generalized to
// angles about the arc's own great-circle axis a.cross(b): x is
// between a and b iff both a->x and x->b agree in rotation sign with
// a->b itself.
//
// `endpoints_inclusive` matters for spherical_difference's own
// edge-vs-edge crossing search: it must accept a candidate that lands
// exactly on the *other* edge's endpoint (its q1/q2), but reject one
// that merely reproduces one of *this* edge's own endpoints (p1/p2) --
// the latter isn't a real transition, just a source vertex that
// happens to sit exactly on the other polygon's boundary. Passing
// `false` for the edge being walked (and leaving the other edge at
// its default `true`) is exactly that asymmetric check; clip_half_space
// always wants the inclusive default, since there the candidate is
// only ever tested against the one edge it came from.
template<Scalar T>
bool on_minor_arc(const Vec<T, 3>& a, const Vec<T, 3>& b, const Vec<T, 3>& x, T radius,
                   bool endpoints_inclusive = true) {
    auto axis = Vec<T, 3>{a.cross(b)};
    auto tol = epsilon<T>() * radius * radius * radius * radius;
    auto s1 = Vec<T, 3>{a.cross(x)}.dot(axis);
    auto s2 = Vec<T, 3>{x.cross(b)}.dot(axis);
    if (endpoints_inclusive) return s1 >= -tol && s2 >= -tol;
    return s1 > tol && s2 > tol;
}

// The one genuinely new piece beyond the planar Sutherland-Hodgman
// case: intersect a great-circle edge (the minor arc prev->curr) with
// the clipping great circle (the plane through the center with normal
// `clip_normal`). Disambiguates the two antipodal candidates by
// picking whichever lands on the prev->curr arc.
template<Scalar T>
Vec<T, 3> great_circle_edge_crossing(const Vec<T, 3>& prev, const Vec<T, 3>& curr,
                                      const Vec<T, 3>& clip_normal, T radius) {
    auto edge_normal = Vec<T, 3>{prev.cross(curr)};
    auto candidate = great_circle_crossing_point(edge_normal, clip_normal, radius);
    if (on_minor_arc(prev, curr, candidate, radius)) return candidate;
    return Vec<T, 3>{-candidate};
}

// Orthonormal basis for the tangent plane at unit vector `n` -- the
// same "seed with a non-parallel coordinate axis, then Gram-Schmidt"
// trick used throughout this codebase (Circle::project(),
// Sphere::log_map()'s antipodal fallback).
template<Scalar T>
std::pair<Vec<T, 3>, Vec<T, 3>> tangent_basis(const Vec<T, 3>& n) {
    Vec<T, 3> seed{};
    for (std::size_t i = 0; i < 3; ++i) {
        if (std::abs(n[i]) < T{0.9}) { seed[i] = T{1}; break; }
    }
    auto e1 = Vec<T, 3>{seed - n * n.dot(seed)}.normalized();
    auto e2 = Vec<T, 3>{n.cross(e1)};
    return {e1, e2};
}

// Gnomonic (central) projection of a sphere point onto the tangent
// plane at `pole`, in that plane's local 2D basis (e1, e2). Gnomonic
// projection's defining property is that it maps every great circle
// to a straight line -- so a spherically-convex region contained in
// the open hemisphere centered at `pole` becomes an ordinary
// planar-convex region under this map. Exact only within that open
// hemisphere (dot(p/radius, pole) > 0); spherical_union checks that
// precondition before calling this.
template<Scalar T>
Vec<T, 2> gnomonic_project(const Vec<T, 3>& p, const Vec<T, 3>& pole,
                            const Vec<T, 3>& e1, const Vec<T, 3>& e2, T radius) {
    auto u = Vec<T, 3>{p * (T{1} / radius)};
    auto denom = u.dot(pole);
    auto proj = Vec<T, 3>{u * (radius / denom)}; // point on the tangent plane
    auto local = Vec<T, 3>{proj - pole * radius};
    return {local.dot(e1), local.dot(e2)};
}

// Inverse of gnomonic_project: a point in the tangent plane's local
// basis, mapped back to the sphere by normalizing its ambient
// direction from the center and rescaling to `radius`.
template<Scalar T>
Vec<T, 3> gnomonic_unproject(const Vec<T, 2>& q, const Vec<T, 3>& pole,
                              const Vec<T, 3>& e1, const Vec<T, 3>& e2, T radius) {
    auto on_plane = Vec<T, 3>{pole * radius + e1 * q[0] + e2 * q[1]};
    return Vec<T, 3>{on_plane.normalized() * radius};
}

// Remove consecutive (wraparound included) near-duplicate points from a
// vertex loop. Needed because a source vertex sitting exactly on a
// clip plane produces a "crossing point" that coincides with the
// vertex itself, and the Sutherland-Hodgman push logic (mirroring
// boolean.hpp's own planar version, which has the same latent gap)
// pushes both the crossing and the vertex separately -- a real
// duplicate-vertex edge case first exposed by testing a clip plane
// that symmetrically bisects a polygon through one of its own
// vertices (see spherical_polygon test: "diagonal cut halves the
// octant"), not a hypothetical. Left undetected, the duplicate
// silently corrupts measure() (Girard's theorem is evaluated over the
// wrong vertex count).
template<Scalar T>
std::vector<Vec<T, 3>> dedup_consecutive(std::vector<Vec<T, 3>> pts, T radius) {
    if (pts.size() < 2) return pts;
    auto tol = epsilon<T>() * radius * radius;
    std::vector<Vec<T, 3>> out;
    out.reserve(pts.size());
    for (auto& p : pts) {
        if (out.empty() || Vec<T, 3>{p - out.back()}.norm_squared() > tol)
            out.push_back(p);
    }
    if (out.size() > 1 && Vec<T, 3>{out.front() - out.back()}.norm_squared() <= tol)
        out.pop_back();
    return out;
}

} // namespace detail

// ── Single half-space clip ─────────────────────────────────────
//
// Clip `subject` to the half-space { p : dot(p, normal) >= 0 }, a
// plane through the sphere's center -- i.e. one great-circle cut.
// This is the sphere analogue of one edge-iteration of planar
// Sutherland-Hodgman (detail::sutherland_hodgman in boolean.hpp):
// same "is this vertex on the inside?" sign test, generalized only in
// that the crossing point (detail::great_circle_edge_crossing) has
// two antipodal candidates instead of one. Public on its own merit --
// clipping an arbitrary spherical region against one coordinate
// hemisphere or bounding half-space is a useful primitive by itself,
// not just an spherical_intersection building block. An empty return
// means "clipped away entirely", not an error -- mirrors
// detail::sutherland_hodgman returning an empty vector for the same
// reason.
template<Scalar T>
SphericalPolygon<T> clip_half_space(const SphericalPolygon<T>& subject, const Vec<T, 3>& normal) {
    auto& in = subject.vertices;
    if (in.empty()) return {};

    auto tol = epsilon<T>() * subject.radius * subject.radius * normal.norm();
    auto is_inside = [&](const Vec<T, 3>& p) { return p.dot(normal) >= -tol; };

    std::vector<Vec<T, 3>> out;
    out.reserve(in.size() + 1);
    auto n = in.size();
    for (std::size_t j = 0; j < n; ++j) {
        auto& curr = in[j];
        auto& prev = in[(j + n - 1) % n];
        bool curr_in = is_inside(curr);
        bool prev_in = is_inside(prev);

        if (curr_in) {
            if (!prev_in)
                out.push_back(detail::great_circle_edge_crossing(prev, curr, normal, subject.radius));
            out.push_back(curr);
        } else if (prev_in) {
            out.push_back(detail::great_circle_edge_crossing(prev, curr, normal, subject.radius));
        }
    }
    out = detail::dedup_consecutive(std::move(out), subject.radius);
    return SphericalPolygon<T>{std::move(out), subject.radius};
}

// ── Boolean operations (convex spherical polygons only) ────────

// Intersection: clip `a` by each of `b`'s edge half-space planes in
// sequence -- exactly boolean.hpp's own multi-edge Sutherland-Hodgman
// composition, with clip_half_space standing in for one edge's clip.
template<Scalar T>
Result<SphericalPolygon<T>> spherical_intersection(const SphericalPolygon<T>& a,
                                                     const SphericalPolygon<T>& b) {
    if (a.vertices.size() < 3 || b.vertices.size() < 3)
        return std::unexpected(Error{ErrorCode::DegenerateInput, "need 3+ vertices"});

    auto result = a;
    auto bn = b.vertices.size();
    for (std::size_t i = 0; i < bn && !result.vertices.empty(); ++i) {
        auto e = b.edge(i);
        auto normal = Vec<T, 3>{e.a.cross(e.b)};
        result = clip_half_space(result, normal);
    }

    if (result.vertices.size() < 3)
        return std::unexpected(Error{ErrorCode::NoIntersection, "no overlap region"});
    return result;
}

// ADL hook so the generic `operator|` in geometry/make.hpp (a | b ==
// intersect(a, b), the project's "|" = intersect/pipe convention)
// picks up spherical polygons for free.
template<Scalar T>
Result<SphericalPolygon<T>> intersect(const SphericalPolygon<T>& a, const SphericalPolygon<T>& b) {
    return spherical_intersection(a, b);
}

// Difference: A \ B. Walks A's edges, keeping vertices outside B and
// inserting A-edge/B-edge great-circle crossing points in arc order --
// the direct generalization of boolean.hpp's difference_region.
// Like that function, this is not a claim of full generality: it
// assumes the removed portion leaves a single simple polygon behind
// (true whenever B clips a single convex "bite" out of A, which is
// the only shape of overlap two convex polygons can produce), not a
// robust arbitrary-topology boolean solver.
//
// CONFIRMED LIMITATION: exact vertex-on-boundary degeneracies where
// *several* coincide on the same handful of edges at once are not
// reliably handled. on_minor_arc's endpoints_inclusive flag (see its
// own comment) resolves the single-coincidence case -- a lone A
// vertex sitting exactly on one B edge, or vice versa -- by requiring
// a genuine crossing to be strictly interior to the edge it was found
// on while still accepting the other edge's own endpoint. That stops
// being sufficient once *multiple* coincidences land together (found
// by hand-tracing two octants related by an exact 45-degree rotation,
// which shares one vertex outright and lands a second exactly at
// another edge's midpoint): the strict/inclusive split can end up
// excluding a vertex that legitimately belongs on the A\B boundary,
// dropping the difference region's own apex. Two convex polygons in
// general position (the common case, and the one exercised in
// tests/test_spherical_polygon.cpp) do not hit this; polygons
// deliberately constructed to share multiple exact boundary points
// might.
template<Scalar T>
Result<SphericalPolygon<T>> spherical_difference(const SphericalPolygon<T>& a,
                                                   const SphericalPolygon<T>& b) {
    if (a.vertices.size() < 3 || b.vertices.size() < 3)
        return std::unexpected(Error{ErrorCode::DegenerateInput, "need 3+ vertices"});

    auto inter = spherical_intersection(a, b);
    if (!inter) return a; // no overlap -- difference is all of A

    if (std::abs(inter->measure() - a.measure()) < epsilon<T>() * a.measure() + epsilon<T>())
        return std::unexpected(Error{ErrorCode::NoIntersection, "A is entirely inside B"});

    std::vector<Vec<T, 3>> result;
    auto an = a.vertices.size();
    auto bn = b.vertices.size();

    for (std::size_t i = 0; i < an; ++i) {
        auto& p1 = a.vertices[i];
        auto& p2 = a.vertices[(i + 1) % an];
        if (!b.contains(p1)) result.push_back(p1);

        auto edge_normal_a = Vec<T, 3>{p1.cross(p2)};
        std::vector<std::pair<T, Vec<T, 3>>> hits;
        for (std::size_t j = 0; j < bn; ++j) {
            auto& q1 = b.vertices[j];
            auto& q2 = b.vertices[(j + 1) % bn];
            auto edge_normal_b = Vec<T, 3>{q1.cross(q2)};
            auto line_dir = Vec<T, 3>{edge_normal_a.cross(edge_normal_b)};
            auto len = line_dir.norm();
            if (len < epsilon<T>() * a.radius * a.radius) continue; // parallel great circles

            auto candidate = Vec<T, 3>{line_dir * (a.radius / len)};
            for (auto& cand : {candidate, Vec<T, 3>{-candidate}}) {
                // Strict on A's own edge (a coincidence with p1/p2
                // isn't a real crossing -- see on_minor_arc's own
                // comment), inclusive on B's edge (landing exactly on
                // one of B's vertices, q1/q2, is a perfectly real
                // crossing).
                if (detail::on_minor_arc(p1, p2, cand, a.radius, false) &&
                    detail::on_minor_arc(q1, q2, cand, a.radius, true)) {
                    using std::acos; using std::clamp;
                    auto cos_t = clamp(p1.dot(cand) / (a.radius * a.radius), T{-1}, T{1});
                    hits.emplace_back(acos(cos_t), cand);
                }
            }
        }
        std::sort(hits.begin(), hits.end(), [](auto& x, auto& y) { return x.first < y.first; });
        for (auto& [t, pt] : hits) result.push_back(pt);
    }

    result = detail::dedup_consecutive(std::move(result), a.radius);
    if (result.size() < 3)
        return std::unexpected(Error{ErrorCode::NoIntersection, "difference is degenerate"});
    return SphericalPolygon<T>{std::move(result), a.radius};
}

// Union, via gnomonic projection reusing the existing planar
// convex_hull(): a spherically-convex point set contained in an open
// hemisphere maps to an ordinary planar-convex point set under
// gnomonic projection (it sends great circles to straight lines), so
// the spherical convex hull of the combined vertices is exactly the
// planar convex hull of their gnomonic images, mapped back -- no
// bespoke spherical convex-hull algorithm needed. Like Polygon's own
// operator+ (also a convex hull of the combined vertices, not the
// literal union), this returns the tightest convex spherical polygon
// containing both inputs, not the (generally non-convex) union region
// itself -- see spherical_union_area() below for the exact union
// *area* without that approximation.
//
// RESTRICTION: fails (Result error) unless every vertex of both
// inputs lies strictly within one open hemisphere, checked against
// the combined centroid direction -- gnomonic projection has no
// finite image outside that hemisphere, so this fails loudly instead
// of silently returning a wrapped-around, wrong hull.
template<Scalar T>
Result<SphericalPolygon<T>> spherical_union(const SphericalPolygon<T>& a,
                                             const SphericalPolygon<T>& b) {
    if (a.vertices.size() < 3 || b.vertices.size() < 3)
        return std::unexpected(Error{ErrorCode::DegenerateInput, "need 3+ vertices"});
    if (std::abs(a.radius - b.radius) > epsilon<T>() * std::max(a.radius, b.radius))
        return std::unexpected(Error{ErrorCode::InvalidArgument, "radii must match"});

    std::vector<Vec<T, 3>> combined = a.vertices;
    combined.insert(combined.end(), b.vertices.begin(), b.vertices.end());

    Vec<T, 3> pole_sum{};
    for (auto& v : combined) pole_sum = Vec<T, 3>{pole_sum + v};
    auto pole_norm = pole_sum.norm();
    if (pole_norm < epsilon<T>())
        return std::unexpected(Error{ErrorCode::DegenerateInput,
            "combined vertices have no well-defined common hemisphere direction"});
    auto pole = Vec<T, 3>{pole_sum * (T{1} / pole_norm)};

    for (auto& v : combined) {
        auto u = Vec<T, 3>{v * (T{1} / a.radius)};
        if (u.dot(pole) <= epsilon<T>())
            return std::unexpected(Error{ErrorCode::OutOfDomain,
                "spherical_union requires both polygons within a single open hemisphere"});
    }

    auto [e1, e2] = detail::tangent_basis(pole);
    std::vector<Vec<T, 2>> planar;
    planar.reserve(combined.size());
    for (auto& v : combined)
        planar.push_back(detail::gnomonic_project(v, pole, e1, e2, a.radius));

    auto hull = convex_hull(std::move(planar));
    if (!hull) return std::unexpected(hull.error());

    std::vector<Vec<T, 3>> sphere_verts;
    sphere_verts.reserve(hull->vertices.size());
    for (auto& q : hull->vertices)
        sphere_verts.push_back(detail::gnomonic_unproject(q, pole, e1, e2, a.radius));

    return SphericalPolygon<T>{std::move(sphere_verts), a.radius};
}

// Exact union area via inclusion-exclusion -- area(A) + area(B) -
// area(A ∩ B) -- sidesteps spherical_union's convex-hull
// approximation and its hemisphere restriction entirely, mirroring
// boolean.hpp's own difference_area()/symmetric_difference_area(),
// which are also area-only precisely because the exact region is the
// hard part.
template<Scalar T>
Result<T> spherical_union_area(const SphericalPolygon<T>& a, const SphericalPolygon<T>& b) {
    auto inter = spherical_intersection(a, b);
    T inter_area = inter ? inter->measure() : T{0};
    return a.measure() + b.measure() - inter_area;
}

// The curved analogue of Polygon::area(): a free-function alias
// matching the naming of spherical_intersection()/spherical_union()/
// spherical_difference() above. Forwards to measure(), never
// reimplements, per docs/conventions.md's measure()-alias rule.
template<Scalar T>
T spherical_area(const SphericalPolygon<T>& poly) { return poly.measure(); }

// ── Operators (convex spherical polygons only) ──────────────────
// Matches the project-wide convention: | = intersect/pipe (see the
// intersect() ADL hook above), & = boolean intersect, + = union,
// - = difference.

template<Scalar T>
Result<SphericalPolygon<T>> operator&(const SphericalPolygon<T>& a, const SphericalPolygon<T>& b) {
    return spherical_intersection(a, b);
}

template<Scalar T>
Result<SphericalPolygon<T>> operator-(const SphericalPolygon<T>& a, const SphericalPolygon<T>& b) {
    return spherical_difference(a, b);
}

template<Scalar T>
Result<SphericalPolygon<T>> operator+(const SphericalPolygon<T>& a, const SphericalPolygon<T>& b) {
    return spherical_union(a, b);
}

} // namespace spatium::geometry

// std::formatter -- CONTRIBUTING.md requires one for any new public
// type with a natural textual representation.
template<spatium::Scalar T>
struct std::formatter<spatium::geometry::SphericalPolygon<T>> {
    constexpr auto parse(auto& ctx) { return ctx.begin(); }
    auto format(const spatium::geometry::SphericalPolygon<T>& p, auto& ctx) const {
        auto out = std::format_to(ctx.out(), "SphericalPolygon(r={}){{", p.radius);
        for (std::size_t i = 0; i < p.vertices.size(); ++i) {
            if (i > 0) out = std::format_to(out, ", ");
            out = std::format_to(out, "{}", p.vertices[i]);
        }
        return std::format_to(out, "}}");
    }
};
