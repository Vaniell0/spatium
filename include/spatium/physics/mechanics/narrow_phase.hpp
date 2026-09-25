#pragma once

// Analytical narrow-phase contact queries.
//
// Wraps the existing `geometry/ray_*` family (Quadric, Torus,
// ParametricSurface) as point-to-surface signed-distance queries,
// the form the IPC barrier (`contact.hpp`) needs:
//
//     d  = unsigned distance from query point to the closest
//          point on the surface
//     n  = outward unit normal at that closest point
//     cp = the closest surface point itself
//     in = true when the query point is on the "inside" of a
//          closed shape (signed distance would be negative)
//
// For analytical sphere / torus the math is closed form. For a
// general `ParametricSurface` we re-use the existing `find_params`
// (grid seed + Newton on (p − S(u,v))·Sᵤ = 0,  ·Sᵥ = 0), which is
// exactly the closest-point optimisation in disguise. The IPC
// barrier is then evaluated on the unsigned distance and rotated
// onto the surface normal to produce a force on the query point.
//
// Sub-pixel contact for *smooth* surfaces — without ever
// triangulating them — is the unique selling point this slice
// unlocks: cloth that drapes over a Klein-bottle pole or a
// billiard ball that rolls along an exact ellipsoid both fall
// out of the same `point_to_*` + `ipc_contact_force` pair.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/geometry/ray_surface.hpp>
#  include <spatium/physics/mechanics/contact.hpp>
#  include <spatium/spaces/parametric.hpp>
#  include <spatium/spaces/sphere.hpp>
#  include <algorithm>
#  include <cmath>
#  include <cstddef>
#  include <functional>
#  include <limits>
#  include <queue>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::physics::mechanics {

// ── Result of a point-to-surface narrow-phase query ────────────
template<Scalar T = double>
struct ContactQuery {
    T distance;                 // ≥ 0: unsigned distance to closest_point
    Vec<T, 3> closest_point;    // foot of the perpendicular on the surface
    Vec<T, 3> normal;           // outward unit normal at closest_point
    bool inside;                // true ⟺ query point is on the inside

    // Signed distance: negative iff inside. Useful for SDF-style debug.
    T signed_distance() const { return inside ? -distance : distance; }
};

// ── Point ↔ sphere (closed form) ───────────────────────────────
// Closest point on the sphere of radius r centred at `center` is
// the projection of (p − center) onto the radial direction.
template<Scalar T>
ContactQuery<T> point_to_sphere(const Vec<T, 3>& p,
                                const Vec<T, 3>& center,
                                T radius)
{
    Vec<T, 3> rel = p - center;
    T len = rel.norm();
    if (len < epsilon<T>()) {
        // Query exactly at centre — pick a canonical normal.
        Vec<T, 3> n{T{1}, T{0}, T{0}};
        return {radius, center + n * radius, n, true};
    }
    Vec<T, 3> n = rel / len;
    Vec<T, 3> cp = center + n * radius;
    using std::abs;
    return {abs(len - radius), cp, n, len < radius};
}

// ── Point ↔ torus (closed form) ────────────────────────────────
// Project to the torus' local frame, then to the centreline circle:
// the foot is the closest point on the circle of radius R, the
// remaining offset (in the meridional plane) of length r gives the
// tube-surface point. Sign of (|offset| − r) tells inside vs out.
template<Scalar T>
ContactQuery<T> point_to_torus(const Vec<T, 3>& p,
                               const ::spatium::geometry::Torus<T>& torus)
{
    Vec<T, 3> w = torus.axis;
    Vec<T, 3> u, v;
    ::spatium::geometry::torus_basis(w, u, v);

    Vec<T, 3> delta = p - torus.center;
    T x = delta.dot(u);
    T y = delta.dot(v);
    // delta·w (axial component) is implicit in (p − centreline) below.

    T planar = std::sqrt(x * x + y * y);
    if (planar < epsilon<T>()) {
        // Point sits on the torus axis — meridional ring is the
        // closest set; pick u as a canonical direction.
        Vec<T, 3> ring_dir = u;
        Vec<T, 3> centre = torus.center + ring_dir * torus.major_radius;
        Vec<T, 3> offset = p - centre;
        T off_len = offset.norm();
        if (off_len < epsilon<T>()) {
            // Ill-defined; fall back to outward radial.
            return {torus.minor_radius,
                    centre + ring_dir * torus.minor_radius,
                    ring_dir,
                    true};
        }
        Vec<T, 3> n = offset / off_len;
        Vec<T, 3> cp = centre + n * torus.minor_radius;
        using std::abs;
        return {abs(off_len - torus.minor_radius),
                cp, n,
                off_len < torus.minor_radius};
    }

    // Closest point on the centreline circle (radius R, in u/v plane).
    Vec<T, 3> centreline = torus.center
                         + u * (torus.major_radius * x / planar)
                         + v * (torus.major_radius * y / planar);

    Vec<T, 3> offset = p - centreline;
    T off_len = offset.norm();
    if (off_len < epsilon<T>()) {
        // p sits exactly on the centreline circle — meridional
        // direction undetermined; pick outward planar.
        Vec<T, 3> n = (u * (x / planar) + v * (y / planar));
        Vec<T, 3> cp = centreline + n * torus.minor_radius;
        return {torus.minor_radius, cp, n, true};
    }
    Vec<T, 3> n = offset / off_len;
    Vec<T, 3> cp = centreline + n * torus.minor_radius;
    using std::abs;
    return {abs(off_len - torus.minor_radius),
            cp, n,
            off_len < torus.minor_radius};
}

// ── IPC barrier evaluated on a contact query ───────────────────
// Energy is the IPC log-barrier on the *unsigned* distance, with
// the configured stiffness κ baked in. Returns 0 once outside the
// activation band [0, d_hat).
template<Scalar T>
T ipc_contact_energy(const ContactQuery<T>& q, T d_hat,
                     T kappa = T{1})
{
    T d = q.distance;
    if (d >= d_hat) return T{0};
    return kappa * ipc_barrier(d, d_hat);
}

// Force on the query point produced by the IPC barrier. The
// gradient ∂B/∂d is negative inside the active band, so −∂B/∂d ≥ 0;
// multiplying by the *outward* surface normal pushes the query
// point away from the surface, which is the physical convention.
//
// For a query already on the inside (signed distance < 0) the
// barrier is mathematically infinite — we return a large but
// finite outward force instead, so the caller can still take a
// repair step without hitting NaN.
template<Scalar T>
Vec<T, 3> ipc_contact_force(const ContactQuery<T>& q, T d_hat,
                            T kappa = T{1})
{
    if (q.inside) {
        // Hard repair direction: push outward along the surface
        // normal with the energy scale at d → 0⁺.
        return q.normal * (kappa * ipc_default_stiffness(d_hat));
    }
    if (q.distance >= d_hat) return Vec<T, 3>{};
    T g = ipc_barrier_grad(q.distance, d_hat);   // ≤ 0 in active band
    return q.normal * (-kappa * g);
}

// ── Concept-driven `point_to` overload set ─────────────────────
// The three analytical paths share the ContactQuery shape but
// otherwise have nothing in common — different math, different
// Big-O, different precision. Concept-overloads let the compiler
// pick the right one statically: `point_to(p, sphere)` resolves
// to the 3.7-ns closed-form path, `point_to(p, torus)` to the
// 17.7-ns one, `point_to(p, parametric)` to the µs Newton solver,
// all without runtime dispatch (no vtables, no tag tests).
//
// Bullet/PhysX style would force one common GJK code path for all
// three at ~100-500 ns each. We trade the boilerplate of multiple
// overloads for an order-of-magnitude speed win on the analytical
// shapes — and any user-supplied Surface type can join the set by
// adding its own `point_to(p, MyShape)` overload (or by satisfying
// the generic Surface fallback below).

// Centred-at-origin Sphere<N, T> from `spaces/`. Spatium convention:
// `Sphere<N>` is the unit N-sphere in R^(N+1), centred at the
// origin, so the only free parameter is `radius`. We restrict the
// physics overload to N == 2 (the surface of a 3D ball — the only
// case where a Vec<T, 3> contact query is well posed); higher N
// would imply a query in R^(N+1) and is outside this header's
// signature. The `requires (N == 2)` constraint surfaces that
// limitation as a compile-time error rather than a silent template
// substitution failure.
template<std::size_t N, Scalar T>
    requires (N == 2)
ContactQuery<T> point_to(const Vec<T, 3>& p,
                         const ::spatium::Sphere<N, T>& s) {
    return point_to_sphere(p, Vec<T, 3>{}, s.radius);
}

// Free torus / parametric overloads — delegate to the
// type-dispatched paths above.
template<Scalar T>
ContactQuery<T> point_to(const Vec<T, 3>& p,
                         const ::spatium::geometry::Torus<T>& t) {
    return point_to_torus(p, t);
}

// Generic fallback: any Surface (in the `spaces` sense) can serve
// as a contact target via its public `project` + `normal`. Picks
// up `ParametricSurface<T>`, `ImplicitSurface<T>`, custom user
// types — anything satisfying the Surface concept — through the
// grid-seed + Newton closest-point logic the surface itself owns.
// Slower than the analytical paths but always available.
template<Surface S>
ContactQuery<typename S::ScalarType> point_to(
    const Vec<typename S::ScalarType, 3>& p, const S& surf)
{
    using T = typename S::ScalarType;
    Vec<T, 3> cp = surf.project(p);
    // Normal is queried at the closest-point `cp`, not at the query
    // point `p`: for a particle far from the surface those two
    // sample sites can disagree, and only the normal at the foot
    // of the perpendicular gives a meaningful inside-test.
    Vec<T, 3> n  = surf.normal(cp);
    Vec<T, 3> rel = p - cp;
    T d = rel.norm();
    bool inside = rel.dot(n) < T{0};
    return {d, cp, n, inside};
}

// ── A distance with a certified floor ──────────────────────────
//
// `point_to` above answers with the distance to the point its search
// found, which is a real point on the surface and so an UPPER bound on
// the distance. Continuous collision needs the other side: a step is safe
// only if it is shorter than a distance the surface provably cannot be
// closer than. Where a local search misses a feature -- a thin spike
// between the seeds of the grid `find_params` starts from -- the upper
// bound can be several times the truth, and a step sized by it goes
// through the feature. `distance_bound` returns both sides:
//
//   lower    no point of the surface is closer than this
//   nearest  a real surface point and its distance, an upper bound
//
// Exact for the shapes with a closed form, where the two coincide. For a
// shape without one, a floor needs something the shape alone does not
// say -- how fast it can move -- and the wrappers below take it as an
// explicit Lipschitz constant. Without one no finite number of samples
// can prove that nothing sits between them, so there is no certified
// overload for a bare `Surface`.
template<Scalar T>
struct DistanceBound {
    T lower;
    ContactQuery<T> nearest;
    std::size_t evaluations = 0;   // what it cost, in evaluations of the shape
};

template<std::size_t N, Scalar T>
    requires (N == 2)
DistanceBound<T> distance_bound(const Vec<T, 3>& p, const ::spatium::Sphere<N, T>& s, T = T{0},
                                std::size_t = 0) {
    const auto q = point_to(p, s);
    return {q.distance, q, 1};
}

template<Scalar T>
DistanceBound<T> distance_bound(const Vec<T, 3>& p, const ::spatium::geometry::Torus<T>& t, T = T{0},
                                std::size_t = 0) {
    const auto q = point_to(p, t);
    return {q.distance, q, 1};
}

// A chart together with a Lipschitz constant of its map:
// |f(a) - f(b)| <= lipschitz * |a - b| for parameters a, b in its domain,
// measured in the Euclidean norm of (u, v). The operator norm of the
// Jacobian, bounded over the domain, is one; for a height field z = h(u, v)
// it is sqrt(1 + max |grad h|^2), for a sphere of radius r in angles it is r.
//
// One constant for the whole domain is loose where a chart degenerates: at
// a sphere's pole a cell spanning every u has a tiny image, while L times
// its parameter size is large, and a search near the pole splits in u far
// past any use. `cell_radius`, when given, replaces L * rho with a caller's
// own bound on how far the image of the cell [u0, u1] x [v0, v1] reaches
// from the image of its centre -- for the sphere, r * max sin(v) * du / 2 +
// r * dv / 2. It is trusted exactly as the constant is.
template<Scalar T>
struct LipschitzChart {
    const ::spatium::ParametricSurface<T>& chart;
    T lipschitz;
    std::function<T(T u0, T u1, T v0, T v1)> cell_radius = {};
};

// Branch and bound over the parameter rectangle. A cell with centre c and
// half-diagonal rho holds no point closer than |p - f(c)| - L * rho, so the
// smallest such bound over a set of cells covering the domain is a floor
// for the whole surface. The cell holding the smallest bound is halved
// until the gap it proves above `floor` is at least half the gap of the
// best real distance found -- enough for an advancing step to make
// progress -- or that distance is already within `floor`, or the budget of
// evaluations is spent. `floor` is where the caller calls it contact: a
// sweep's tolerance plus the radius of the ball it moves. Every exit
// returns a valid lower bound; the budget only decides how tight.
template<Scalar T>
DistanceBound<T> distance_bound(const Vec<T, 3>& p, const LipschitzChart<T>& c, T floor = T{0},
                                std::size_t budget = std::size_t{1} << 18) {
    using std::sqrt;
    struct Cell { T u0, u1, v0, v1, lower; };
    const auto greater = [](const Cell& a, const Cell& b) { return a.lower > b.lower; };
    std::priority_queue<Cell, std::vector<Cell>, decltype(greater)> open(greater);

    const auto& dom = c.chart.domain();
    T best = std::numeric_limits<T>::infinity();
    T best_u = dom.u_min, best_v = dom.v_min;
    Vec<T, 3> best_point{};
    std::size_t evaluations = 0;

    // Rounding in |p - f(c)| is a few ulps of the coordinates involved;
    // the margin keeps the floor a floor in floating point too.
    const T ulp = std::numeric_limits<T>::epsilon() * T{8};
    const auto reach = [&c](T u0, T u1, T v0, T v1) {
        if (c.cell_radius) return c.cell_radius(u0, u1, v0, v1);
        const T du = u1 - u0, dv = v1 - v0;
        return c.lipschitz * sqrt(du * du + dv * dv) / T{2};
    };
    const auto push = [&](T u0, T u1, T v0, T v1) {
        const T u = (u0 + u1) / T{2}, v = (v0 + v1) / T{2};
        const Vec<T, 3> x = c.chart(u, v);
        ++evaluations;
        const T d = Vec<T, 3>{p - x}.norm();
        if (d < best) { best = d; best_u = u; best_v = v; best_point = x; }
        open.push(Cell{u0, u1, v0, v1, d - reach(u0, u1, v0, v1) - ulp * (p.norm() + x.norm())});
    };

    push(dom.u_min, dom.u_max, dom.v_min, dom.v_max);
    T lower = open.top().lower;
    while (true) {
        const Cell cell = open.top();
        lower = cell.lower;
        if (best <= floor || lower - floor >= (best - floor) / T{2} || evaluations + 2 > budget)
            break;
        open.pop();
        // Halve the side whose halving shrinks the reach more. Splitting
        // both at once wastes most of the work where a chart degenerates:
        // near a sphere's pole only v matters, and a four-way split would
        // double the cells in u at every level for nothing.
        const T um = (cell.u0 + cell.u1) / T{2}, vm = (cell.v0 + cell.v1) / T{2};
        if (reach(cell.u0, um, cell.v0, cell.v1) <= reach(cell.u0, cell.u1, cell.v0, vm)) {
            push(cell.u0, um, cell.v0, cell.v1);
            push(um, cell.u1, cell.v0, cell.v1);
        } else {
            push(cell.u0, cell.u1, cell.v0, vm);
            push(cell.u0, cell.u1, vm, cell.v1);
        }
    }
    Vec<T, 3> n = c.chart.normal_at(best_u, best_v);
    return {std::max(lower, T{0}), ContactQuery<T>{best, best_point, n, false}, evaluations};
}

// An implicit surface f = 0 with a Lipschitz constant of f: |f(p)| / L is
// a floor on the distance to the zero set, the bound a sphere tracer
// steps by. It names no surface point, so `nearest` carries the floor as
// its distance too and a contact is declared once that floor is small.
template<Scalar T>
struct LipschitzImplicit {
    std::function<T(const Vec<T, 3>&)> f;
    T lipschitz;
};

template<Scalar T>
DistanceBound<T> distance_bound(const Vec<T, 3>& p, const LipschitzImplicit<T>& s, T = T{0},
                                std::size_t = 0) {
    using std::abs;
    const T value = s.f(p);
    const T lower = abs(value) / s.lipschitz;
    return {lower, ContactQuery<T>{lower, p, Vec<T, 3>{}, value < T{0}}, 1};
}

// `floor` is where a search may stop because contact is already found;
// `budget` caps its evaluations, and a shape with a closed form ignores both.
template<typename S, typename T>
concept HasDistanceBound = requires(const S& s, const Vec<T, 3>& p, T floor, std::size_t budget) {
    { distance_bound(p, s, floor, budget) } -> std::same_as<DistanceBound<T>>;
};

// The `ContactSurface<S>` concept is satisfied by any type for
// which `point_to(p, s)` is callable and returns a ContactQuery.
// Add a `point_to` overload for your own type and it joins the set.
template<typename S, typename T = double>
concept ContactSurface = requires(const S& s, const Vec<T, 3>& p) {
    { point_to(p, s) } -> std::same_as<ContactQuery<T>>;
};

// Concept-constrained convenience: query + force in one call.
template<typename S, Scalar T>
    requires ContactSurface<S, T>
Vec<T, 3> ipc_contact_force_on(const Vec<T, 3>& p, const S& surf,
                               T d_hat, T kappa = T{1}) {
    return ipc_contact_force(point_to(p, surf), d_hat, kappa);
}

} // namespace spatium::physics::mechanics
