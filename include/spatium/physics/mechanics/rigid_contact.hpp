#pragma once

// Pairwise narrow-phase contact queries between two DYNAMIC rigid
// bodies, plus the O(n^2) broad-phase AABB sweep that feeds them.
//
// `narrow_phase.hpp` answers "how far is this point from a STATIC
// analytical surface" (a cloth vertex against a fixed sphere/torus).
// That framing doesn't fit rigid-rigid collision: neither sphere in a
// sphere-sphere pair, and neither the sphere nor the box in a
// sphere-box pair, is a fixed shape that defines the coordinate frame
// the other lives in — both operands can move every substep. This
// header ships the pairwise closed-form queries that case needs:
//
//   sphere_sphere_contact(center_a, radius_a, center_b, radius_b)
//   sphere_aabb_contact(sphere_center, radius, box_min, box_max)
//   sweep_sphere_sphere(center_a, radius_a, disp_a, center_b, radius_b, disp_b)
//
// The first two are discrete (snapshot) queries — do these two shapes
// touch right now. `sweep_sphere_sphere` is the continuous counterpart
// for sphere-sphere: given each body's motion over the step, does the
// swept path cross the other body at all, closing the tunnelling gap a
// discrete-only check leaves open for fast movers (see that function's
// own comment for the reduction to `geometry::ray_quadric`'s existing
// closed-form root solve, and for exactly what is/isn't covered).
//
// The discrete queries return `ContactQuery<T>` (reused verbatim from narrow_phase.hpp
// — `distance`/`normal`/`inside` mean exactly what they mean there,
// `closest_point` is documented per-query below since "foot of the
// perpendicular from a single point" doesn't quite apply when both
// operands are extended bodies).
//
// Scope, deliberately: two spheres, and one sphere against one
// axis-aligned box. No box-box, no oriented boxes, no capsules — the
// project's documented cost/benefit call is that a full narrow-phase
// library (the union-of-convex-shapes generality ipc-toolkit or
// Bullet/PhysX provide) is out of scope; this is the small, closed-
// form slice that turns "particle vs. static surface" into "body vs.
// body" for the shapes actually needed by a falling-spheres demo.
//
// Broad phase: `broad_phase_aabb_pairs` is a plain O(n^2) all-pairs
// AABB overlap sweep over `geometry::Box<3, T>`. For the tens-of-
// bodies scale this header targets, that is the right amount of
// machinery — `spatial/bvh.hpp` is built and tuned for many static
// primitives under repeated single queries (ray casts, nearest-point),
// not for an all-pairs sweep over a handful of moving bodies that is
// rebuilt every substep; reaching for it here would trade a five-line
// double loop for a tree-build most of whose cost never pays back at
// this body count.
//
// `swept_sphere_aabb` is the broad-phase counterpart to
// `sweep_sphere_sphere` above: the union of a sphere's start- and
// end-of-step boxes, so a fast mover's swept path isn't missed by the
// broad phase either, one level up from where `sweep_sphere_sphere`
// closes the same tunnelling gap in the narrow phase.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <spatium/geometry/box.hpp>
#  include <spatium/geometry/ray_surface.hpp>
#  include <spatium/physics/mechanics/narrow_phase.hpp>
#  include <spatium/physics/mechanics/xpbd.hpp>
#  include <algorithm>
#  include <cmath>
#  include <cstddef>
#  include <limits>
#  include <utility>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::physics::mechanics {

// ── Sphere ↔ sphere (closed form) ───────────────────────────────
// `distance`: unsigned gap between the two surfaces along the line
// joining the centers — |d - (radius_a + radius_b)| where d is the
// center-to-center distance. Zero when the spheres are exactly
// touching, growing whether they separate further or overlap deeper
// (mirrors `point_to_sphere`'s |len - radius| convention).
// `normal`: unit vector from A's center toward B's center — order-
// dependent by construction (this is a pair query, not a surface
// query with a fixed "outward" side); the caller applies +normal to
// push A away from B and -normal to push B away from A, exactly the
// convention `xpbd_solve_sphere_collision` below assumes.
// `closest_point`: midpoint of the two spheres' surface points along
// the center line — the natural, order-independent contact location.
// It coincides with both surface points exactly when touching and is
// still well defined (approximate, not "the" contact point) when
// overlapping — same spirit as narrow_phase.hpp's own "closest point"
// naming, generalized from one static surface to two moving ones.
// `inside`: true when the spheres overlap (this is not a containment
// test — it reuses the field name for "penetrating," matching
// `point_to_sphere`'s use of `inside` for "closer than the surface").
template<Scalar T>
ContactQuery<T> sphere_sphere_contact(const Vec<T, 3>& center_a, T radius_a,
                                      const Vec<T, 3>& center_b, T radius_b)
{
    Vec<T, 3> rel = center_b - center_a;
    T d = rel.norm();
    T rest = radius_a + radius_b;

    Vec<T, 3> n;
    if (d < epsilon<T>()) {
        // Coincident centers: direction is undefined, pick a canonical
        // axis (same fallback `point_to_sphere` uses at its own
        // degenerate point).
        n = Vec<T, 3>{T{1}, T{0}, T{0}};
    } else {
        n = rel / d;
    }

    Vec<T, 3> surf_a = center_a + n * radius_a;
    Vec<T, 3> surf_b = center_b - n * radius_b;
    Vec<T, 3> cp = Vec<T, 3>{(surf_a + surf_b) * T{0.5}};

    using std::abs;
    return {abs(d - rest), cp, n, d < rest};
}

// ── Swept (continuous) sphere ↔ sphere contact ──────────────────
// `sphere_sphere_contact` above is a snapshot: exactly the wrong
// question for a body moving fast enough to cross another between one
// substep's start and end position while being separated at both.
// Relative to A, the question "does B's swept motion cross the sphere
// of radius (radius_a + radius_b) centred on A" is *exactly* the
// ray-vs-sphere problem `geometry::ray_quadric` already solves in
// closed form (used throughout the ray-tracing examples) — the "ray"
// here is just the relative position/displacement instead of a camera
// ray. Reusing it means the time-of-impact comes from the same tested
// quadratic root solve the renderer trusts, not a second hand-rolled
// one.
//
// `disp_a`/`disp_b` are each body's full displacement over the step
// under test (velocity * dt, or simply x_end - x_start) — not a
// per-unit-time velocity, so there is no separate `dt` parameter, and
// `toi` comes back as a fraction of that same displacement, in [0, 1].
//
// Scope, deliberately: exact for sphere-sphere (and, as the zero-
// velocity special case, sphere vs. a STATIC sphere obstacle). NOT
// extended to sphere-vs-general-quadric (ellipsoid/cylinder/cone) or
// sphere-vs-torus: the Minkowski sum of a ball with a non-spherical
// quadric is not itself a quadric (offsetting a curved surface by a
// constant radius raises its algebraic degree), so those cases don't
// reduce to one `ray_quadric` call the way sphere-sphere does — an
// exact treatment would need a quartic (torus) or higher-degree solve
// this change does not attempt; deferred rather than shipped as a
// silently-approximate "close enough" version. Sphere-AABB continuous
// detection is deferred for the same shape of reason — a box grown by
// a ball is a *rounded* box, not a plain one — so `sphere_aabb_contact`
// above stays a discrete, per-substep query, with the same tunnelling
// exposure a fast enough body would have against any discrete check.
template<Scalar T>
struct SweptContact {
    bool hit;                // true iff a first-touch time exists in [0, 1]
    T toi;                   // fraction of the step at first contact (0 if
                             // already touching at the start of the step;
                             // 1, with hit=false, if they never touch)
    ContactQuery<T> contact; // geometry at time `toi` (reuses sphere_sphere_contact)
};

template<Scalar T>
SweptContact<T> sweep_sphere_sphere(const Vec<T, 3>& center_a0, T radius_a,
                                    const Vec<T, 3>& disp_a,
                                    const Vec<T, 3>& center_b0, T radius_b,
                                    const Vec<T, 3>& disp_b)
{
    // Already touching/penetrating at the start of the step: report
    // immediate contact rather than solving the sweep. Not just an
    // optimization — if the pair starts overlapping, the relative
    // "ray" below starts INSIDE the target sphere, and the one
    // non-negative root `ray_quadric` would return is the *exit*
    // point (where they'd separate back to exactly touching), not
    // "touching right now," so running the sweep in this branch would
    // give the wrong answer, not just a slower one.
    auto q0 = sphere_sphere_contact(center_a0, radius_a, center_b0, radius_b);
    if (q0.inside || q0.distance <= epsilon<T>())
        return {true, T{0}, q0};

    Vec<T, 3> rel0 = center_b0 - center_a0;
    Vec<T, 3> rel_disp = disp_b - disp_a;
    T rest = radius_a + radius_b;

    // Degenerate: no relative motion this step. Already known
    // separated (the branch above would have caught touching or
    // overlapping), so they stay separated — this also sidesteps
    // handing `ray_quadric`/`solve_quadratic` a zero direction, which
    // this header doesn't rely on being safe to divide by.
    if (rel_disp.norm() < epsilon<T>()) {
        Vec<T, 3> a1 = Vec<T, 3>{center_a0 + disp_a};
        Vec<T, 3> b1 = Vec<T, 3>{center_b0 + disp_b};
        return {false, T{1}, sphere_sphere_contact(a1, radius_a, b1, radius_b)};
    }

    auto hits = ::spatium::geometry::ray_quadric(
        ::spatium::geometry::Ray<3, T>{rel0, rel_disp},
        ::spatium::geometry::Quadric<T>::sphere(rest));

    for (auto& h : hits) {
        if (h.t >= T{0} && h.t <= T{1}) {
            T toi = h.t;
            Vec<T, 3> pa = Vec<T, 3>{center_a0 + disp_a * toi};
            Vec<T, 3> pb = Vec<T, 3>{center_b0 + disp_b * toi};
            return {true, toi, sphere_sphere_contact(pa, radius_a, pb, radius_b)};
        }
    }

    Vec<T, 3> a1 = Vec<T, 3>{center_a0 + disp_a};
    Vec<T, 3> b1 = Vec<T, 3>{center_b0 + disp_b};
    return {false, T{1}, sphere_sphere_contact(a1, radius_a, b1, radius_b)};
}

// ── Sphere ↔ axis-aligned box (closed form) ─────────────────────
// Shallow case (sphere center outside the box): the closest point on
// the box is `Box::project`'s clamp, and the sphere-to-box gap is the
// distance from the center to that point minus the radius — same
// shape as `sphere_sphere_contact` above with the box's "radius" being
// zero at the clamped point.
// Deep case (sphere center at or inside the box, `distance == 0` from
// the clamp): the clamp alone can't tell which face is "the" contact,
// so this falls back to the minimum-translation axis — the face whose
// plane is nearest the center — the standard AABB deep-penetration
// rule (Ericson, "Real-Time Collision Detection" §5.2.5). Only matters
// for a sphere that has tunnelled fully inside the box in one substep;
// a per-substep XPBD/broad-phase loop with reasonable step sizes
// should rarely hit it, but the fallback keeps the query total instead
// of returning a zero/undefined normal.
template<Scalar T>
ContactQuery<T> sphere_aabb_contact(const Vec<T, 3>& sphere_center, T radius,
                                    const Vec<T, 3>& box_min, const Vec<T, 3>& box_max)
{
    using std::clamp;
    Vec<T, 3> closest{clamp(sphere_center[0], box_min[0], box_max[0]),
                       clamp(sphere_center[1], box_min[1], box_max[1]),
                       clamp(sphere_center[2], box_min[2], box_max[2])};
    Vec<T, 3> delta = sphere_center - closest;
    T dist = delta.norm();

    using std::abs;
    if (dist > epsilon<T>()) {
        Vec<T, 3> n = delta / dist;
        T gap = dist - radius;
        return {abs(gap), closest, n, gap < T{0}};
    }

    // Deep case: center is inside (or exactly on the boundary of) the
    // box. Find the nearest face by scanning the six candidate
    // penetration depths along each axis.
    T best = std::numeric_limits<T>::max();
    Vec<T, 3> n{T{0}, T{0}, T{1}};
    Vec<T, 3> cp = sphere_center;
    for (std::size_t i = 0; i < 3; ++i) {
        T to_max = box_max[i] - sphere_center[i];
        T to_min = sphere_center[i] - box_min[i];
        if (to_max < best) {
            best = to_max;
            n = Vec<T, 3>{};
            n[i] = T{1};
            cp = sphere_center;
            cp[i] = box_max[i];
        }
        if (to_min < best) {
            best = to_min;
            n = Vec<T, 3>{};
            n[i] = T{-1};
            cp = sphere_center;
            cp[i] = box_min[i];
        }
    }
    T gap = -(best + radius);
    return {abs(gap), cp, n, true};
}

// ── XPBD collision constraint (sphere ↔ sphere) ─────────────────
// Same compliance/Lagrange-multiplier machinery as
// `XpbdDistanceConstraint`/`xpbd_solve_distance` in xpbd.hpp, applied
// to a unilateral (inequality) constraint instead of an equality one:
// C(x_i, x_j) = |x_i - x_j| - (radius_i + radius_j) must stay >= 0,
// so the projection only fires while C < 0 (penetrating) and is a
// no-op otherwise — a non-penetration constraint may push bodies
// apart, never pull them together, which is exactly what skipping the
// projection at C >= 0 encodes. `radius_i`/`radius_j` live on the
// constraint (not on `XpbdParticle`, which stays shape-agnostic) so
// the same particle can take part in constraints with different
// effective radii if a caller ever needs that.
template<std::size_t N, Scalar T = double>
struct XpbdSphereCollisionConstraint {
    std::size_t i;
    std::size_t j;
    T radius_i;
    T radius_j;
    T compliance{T{0}};    // 0 = perfectly rigid contact (hard PBD limit)
    T lambda{T{0}};        // accumulated Lagrange multiplier

    void reset() { lambda = T{0}; }
};

// Project one sphere-sphere collision constraint by Gauss-Seidel —
// identical math to `xpbd_solve_distance`, gated by C < 0.
template<std::size_t N, Scalar T = double>
void xpbd_solve_sphere_collision(XpbdSphereCollisionConstraint<N, T>& c,
                                 std::vector<XpbdParticle<N, T>>& parts,
                                 T dt)
{
    auto& pi = parts[c.i];
    auto& pj = parts[c.j];

    Vec<T, N> diff = Vec<T, N>{pi.x - pj.x};
    using std::sqrt;
    T len = sqrt(diff.dot(diff));
    if (len < epsilon<T>()) return;

    T rest = c.radius_i + c.radius_j;
    T C = len - rest;
    if (C >= T{0}) return;   // not penetrating — inequality constraint inactive

    Vec<T, N> grad = Vec<T, N>{diff * (T{1} / len)};
    T w_sum = pi.w + pj.w;
    if (w_sum < epsilon<T>()) return;                        // both pinned

    T alpha_tilde = c.compliance / (dt * dt);
    T denom = w_sum + alpha_tilde;
    T dlambda = -(C + alpha_tilde * c.lambda) / denom;
    c.lambda += dlambda;

    Vec<T, N> correction = Vec<T, N>{grad * dlambda};
    pi.x = Vec<T, N>{pi.x + correction * pi.w};
    pj.x = Vec<T, N>{pj.x - correction * pj.w};
}

// ── Broad phase: O(n^2) AABB overlap sweep ──────────────────────
// See file header for why a full spatial structure is the wrong tool
// at this body count. `sphere_aabb` is the one-line helper every
// caller needs to turn a moving sphere into the `Box<3, T>` this sweep
// (and `Box::intersects`) consume.
template<Scalar T>
::spatium::geometry::Box<3, T> sphere_aabb(const Vec<T, 3>& center, T radius) {
    return ::spatium::geometry::Box<3, T>::from_center_half_extents(
        center, Vec<T, 3>{radius, radius, radius});
}

// Swept AABB for one sphere over a step: the union of its start- and
// end-of-step bounding boxes. Feed this to `broad_phase_aabb_pairs`
// instead of `sphere_aabb` at just the end position when bodies can
// move fast enough that their end-of-step boxes might not even
// overlap despite crossing paths mid-step — the same tunnelling gap
// `sweep_sphere_sphere` closes in the narrow phase, one level up in
// the broad phase.
template<Scalar T>
::spatium::geometry::Box<3, T> swept_sphere_aabb(const Vec<T, 3>& center0, T radius,
                                                 const Vec<T, 3>& disp) {
    auto start = sphere_aabb(center0, radius);
    auto end = sphere_aabb(Vec<T, 3>{center0 + disp}, radius);
    return start.union_with(end);
}

template<Scalar T>
std::vector<std::pair<std::size_t, std::size_t>>
broad_phase_aabb_pairs(const std::vector<::spatium::geometry::Box<3, T>>& aabbs)
{
    std::vector<std::pair<std::size_t, std::size_t>> pairs;
    std::size_t n = aabbs.size();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i + 1; j < n; ++j)
            if (aabbs[i].intersects(aabbs[j]))
                pairs.emplace_back(i, j);
    return pairs;
}

} // namespace spatium::physics::mechanics
