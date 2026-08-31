#pragma once

// Extended Position-Based Dynamics (XPBD) for cloth and rods.
//
// Macklin, Müller, Chentanez 2016. "XPBD: Position-Based Simulation of
// Compliant Constrained Dynamics."
//
// XPBD treats stiffness as a Lagrange multiplier rather than a spring
// constant, decoupling time-step size from material rigidity. Each
// constraint C(x) = 0 contributes a corrective displacement that's
// proportional to its violation but inversely proportional to its
// "compliance" α = 1/k. For a stiff edge (k → ∞, α → 0) the constraint
// is enforced exactly per substep regardless of dt; for a soft edge
// (small k) the response decays gracefully.
//
// What's here:
//   - `XpbdParticle` — point with position, prediction, and inverse mass.
//   - `XpbdDistanceConstraint` — keeps two particles at a target distance
//     (cloth edge, rope segment).
//   - `xpbd_solve_distance` — single Gauss-Seidel pass projecting one
//     constraint, applies the lagrangian correction.
//   - `xpbd_step` — one substep: predict, iterate constraints `n_iter`
//     times, update velocities from corrected positions.
//
// Bending and other constraints follow the same template; this slice
// ships the canonical distance constraint and the distance-bending
// helper. Cloth-on-obstacle demos live in `examples/` and are slated
// to return on top of the implicit IPC backend (see
// `examples/CMakeLists.txt` for the current state).

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <spatium/core/concepts.hpp>
#  include <spatium/core/epsilon.hpp>
#  include <array>
#  include <cstddef>
#  include <cstdint>
#  include <unordered_map>
#  include <unordered_set>
#  include <utility>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::physics::mechanics {

template<std::size_t N, Scalar T = double>
struct XpbdParticle {
    Vec<T, N> x{};         // current position
    Vec<T, N> x_prev{};    // previous step (used to recover velocity)
    T         w{T{1}};     // inverse mass (0 ⇒ pinned)
};

// Distance constraint between particles `i` and `j`:
//   C(x_i, x_j) = |x_i − x_j| − rest
// XPBD compliance α (= 1/k) is per-edge; lagrangian λ accumulates
// across iterations within one timestep.
template<std::size_t N, Scalar T = double>
struct XpbdDistanceConstraint {
    std::size_t i;
    std::size_t j;
    T rest;
    T compliance{T{0}};    // 0 = perfectly stiff (PBD limit)
    T lambda{T{0}};        // accumulated Lagrange multiplier

    void reset() { lambda = T{0}; }
};

// Project one distance constraint by Gauss-Seidel.
//
//   Δλ = -(C + α̃·λ) / ((w_i + w_j) + α̃)
//   Δx_i = +Δλ · w_i · ∇C_i
//   Δx_j = +Δλ · w_j · ∇C_j
//
// where α̃ = α / dt² is the time-step-rescaled compliance. With α = 0 this
// reduces to PBD's hard projection.
template<std::size_t N, Scalar T = double>
void xpbd_solve_distance(XpbdDistanceConstraint<N, T>& c,
                         std::vector<XpbdParticle<N, T>>& parts,
                         T dt)
{
    auto& pi = parts[c.i];
    auto& pj = parts[c.j];

    Vec<T, N> diff = Vec<T, N>{pi.x - pj.x};
    using std::sqrt;
    T len = sqrt(diff.dot(diff));
    if (len < epsilon<T>()) return;
    T C = len - c.rest;

    Vec<T, N> grad = Vec<T, N>{diff * (T{1} / len)};         // ∇C w.r.t. x_i
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

// ── Distance-based bending constraint ─────────────────────────
// For each *interior* edge of a triangle mesh (an edge shared by
// exactly two triangles), constrain the distance between the two
// vertices that lie *opposite* the edge in the two triangles. A
// flat sheet has a fixed rest-distance; folding the sheet brings
// the opposite vertices closer, so the constraint resists folding.
//
// This is the "distance bending" model used by Houdini Vellum's
// `Stiffness Mode = Distance` and by many XPBD cloth implementations
// as a simpler, singularity-free alternative to the dihedral-angle
// formulation (Müller 2007 §4.5). It is mechanically just an
// XpbdDistanceConstraint between non-adjacent vertices, but the
// helper `build_bending_distance_constraints` makes its purpose
// explicit and isolates it from the structural distance edges so
// callers can tune compliance separately ("how stiff is the cloth's
// resistance to bending?" vs "how stretchy are the threads?").
//
// We reuse XpbdDistanceConstraint and xpbd_solve_distance verbatim
// — the type alias documents intent, the build helper does the
// edge → opposite-vertex pairing.
template<std::size_t N, Scalar T = double>
using XpbdBendingDistanceConstraint = XpbdDistanceConstraint<N, T>;

// One XPBD substep:
//   1. Reset Lagrange multipliers.
//   2. Predict positions from current velocities + external accelerations.
//   3. Iterate constraints `n_iter` times (Gauss-Seidel).
//   4. Recover velocities as (x − x_prev) / dt.
//
// External per-particle acceleration is supplied through `external_accel`,
// a callable `(particle, dt) → Vec<T, N>` (typically gravity). The name
// reflects what the callable returns — the *acceleration*, not a side
// effect — so the loop body below stays self-explanatory.
template<std::size_t N, Scalar T, typename ExternalAccel, typename ConstraintIter>
void xpbd_step(std::vector<XpbdParticle<N, T>>& parts,
               ConstraintIter constraints_begin,
               ConstraintIter constraints_end,
               ExternalAccel&& external_accel,
               T dt, int n_iter = 4)
{
    // Predict.
    for (auto& p : parts) {
        if (p.w <= T{0}) { p.x_prev = p.x; continue; }     // pinned
        Vec<T, N> a = external_accel(p, dt);
        Vec<T, N> v = Vec<T, N>{(p.x - p.x_prev) / dt};
        Vec<T, N> v_pred = Vec<T, N>{v + a * dt};
        p.x_prev = p.x;
        p.x = Vec<T, N>{p.x + v_pred * dt};
    }

    // Reset λ for this step.
    for (auto it = constraints_begin; it != constraints_end; ++it)
        it->reset();

    // Iterative projection.
    for (int iter = 0; iter < n_iter; ++iter) {
        for (auto it = constraints_begin; it != constraints_end; ++it)
            xpbd_solve_distance(*it, parts, dt);
    }

    // Velocities are implicit in (x - x_prev) / dt for the next step;
    // no explicit update needed.
}

// internal — do not use, no API stability. EdgeKey/EdgeKeyHash
// and edge_key() back the unique-edge enumeration in
// build_distance_constraints(); rely on that public helper, not
// on these.
namespace detail {

struct EdgeKey {
    std::uint32_t a, b;
    constexpr bool operator==(const EdgeKey&) const noexcept = default;
};

struct EdgeKeyHash {
    std::size_t operator()(const EdgeKey& k) const noexcept {
        return std::hash<std::uint64_t>{}(
            (static_cast<std::uint64_t>(k.a) << 32) | k.b);
    }
};

inline EdgeKey edge_key(std::uint32_t a, std::uint32_t b) {
    return (a < b) ? EdgeKey{a, b} : EdgeKey{b, a};
}

} // namespace detail

// Helper for typical cloth setup — build distance constraints for every
// edge of a triangle mesh (each pair of vertices sharing an edge gets
// one constraint, with rest length = current Euclidean distance).
template<std::size_t N, Scalar T = double>
std::vector<XpbdDistanceConstraint<N, T>>
build_distance_constraints(const std::vector<XpbdParticle<N, T>>& parts,
                           const std::vector<std::array<std::uint32_t, 3>>& faces,
                           T compliance = T{0})
{
    using std::sqrt;
    std::vector<XpbdDistanceConstraint<N, T>> out;
    out.reserve(faces.size() * 3 / 2);
    std::unordered_set<detail::EdgeKey, detail::EdgeKeyHash> seen;
    seen.reserve(faces.size() * 3);

    auto add_edge = [&](std::uint32_t a, std::uint32_t b) {
        auto k = detail::edge_key(a, b);
        if (!seen.insert(k).second) return;
        Vec<T, N> diff = Vec<T, N>{parts[k.a].x - parts[k.b].x};
        T rest = sqrt(diff.dot(diff));
        out.push_back({k.a, k.b, rest, compliance, T{0}});
    };

    for (auto& f : faces) {
        add_edge(f[0], f[1]);
        add_edge(f[1], f[2]);
        add_edge(f[2], f[0]);
    }
    return out;
}

// Build distance-bending constraints — one per interior edge
// (edge shared by exactly two triangles). Each constraint links
// the two vertices that lie opposite the shared edge in the
// adjacent triangles; rest distance comes from their initial
// positions, so a flat sheet starts at C = 0 and folding shrinks
// |p_k − p_l| below rest, raising C in magnitude until the
// constraint's compliance allows it (or projects exactly when α=0).
template<std::size_t N, Scalar T>
std::vector<XpbdBendingDistanceConstraint<N, T>>
build_bending_distance_constraints(
    const std::vector<XpbdParticle<N, T>>& parts,
    const std::vector<std::array<std::uint32_t, 3>>& faces,
    T compliance = T{0})
{
    using std::sqrt;
    std::vector<XpbdBendingDistanceConstraint<N, T>> out;

    // For each undirected edge (a,b), record the first opposite
    // vertex; the second face's opposite vertex closes the bending
    // pair (k, l) and the constraint between them is appended once.
    // Non-manifold edges (>2 incident faces) are not handled
    // specially: the third face's opposite vertex would silently
    // get matched with the first one, producing a duplicate
    // constraint. Cloth meshes are manifold by construction; if
    // a future caller feeds a soup, this helper should be revisited.
    std::unordered_map<detail::EdgeKey, std::uint32_t, detail::EdgeKeyHash> first_opp;
    first_opp.reserve(faces.size() * 3);

    auto add_edge = [&](std::uint32_t a, std::uint32_t b, std::uint32_t opp) {
        auto k = detail::edge_key(a, b);
        auto [it, inserted] = first_opp.try_emplace(k, opp);
        if (inserted) return;
        std::uint32_t kk = it->second, l = opp;
        Vec<T, N> diff = Vec<T, N>{parts[kk].x - parts[l].x};
        T rest = sqrt(diff.dot(diff));
        out.push_back({kk, l, rest, compliance, T{0}});
    };

    for (auto& f : faces) {
        add_edge(f[0], f[1], f[2]);
        add_edge(f[1], f[2], f[0]);
        add_edge(f[2], f[0], f[1]);
    }
    return out;
}

} // namespace spatium::physics::mechanics
