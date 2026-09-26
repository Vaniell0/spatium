#pragma once
// Continuous collision as a chain of primitives, so a search can put
// together the algorithm a class of queries wants instead of one general
// method paying for every case.
//
// The reason, measured in physics/mechanics/rigid_contact.hpp: advancement
// alone is cheap far from a surface and crawls at a graze; the search over
// parameters x time is tight at a graze and wasteful where advancement
// would have cleared the step in two moves; a closed form, where there is
// one, beats both. No single method is the cheap one everywhere, which is
// the whole problem a general CCD -- or IPC -- has with performance.
//
// A query is a ball moving by `disp` over the step past one obstacle. A
// chain runs its primitives in order over one State -- how far the step is
// proved clear, whether it is decided, what it has spent -- so a primitive
// starts where the one before it stopped, and the obstacle's cell tree is
// shared, so nothing is split twice:
//
//   CF       the closed form: a quadric's roots, a torus's quartic.
//            Decides the query; not available for a chart.
//   A(k)     k moves of advancement by the floor on distance.
//   S(b)     the search in parameters x time from where the step is proved
//            clear, with a budget of b new cells.
//   P(dt)    a speculative jump: [t, t + dt] is proved clear by one floor
//            at its middle, or nothing changes.
//
// Every primitive keeps the one-sided rule: none moves t past a moment it
// has not proved clear. A chain that ends undecided answers "contact at
// t" -- early, never late.

#include <spatium/geometry/ray_surface.hpp>
#include <spatium/physics/mechanics/rigid_contact.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace rsc::ccd {

using V3 = spatium::Vec<double, 3>;
namespace mech = spatium::physics::mechanics;

enum class Shape : std::uint8_t { Sphere, Torus, Chart };

struct Obstacle {
    Shape shape = Shape::Sphere;
    double radius = 1.0;                        // Sphere
    spatium::geometry::Torus<double> torus{};   // Torus
    // Every obstacle has a chart too, for the primitives that search one.
    std::shared_ptr<const spatium::ParametricSurface<double>> chart;
    mech::LipschitzChart<double> bound() const {
        if (shape == Shape::Sphere) {
            const double r = radius;
            return {*chart, r, [r](double u0, double u1, double v0, double v1) {
                        const double h = std::numbers::pi / 2;
                        const double s = (v0 <= h && h <= v1) ? 1.0 : std::max(std::sin(v0), std::sin(v1));
                        return r * s * (u1 - u0) / 2 + r * (v1 - v0) / 2;
                    }};
        }
        if (shape == Shape::Torus) {
            const double R = torus.major_radius, r = torus.minor_radius;
            return {*chart, R + r, [R, r](double u0, double u1, double v0, double v1) {
                        return (R + r) * (u1 - u0) / 2 + r * (v1 - v0) / 2;
                    }};
        }
        return {*chart, lipschitz, {}};
    }
    double lipschitz = 1.0;                     // Chart
};

struct Query {
    V3 p0, disp;
    double radius = 0;
    const Obstacle* obstacle = nullptr;
};

struct Answer {
    bool hit = false;
    double toi = 1.0;
    std::size_t cost = 0;   // evaluations: chart cells, distances, closed forms
};

struct State {
    double t = 0;           // [0, t) is proved clear
    bool done = false;
    Answer answer;
    std::unique_ptr<mech::ChartCellTree<double>> tree;
};

enum class Op : std::uint8_t { CF, A, S, P };
struct Step {
    Op op;
    double arg = 0;         // A: moves, S: budget, P: dt
};
using Chain = std::vector<Step>;

inline std::string name(const Step& s) {
    switch (s.op) {
        case Op::CF: return "CF";
        case Op::A: return "A" + std::to_string(static_cast<int>(s.arg));
        case Op::S: return "S" + std::to_string(static_cast<int>(std::log2(s.arg)));
        case Op::P: return "P" + std::to_string(static_cast<int>(std::lround(s.arg * 100)));
    }
    return "?";
}
inline std::string name(const Chain& c) {
    std::string s;
    for (const auto& st : c) s += (s.empty() ? "" : " ") + name(st);
    return s;
}

namespace detail {

inline double tolerance(const Query& q) {
    const double scale = std::max({q.p0.norm(), V3{q.p0 + q.disp}.norm(), q.radius});
    return std::sqrt(std::numeric_limits<double>::epsilon()) * scale;
}

// A floor on the distance from p to the obstacle's surface, and the best
// real distance, from the closed form where there is one.
inline std::pair<double, double> floor_at(const Query& q, State& s, const V3& p) {
    const auto& ob = *q.obstacle;
    if (ob.shape == Shape::Sphere) {
        ++s.answer.cost;
        const double d = std::abs(p.norm() - ob.radius);
        return {d, d};
    }
    if (ob.shape == Shape::Torus) {
        ++s.answer.cost;
        const double d = mech::point_to(p, ob.torus).distance;
        return {d, d};
    }
    if (!s.tree) s.tree = std::make_unique<mech::ChartCellTree<double>>(ob.bound());
    const auto before = s.tree->size();
    const auto r = mech::detail::tree_distance_bound(p, *s.tree, detail::tolerance(q) + q.radius, 1u << 20);
    s.answer.cost += s.tree->size() - before + 1;
    return r;
}

inline void decide(State& s, bool hit, double toi) {
    s.done = true;
    s.answer.hit = hit;
    s.answer.toi = toi;
}

}  // namespace detail

// One primitive, on the query's state.
inline void apply(const Query& q, State& s, const Step& st) {
    if (s.done) return;
    const double speed = q.disp.norm();
    const double tol = detail::tolerance(q);
    const auto& ob = *q.obstacle;
    switch (st.op) {
        case Op::CF: {
            if (ob.shape == Shape::Chart) return;   // no closed form: nothing to do
            ++s.answer.cost;
            const V3 a{q.p0 + q.disp * s.t};
            const V3 rest{q.disp * (1.0 - s.t)};
            const double len = rest.norm();
            if (len <= 0) { detail::decide(s, false, 1.0); return; }
            const spatium::geometry::Ray<3, double> ray{a, V3{rest * (1.0 / len)}};
            double first = std::numeric_limits<double>::infinity();
            if (ob.shape == Shape::Sphere) {
                const auto sphere = spatium::geometry::Quadric<double>::sphere(ob.radius + q.radius);
                for (const auto& h : spatium::geometry::ray_quadric(ray, sphere))
                    if (h.t >= 0) first = std::min(first, h.t);
            } else {
                auto grown = ob.torus;
                grown.minor_radius += q.radius;
                for (const auto& h : spatium::geometry::ray_torus(ray, grown))
                    if (h.t >= 0) first = std::min(first, h.t);
            }
            // The roots solve the surface exactly; to answer early, never
            // late, the contact is reported the tolerance before the root.
            if (first <= len) detail::decide(s, true, std::max(s.t, s.t + (first - tol) / len * (1.0 - s.t)));
            else detail::decide(s, false, 1.0);
            return;
        }
        case Op::A: {
            for (int k = 0; k < static_cast<int>(st.arg) && !s.done; ++k) {
                const V3 p{q.p0 + q.disp * s.t};
                const auto [lower, best] = detail::floor_at(q, s, p);
                if (best - q.radius <= tol) { detail::decide(s, true, s.t); return; }
                const double step = (lower - q.radius) / std::max(speed, 1e-300);
                if (step <= 0) return;
                if (s.t + step >= 1.0) { detail::decide(s, false, 1.0); return; }
                s.t += step;
            }
            return;
        }
        case Op::S: {
            if (!s.tree) s.tree = std::make_unique<mech::ChartCellTree<double>>(ob.bound());
            const V3 a{q.p0 + q.disp * s.t};
            const V3 rest{q.disp * (1.0 - s.t)};
            const auto r = mech::first_contact(a, q.radius, rest, *s.tree, static_cast<std::size_t>(st.arg), 0);
            s.answer.cost += r.evaluations;
            detail::decide(s, r.hit, s.t + r.toi * (1.0 - s.t));
            return;
        }
        case Op::P: {
            const double dt = std::min(st.arg, 1.0 - s.t);
            if (dt <= 0) return;
            const V3 mid{q.p0 + q.disp * (s.t + 0.5 * dt)};
            const auto [lower, best] = detail::floor_at(q, s, mid);
            (void)best;
            // Nothing within reach of the half-interval's travel: clear.
            if (lower - q.radius > speed * 0.5 * dt + tol) {
                s.t += dt;
                if (s.t >= 1.0) detail::decide(s, false, 1.0);
            }
            return;
        }
    }
}

// A chain run to its end: undecided means "contact at t", early and safe.
inline Answer run(const Query& q, const Chain& chain) {
    State s;
    for (const auto& st : chain) apply(q, s, st);
    if (!s.done) detail::decide(s, true, s.t);
    return s.answer;
}

}  // namespace rsc::ccd
