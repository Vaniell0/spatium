// Does a chain found on one space still work on another?
//
// This is the claim the whole RSC-as-search direction rests on, and the
// one that is not a compute question: a chain of concept-constrained
// operations is a program over the *concept*, so running it on a sphere
// and on a torus should be two instantiations of one algorithm rather
// than two algorithms. AlphaTensor's discoveries are bound to the ring
// they were found in and do not generalise themselves; this is the axis
// nobody else has, because nobody else has one interface through which
// the generalisation could even be checked.
//
// Checked here by compiler and by measurement rather than by argument.
// The operations below are templates over `Surface S`; the search runs on
// `Sphere<2>`; the chain it returns is then replayed, unchanged, on
// `Euclidean<3>` and on a `ParametricSurface` torus. Those three behave
// genuinely differently under the same operation -- subdivision projects
// new vertices onto the sphere, does nothing of the kind in Euclidean
// space, and projects onto a torus in the third -- so a chain that
// survives all three survived something.
//
// The prediction, from the design rather than from hope: **structure
// transfers, constants do not.** A specification in absolute units is
// about the space it was written for, and should fail elsewhere. A
// relative one should hold. That is the same dispatch/calibration split
// the RSC design already draws -- dispatch is learned and space-agnostic,
// calibration is classical and re-fitted per deployment -- arriving here
// on its own, which is the sort of agreement worth taking seriously.
//
// Own generic op table rather than mesh_toy.hpp's: that one is deliberately
// concrete, because the collision study needs a single fixed space to
// count states in. Here being generic is the entire point.
//
// Build: part of the rsc tools target.

#include <spatium/mesh/mesh.hpp>
#include <spatium/mesh/operations.hpp>
#include <spatium/mesh/primitives.hpp>
#include <spatium/mesh/subdivision.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <spatium/spaces/parametric.hpp>
#include <spatium/spaces/sphere.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <print>
#include <string>
#include <vector>

namespace {

// The action set, over any Surface. Nothing here names a space.
template<spatium::Surface S>
struct Ops {
    using M = spatium::mesh::Mesh<S>;
    using Fn = M (*)(const M&, const S&);

    static M subdivide(const M& m, const S& s) { return spatium::mesh::subdivide_once(m, s); }
    static M scale_up(const M& m, const S&) {
        return spatium::mesh::transform(m, [](const auto& p) { return decltype(p){p * 1.5}; });
    }
    static M scale_down(const M& m, const S&) {
        return spatium::mesh::transform(m, [](const auto& p) {
            return decltype(p){p * (1.0 / 1.5)};
        });
    }
    static M shift_x(const M& m, const S&) {
        return spatium::mesh::transform(m, [](const auto& p) {
            auto q = p; q[0] += 0.25; return q;
        });
    }
    static M flip(const M& m, const S&) { return spatium::mesh::flip_normals(m); }

    static const std::vector<std::pair<const char*, Fn>>& all() {
        static const std::vector<std::pair<const char*, Fn>> v{
            {"subdivide", &subdivide}, {"scale_up", &scale_up},
            {"scale_down", &scale_down}, {"shift_x", &shift_x}, {"flip", &flip}};
        return v;
    }
};

template<spatium::Surface S>
double longest_edge(const spatium::mesh::Mesh<S>& m) {
    double worst = 0.0;
    for (const auto& f : m.faces)
        for (int k = 0; k < 3; ++k) {
            const auto& a = m.vertices[f[static_cast<std::size_t>(k)]];
            const auto& b = m.vertices[f[static_cast<std::size_t>((k + 1) % 3)]];
            double d2 = 0.0;
            for (int c = 0; c < 3; ++c) {
                const double t = a[static_cast<std::size_t>(c)] - b[static_cast<std::size_t>(c)];
                d2 += t * t;
            }
            worst = std::max(worst, std::sqrt(d2));
        }
    return worst;
}

// Absolute: the numbers mean what they say. Relative: they are multiples
// of whatever the starting mesh happens to be, so the same words describe
// a different requirement on every space.
struct Spec {
    bool relative;
    double face_factor;   // relative: multiple of start faces; absolute: a count
    double edge_lo, edge_hi;
};

template<spatium::Surface S>
bool satisfies(const spatium::mesh::Mesh<S>& m, const spatium::mesh::Mesh<S>& start,
               const Spec& spec) {
    const double e = longest_edge<S>(m);
    if (spec.relative) {
        const double e0 = longest_edge<S>(start);
        return static_cast<double>(m.face_count()) >=
                   spec.face_factor * static_cast<double>(start.face_count()) &&
               e >= spec.edge_lo * e0 && e <= spec.edge_hi * e0;
    }
    return static_cast<double>(m.face_count()) >= spec.face_factor &&
           e >= spec.edge_lo && e <= spec.edge_hi;
}

template<spatium::Surface S>
std::vector<std::string> search(const spatium::mesh::Mesh<S>& start, const S& space,
                                const Spec& spec, std::size_t max_steps) {
    struct Node { spatium::mesh::Mesh<S> mesh; std::vector<std::string> path; };
    std::vector<Node> frontier{Node{start, {}}};
    for (std::size_t d = 1; d <= max_steps; ++d) {
        std::vector<Node> next;
        for (const auto& n : frontier)
            for (const auto& [name, fn] : Ops<S>::all()) {
                auto m = fn(n.mesh, space);
                if (m.vertex_count() > 6000) continue;
                auto path = n.path;
                path.push_back(name);
                if (satisfies<S>(m, start, spec)) return path;
                next.push_back(Node{std::move(m), std::move(path)});
            }
        frontier = std::move(next);
        if (frontier.empty()) break;
    }
    return {};
}

// Replay, not re-search. The chain is fixed; only the space changes.
template<spatium::Surface S>
void replay(const char* space_name, const spatium::mesh::Mesh<S>& start, const S& space,
            const std::vector<std::string>& chain, const Spec& spec) {
    auto m = start;
    for (const auto& step : chain)
        for (const auto& [name, fn] : Ops<S>::all())
            if (step == name) { m = fn(m, space); break; }

    const bool ok = satisfies<S>(m, start, spec);
    std::println("    {:<22} start {:4} faces / edge {:7.4f}  ->  {:5} faces / edge {:7.4f}   {}",
                 space_name, start.face_count(), longest_edge<S>(start),
                 m.face_count(), longest_edge<S>(m), ok ? "holds" : "FAILS");
}

spatium::ParametricSurface<double> torus() {
    using T = double;
    return spatium::ParametricSurface<T>(
        [](T u, T v) {
            const T R = 1.0, r = 0.35;
            return spatium::Vec<T, 3>{(R + r * std::cos(v)) * std::cos(u),
                                      (R + r * std::cos(v)) * std::sin(u), r * std::sin(v)};
        },
        {0.0, 2.0 * std::numbers::pi, 0.0, 2.0 * std::numbers::pi}, true, true);
}

}  // namespace

int main() {
    using Sph = spatium::Sphere<2, double>;
    using Euc = spatium::Euclidean<3, double>;
    using Par = spatium::ParametricSurface<double>;

    const Sph sphere{};
    const Euc euclid{};
    const Par tor = torus();

    const auto sphere_start = spatium::mesh::icosahedron(sphere);
    const auto euclid_start = spatium::mesh::box_mesh<double>();
    const auto torus_start = spatium::mesh::parametric_mesh(tor, 12, 8);

    std::println("Transfer: search on one space, replay the chain on two others");
    std::println("");

    // Relative specification -- the same sentence on every space.
    const Spec rel{true, 10.0, 0.0, 0.45};
    std::println("A. relative spec: at least 10x the starting faces, longest edge at most");
    std::println("   0.45x the starting longest edge");
    const auto chain_rel = search<Sph>(sphere_start, sphere, rel, 4);
    if (chain_rel.empty()) {
        std::println("    no chain found on the sphere; nothing to transfer");
    } else {
        std::string s;
        for (const auto& c : chain_rel) s += (s.empty() ? "" : " -> ") + c;
        std::println("    found on Sphere<2>: {}", s);
        replay<Sph>("Sphere<2>", sphere_start, sphere, chain_rel, rel);
        replay<Euc>("Euclidean<3>", euclid_start, euclid, chain_rel, rel);
        replay<Par>("ParametricSurface", torus_start, tor, chain_rel, rel);
    }
    std::println("");

    // Absolute specification -- the numbers describe the sphere.
    const Spec abs{false, 200.0, 0.40, 0.50};
    std::println("B. absolute spec: at least 200 faces, longest edge in [0.40, 0.50]");
    const auto chain_abs = search<Sph>(sphere_start, sphere, abs, 4);
    if (chain_abs.empty()) {
        std::println("    no chain found on the sphere; nothing to transfer");
    } else {
        std::string s;
        for (const auto& c : chain_abs) s += (s.empty() ? "" : " -> ") + c;
        std::println("    found on Sphere<2>: {}", s);
        replay<Sph>("Sphere<2>", sphere_start, sphere, chain_abs, abs);
        replay<Euc>("Euclidean<3>", euclid_start, euclid, chain_abs, abs);
        replay<Par>("ParametricSurface", torus_start, tor, chain_abs, abs);
    }

    std::println("");
    std::println("The three spaces do different things under the same operation --");
    std::println("subdivision projects onto the sphere, projects onto nothing in");
    std::println("Euclidean space, and projects onto the torus -- so a chain holding on");
    std::println("all three held against something rather than against a relabelling.");
    return 0;
}
