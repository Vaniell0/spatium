// Synthesis on the substrate: find a sequence of mesh operations that
// produces a mesh meeting a specification.
//
// The same search as `synthesize.cpp`, moved from numbers to geometry, and
// the move is justified by measurement rather than by preference. A chain
// over a scalar accumulator collapses -- 9 765 625 paths at depth 10 fold
// into 13 836 distinct states, an effective branching factor of 2.13 out
// of 5 -- so the space is exhausted for nothing and there is nothing for a
// search to be good at. Meshes collide too, since transforms compose, but
// about four times less: 4.66 out of 6, converged, measured in
// `mesh_collisions.cpp`. That is the difference between a space worth
// searching and one that is not.
//
// One thing genuinely changes with the move, and it is not a detail.
// **Over numbers the target is a value; over geometry it has to be a
// specification.** Nobody can hand a search the exact mesh they want --
// if they had it, they would not be searching. So the goal is stated as
// properties, and the oracle checks them by running the candidate. That is
// still free and still unfakeable: the properties are measured off the
// actual mesh, not asserted.
//
// It also makes refusals sharper. "No sequence of these operations, in
// this many steps, produces a mesh with those properties" is a statement
// about the operations, which is a more useful thing to learn than that a
// particular number was out of reach.
//
// Build: part of the rsc tools target.

#include "mesh_toy.hpp"

#include <format>
#include <limits>
#include <optional>
#include <print>
#include <set>
#include <string>
#include <vector>

namespace {

// What the caller wants, as properties rather than as a mesh.
struct Spec {
    std::string label;
    std::size_t min_faces = 0;
    double min_edge = 0.0;
    double max_edge = std::numeric_limits<double>::infinity();
};

double longest_edge(const meshtoy::MeshT& m) {
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

bool satisfies(const meshtoy::MeshT& m, const Spec& s) {
    if (m.face_count() < s.min_faces) return false;
    const double e = longest_edge(m);
    return e >= s.min_edge && e <= s.max_edge;
}

struct Node {
    meshtoy::MeshT mesh;
    std::vector<std::string> path;
};

struct Found {
    bool ok = false;
    std::vector<std::string> path;
    std::size_t expanded = 0;
    std::size_t faces = 0;
    double edge = 0.0;
};

Found search(const meshtoy::MeshT& start, const meshtoy::Space& space, const Spec& spec,
             const std::vector<std::string>& allowed, std::size_t max_steps) {
    Found res;
    std::vector<Node> frontier{Node{start, {}}};
    std::set<std::string> seen{meshtoy::fingerprint(start, 9)};

    if (satisfies(start, spec)) {
        res.ok = true;
        res.faces = start.face_count();
        res.edge = longest_edge(start);
        return res;
    }

    for (std::size_t depth = 1; depth <= max_steps; ++depth) {
        std::vector<Node> next;
        for (const auto& n : frontier) {
            for (const auto& op : meshtoy::ops()) {
                if (!allowed.empty() &&
                    std::find(allowed.begin(), allowed.end(), op.name) == allowed.end())
                    continue;
                ++res.expanded;

                meshtoy::MeshT m = op.apply(n.mesh, space);
                if (m.vertex_count() > meshtoy::kVertexCap) continue;

                if (satisfies(m, spec)) {
                    res.ok = true;
                    res.path = n.path;
                    res.path.push_back(op.name);
                    res.faces = m.face_count();
                    res.edge = longest_edge(m);
                    return res;
                }
                // Deduplicated at 1e-9 rather than bit-exactly: two meshes
                // that agree to that tolerance satisfy the same
                // specifications, so treating them as one is right here
                // even though the collision study reports both.
                if (!seen.insert(meshtoy::fingerprint(m, 9)).second) continue;

                Node child{std::move(m), n.path};
                child.path.push_back(op.name);
                next.push_back(std::move(child));
            }
        }
        frontier = std::move(next);
        if (frontier.empty()) break;
    }
    return res;
}

void run(const meshtoy::MeshT& start, const meshtoy::Space& space, const Spec& spec,
         const std::vector<std::string>& allowed, std::size_t max_steps) {
    std::string names;
    for (const auto& op : meshtoy::ops())
        if (allowed.empty() ||
            std::find(allowed.begin(), allowed.end(), op.name) != allowed.end())
            names += (names.empty() ? "" : ", ") + std::string(op.name);

    std::println("{}", spec.label);
    std::println("  allowed: {{{}}}   max steps: {}", names, max_steps);

    const Found f = search(start, space, spec, allowed, max_steps);
    if (f.ok) {
        std::string chain;
        for (const auto& s : f.path) chain += (chain.empty() ? "" : " -> ") + s;
        std::println("  found in {} steps, {} meshes built:\n    {}",
                     f.path.size(), f.expanded, chain.empty() ? "(already satisfied)" : chain);
        std::println("    result: {} faces, longest edge {:.4f}", f.faces, f.edge);
    } else {
        std::println("  no sequence of these operations reaches it within {} steps "
                     "({} meshes built)", max_steps, f.expanded);
    }
    std::println("");
}

}  // namespace

int main() {
    const meshtoy::Space space{};
    const meshtoy::MeshT start = spatium::mesh::icosahedron(space);

    std::println("Synthesis on the substrate: a specification, not a target mesh");
    std::println("start: icosahedron, {} faces, longest edge {:.4f}",
                 start.face_count(), longest_edge(start));
    std::println("");

    // The band is chosen from the measured subdivision ladder rather than
    // from an estimate, after an estimate got it wrong once: subdividing
    // an icosahedron on a sphere does not halve the longest edge, because
    // the new vertices are projected back onto the surface. Measured, the
    // ladder is 1.0515 -> 0.6180 -> 0.3249 -> 0.1646.
    //
    // So [0.40, 0.50] falls between two rungs and no amount of pure
    // subdivision lands in it, while 0.3249 * 1.5 = 0.4874 does.
    //
    // The restricted case below was written expecting a refusal, on the
    // argument that shifts and normal flips do not change edge lengths.
    // The search found `subdivide -> shift_x -> subdivide` instead, and
    // the argument was wrong: a shift does not change edge lengths by
    // itself, but `subdivide_once` projects its new vertices onto the
    // space, and a shifted mesh is no longer concentric with the sphere,
    // so the next subdivision projects differently. The effect lives in
    // the composition, not in either operation -- which is the kind of
    // thing a search finds and an argument does not. Left in place rather
    // than tuned away, because it is the more useful result.
    const Spec band{"1. at least 200 faces, longest edge in [0.40, 0.50]", 200, 0.40, 0.50};
    run(start, space, band, {}, 4);

    Spec band_restricted = band;
    band_restricted.label = "2. the same, with no scaling available";
    run(start, space, band_restricted, {"subdivide", "shift_x", "shift_y", "flip"}, 4);

    // The control. Without it a refusal proves nothing, since a harness
    // that never succeeds refuses everything.
    const Spec reachable{"3. the control: at least 200 faces, longest edge at most 0.60",
                         200, 0.0, 0.60};
    run(start, space, reachable, {"subdivide", "shift_x", "shift_y", "flip"}, 4);

    // A refusal whose closure can be proved rather than argued, after two
    // arguments in this file turned out wrong. `flip` reverses winding and
    // touches nothing else, so the face count is invariant under it: from
    // a 20-face icosahedron, "at least 200 faces" is unreachable by any
    // number of flips, and the search reports that rather than running out
    // of budget.
    const Spec impossible{"4. at least 200 faces, flips only -- unreachable by construction",
                          200, 0.0, std::numeric_limits<double>::infinity()};
    run(start, space, impossible, {"flip"}, 4);

    std::println("Case 2 was written expecting a refusal and is the most useful line here:");
    std::println("the search found a route the author had argued was impossible, because a");
    std::println("shift changes what the *next* subdivision projects. Case 4 is the refusal,");
    std::println("and its closure is provable rather than argued -- flips cannot change a");
    std::println("face count. A refusal only means anything next to a control that succeeds.");
    return 0;
}
