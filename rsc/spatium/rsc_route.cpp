// rsc_route <type> <operation> [--depth N] [--graph graph.json]
// rsc_route --all --json [--depth N] [--graph graph.json]
//
// The routes to an operation on a type of the library, read off the typed graph
// (rsc/spatium/graph.json): every distinct way, shortest first, each edge with the level the
// connectivity matrix verified it at. Exit 0 if there is a route, 1 if there is none, 2 on a bad
// call or an unreadable graph.
//
// --json prints the routes as data (edge indices into the graph's edge list), and --all does it
// for every operation of every type: that is what rsc/spatium/rsc_find.py reads, so the question
// "by what routes" has one implementation, in rsc/core.

#include <rsc/core/graph.hpp>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <string_view>

#ifndef RSC_GRAPH_JSON
#  define RSC_GRAPH_JSON "rsc/spatium/graph.json"
#endif

namespace {

int usage() {
    std::fputs("usage: rsc_route <type> <operation> [--depth N] [--graph graph.json]\n"
               "  e.g. rsc_route 'Sphere<2>' distance\n", stderr);
    return 2;
}

std::string joined(const std::vector<std::string>& v, std::size_t from) {
    std::string s;
    for (std::size_t i = from; i < v.size(); ++i) s += (s.empty() ? "" : ", ") + v[i];
    return s;
}

}  // namespace

int main(int argc, char** argv) {
    std::string graph_path = RSC_GRAPH_JSON;
    std::size_t depth = 6;
    bool as_json = false, all = false;
    std::vector<std::string> positional;
    for (int i = 1; i < argc; ++i) {
        const std::string_view a = argv[i];
        if (a == "--depth" && i + 1 < argc) depth = std::strtoul(argv[++i], nullptr, 10);
        else if (a == "--graph" && i + 1 < argc) graph_path = argv[++i];
        else if (a == "--json") as_json = true;
        else if (a == "--all") all = true;
        else if (a.starts_with("--")) return usage();
        else positional.emplace_back(a);
    }
    if (positional.size() != (all ? 0u : 2u) || (all && !as_json)) return usage();

    auto g = rsc::Graph::load(graph_path);
    if (!g) {
        std::fprintf(stderr, "rsc_route: %s\n", g.error().c_str());
        return 2;
    }
    if (all) {
        // every (type, quantity) that has a route: {"depth": k, "routes": [{"type", "op", "routes": [[edge...]...]}...]}
        std::printf("{\"depth\": %zu, \"routes\": [", depth);
        bool first = true;
        for (const rsc::Node& n : g->nodes()) {
            if (n.kind != rsc::NodeKind::Type) continue;
            std::vector<std::string> ops;
            for (const rsc::Edge& e : g->edges())
                if (e.kind == rsc::EdgeKind::Operation && e.inputs[0] == n.id && std::ranges::find(ops, e.output) == ops.end())
                    ops.push_back(e.output);
            for (const std::string& op : ops) {
                const auto routes = g->routes(n.id, op, depth);
                if (routes.empty()) continue;
                std::printf("%s\n  {\"type\": \"%s\", \"op\": \"%s\", \"routes\": [", first ? "" : ",", n.id.c_str(), op.c_str());
                first = false;
                for (std::size_t r = 0; r < routes.size(); ++r) {
                    std::printf("%s[", r ? ", " : "");
                    for (std::size_t k = 0; k < routes[r].size(); ++k) std::printf("%s%zu", k ? ", " : "", routes[r][k]);
                    std::printf("]");
                }
                std::printf("]}");
            }
        }
        std::printf("\n]}\n");
        return 0;
    }
    const rsc::Node* type = g->find(positional[0]);
    if (!type || type->kind == rsc::NodeKind::Concept) {
        std::fprintf(stderr, "rsc_route: no type named '%s' in the graph\n", positional[0].c_str());
        return 2;
    }

    const auto routes = g->routes(type->id, positional[1], depth);
    if (as_json) {
        std::printf("{\"type\": \"%s\", \"op\": \"%s\", \"depth\": %zu, \"routes\": [", type->id.c_str(), positional[1].c_str(), depth);
        for (std::size_t r = 0; r < routes.size(); ++r) {
            std::printf("%s[", r ? ", " : "");
            for (std::size_t k = 0; k < routes[r].size(); ++k) std::printf("%s%zu", k ? ", " : "", routes[r][k]);
            std::printf("]");
        }
        std::printf("]}\n");
        return routes.empty() ? 1 : 0;
    }
    std::printf("%s -> %s: %zu route%s (depth <= %zu)\n", type->id.c_str(), positional[1].c_str(), routes.size(),
                routes.size() == 1 ? "" : "s", depth);
    for (std::size_t r = 0; r < routes.size(); ++r) {
        std::printf("route %zu  level %d  %zu step%s\n", r + 1, g->route_level(routes[r]), routes[r].size(),
                    routes[r].size() == 1 ? "" : "s");
        for (const std::size_t i : routes[r]) {
            const rsc::Edge& e = g->edges()[i];
            std::printf("  %-12s [%s]  L%d", e.output.c_str(), e.variant.c_str(), e.level);
            if (e.inputs[0] != type->id) std::printf("  on %s", e.inputs[0].c_str());
            if (!e.levels.empty()) {
                std::printf("  (");
                bool first = true;
                for (const auto& [scalar, lv] : e.levels) {
                    std::printf("%s%s L%d", first ? "" : ", ", scalar.c_str(), lv);
                    first = false;
                }
                std::printf(")");
            }
            if (e.inputs.size() > 1) std::printf("  needs %s", joined(e.inputs, 1).c_str());
            if (e.entry) std::printf("  (entry)");
            std::printf("\n");
        }
    }
    return routes.empty() ? 1 : 0;
}
