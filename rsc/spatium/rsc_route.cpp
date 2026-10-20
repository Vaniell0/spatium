// rsc_route <type> <operation> [--depth N] [--graph graph.json]
//
// The routes to an operation on a type of the library, read off the typed graph
// (rsc/spatium/graph.json): every distinct way, shortest first, each edge with the level the
// connectivity matrix verified it at. Exit 0 if there is a route, 1 if there is none, 2 on a bad
// call or an unreadable graph.

#include <rsc/core/graph.hpp>

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
    std::vector<std::string> positional;
    for (int i = 1; i < argc; ++i) {
        const std::string_view a = argv[i];
        if (a == "--depth" && i + 1 < argc) depth = std::strtoul(argv[++i], nullptr, 10);
        else if (a == "--graph" && i + 1 < argc) graph_path = argv[++i];
        else if (a.starts_with("--")) return usage();
        else positional.emplace_back(a);
    }
    if (positional.size() != 2) return usage();

    auto g = rsc::Graph::load(graph_path);
    if (!g) {
        std::fprintf(stderr, "rsc_route: %s\n", g.error().c_str());
        return 2;
    }
    const rsc::Node* type = g->find(positional[0]);
    if (!type || type->kind == rsc::NodeKind::Concept) {
        std::fprintf(stderr, "rsc_route: no type named '%s' in the graph\n", positional[0].c_str());
        return 2;
    }

    const auto routes = g->routes(type->id, positional[1], depth);
    std::printf("%s -> %s: %zu route%s (depth <= %zu)\n", type->id.c_str(), positional[1].c_str(), routes.size(),
                routes.size() == 1 ? "" : "s", depth);
    for (std::size_t r = 0; r < routes.size(); ++r) {
        std::printf("route %zu  level %d  %zu step%s\n", r + 1, g->route_level(routes[r]), routes[r].size(),
                    routes[r].size() == 1 ? "" : "s");
        for (const std::size_t i : routes[r]) {
            const rsc::Edge& e = g->edges()[i];
            std::printf("  %-12s [%s]  L%d", e.output.c_str(), e.variant.c_str(), e.level);
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
            std::printf("\n");
        }
    }
    return routes.empty() ? 1 : 0;
}
