// The typed graph of the library (rsc/spatium/graph.json, made by rsc/spatium/gen_graph.py) and the
// first thing that reads it: a query for the routes to an operation.
//
// The graph is generated and committed, so it can go stale; the pair that catches that is declared
// against observed -- what the JSON says a type satisfies, against the compiler's own answer for the
// same type and concept, asked here, in this translation unit, not in the generator's.

#include <catch2/catch_test_macros.hpp>

#include <rsc/core/graph.hpp>

#include <spatium/algebra/dual.hpp>
#include <spatium/algebra/groups/so3.hpp>
#include <spatium/core/access.hpp>
#include <spatium/core/concepts.hpp>
#include <spatium/discrete/finite_set.hpp>
#include <spatium/geometry/geometry.hpp>
#include <spatium/spaces/euclidean.hpp>
#include <spatium/spaces/hyperbolic.hpp>
#include <spatium/spaces/lie_group.hpp>
#include <spatium/spaces/metric_chart.hpp>
#include <spatium/spaces/minkowski.hpp>
#include <spatium/spaces/product.hpp>
#include <spatium/spaces/sphere.hpp>
#if SPATIUM_HAS_BOOST_MULTIPRECISION
#  include <spatium/core/precision.hpp>
#endif

#include <algorithm>
#include <set>
#include <string>
#include <vector>

#ifndef RSC_GRAPH_JSON
#  error "RSC_GRAPH_JSON must name rsc/spatium/graph.json (tests/CMakeLists.txt sets it)"
#endif

namespace {

using namespace spatium;

const rsc::Graph& graph() {
    static const rsc::Graph g = [] {
        auto r = rsc::Graph::load(RSC_GRAPH_JSON);
        INFO((r.has_value() ? std::string() : r.error()));
        REQUIRE(r.has_value());
        return std::move(*r);
    }();
    return g;
}

struct Observation {
    std::string concept_id;
    bool observed;
};

// One concept the test asks the compiler about, under the id the graph gives it.
#define OBS(ID, CONCEPT, T) Observation{ID, CONCEPT<T>}

template<class T>
std::vector<Observation> observe() {
    return {
        OBS("Scalar", spatium::Scalar, T),
        OBS("Set", spatium::Set, T),
        OBS("TopologicalSpace", spatium::TopologicalSpace, T),
        OBS("MetricSpace", spatium::MetricSpace, T),
        OBS("VectorSpace", spatium::VectorSpace, T),
        OBS("NormedSpace", spatium::NormedSpace, T),
        OBS("InnerProductSpace", spatium::InnerProductSpace, T),
        OBS("HilbertSpace", spatium::HilbertSpace, T),
        OBS("EuclideanSpace", spatium::EuclideanSpace, T),
        OBS("Manifold", spatium::Manifold, T),
        OBS("RiemannianManifold", spatium::RiemannianManifold, T),
        OBS("LorentzianManifold", spatium::LorentzianManifold, T),
        OBS("Surface", spatium::Surface, T),
        OBS("Shape", spatium::geometry::Shape, T),
        OBS("ClosedShape", spatium::geometry::ClosedShape, T),
        OBS("Measurable", spatium::geometry::Measurable, T),
        OBS("Bounded", spatium::geometry::Bounded, T),
        OBS("Group", spatium::algebra::Group, T),
        OBS("Exponential", spatium::spaces::Exponential, T),
        OBS("Logarithmic", spatium::spaces::Logarithmic, T),
        OBS("Metrized", spatium::spaces::Metrized, T),
        OBS("Riemannian", spatium::spaces::Riemannian, T),
        OBS("Distanced", spatium::spaces::Distanced, T),
    };
}

template<class T>
void agree(const char* type_id) {
    const rsc::Graph& g = graph();
    INFO("type " << type_id);
    REQUIRE(g.find(type_id) != nullptr);
    std::size_t asked = 0;
    for (const Observation& o : observe<T>()) {
        INFO("concept " << o.concept_id);
        REQUIRE(g.find(o.concept_id) != nullptr);
        CHECK(g.satisfies(type_id, o.concept_id) == o.observed);
        ++asked;
    }
    CHECK(asked > 20);
}

const rsc::Edge* edge_with(const rsc::Graph& g, const rsc::Route& r, std::string_view output,
                           std::string_view variant) {
    for (const std::size_t i : r) {
        const rsc::Edge& e = g.edges()[i];
        if (e.output == output && e.variant == variant) return &e;
    }
    return nullptr;
}

}  // namespace

TEST_CASE("rsc graph: the library's concepts, types and edges are in it", "[rsc]") {
    const rsc::Graph& g = graph();
    CHECK(g.count(rsc::NodeKind::Concept) >= 70);
    CHECK(g.count(rsc::NodeKind::Type) >= 40);
    CHECK(g.count(rsc::NodeKind::Scalar) == 4);
    CHECK(g.count(rsc::EdgeKind::Refines) > 0);
    CHECK(g.count(rsc::EdgeKind::Satisfies) > 0);
    CHECK(g.count(rsc::EdgeKind::Operation) > 0);
    CHECK(g.find("no such node") == nullptr);
    CHECK(g.find("Sphere<2>") == g.find("S2"));   // an alias is the same node
}

TEST_CASE("rsc graph: what the JSON says a type satisfies is what the compiler says", "[rsc]") {
    agree<Euclidean<3>>("Euclidean<3>");
    agree<Sphere<2>>("Sphere<2>");
    agree<Hyperbolic<2>>("Hyperbolic<2>");
    agree<Minkowski<4>>("Minkowski<4>");
    agree<ProductSpace<Sphere<2>, Euclidean<1>>>("S2xE1");
    agree<LieGroupManifold<algebra::SO3<>>>("LieGroupManifold<SO3>");
    agree<geometry::Triangle<3>>("Triangle<3>");
    agree<geometry::Box<3>>("Box<3>");
    agree<geometry::Line<3>>("Line<3>");
    agree<discrete::FiniteSet<int>>("FiniteSet<int>");
    agree<double>("double");
    agree<float>("float");
    agree<Dual<double>>("Dual<double>");
#if SPATIUM_HAS_BOOST_MULTIPRECISION
    agree<Real50>("Real50");
#endif
}

TEST_CASE("rsc graph: a type satisfies everything its concepts refine", "[rsc]") {
    // declared against declared: the refinements are parsed from the headers, the satisfactions
    // are compiled, and a type that satisfies X must satisfy every concept X is built from.
    const rsc::Graph& g = graph();
    std::size_t checked = 0;
    for (const rsc::Edge& s : g.edges()) {
        if (s.kind != rsc::EdgeKind::Satisfies) continue;
        for (const rsc::Edge& r : g.edges()) {
            if (r.kind != rsc::EdgeKind::Refines || r.from != s.to) continue;
            INFO(s.from << " satisfies " << s.to << " which refines " << r.to);
            CHECK(g.satisfies(s.from, r.to));
            ++checked;
        }
    }
    CHECK(checked > 100);
}

TEST_CASE("rsc graph: distance on Sphere<2> has a closed form and a derived route, both verified", "[rsc]") {
    const rsc::Graph& g = graph();
    const std::vector<rsc::Route> routes = g.routes("Sphere<2>", "distance", 6);
    REQUIRE(routes.size() >= 2);

    const rsc::Route* member = nullptr;
    const rsc::Route* derived = nullptr;
    for (const rsc::Route& r : routes) {
        if (edge_with(g, r, "distance", "member")) member = &r;
        if (edge_with(g, r, "distance", "derived")) derived = &r;
    }
    REQUIRE(member != nullptr);
    REQUIRE(derived != nullptr);
    CHECK(member->size() == 1);                       // the member alone
    CHECK(derived->size() == 3);                      // log_map, metric, then |log|_g
    CHECK(g.route_level(*member) == 3);               // L3 in the connectivity matrix
    CHECK(g.route_level(*derived) == 3);

    // the routes are distinct, and each one ends in the goal
    std::set<rsc::Route> distinct(routes.begin(), routes.end());
    CHECK(distinct.size() == routes.size());
    for (const rsc::Route& r : routes) CHECK(g.edges()[r.back()].output == "distance");
}

namespace {
// The metric of the sphere in (theta, phi): a space given by its metric alone, which has no
// member distance (exp, log and distance are all derived by metric_chart).
struct SphereMetric {
    template<class S> Matrix<S, 2, 2> operator()(const Vec<S, 2>& x) const {
        using std::sin;
        Matrix<S, 2, 2> m{};
        m(0, 0) = S(1);
        m(1, 1) = sin(x[0]) * sin(x[0]);
        return m;
    }
};
template<class S>
concept HasMemberDistance = requires(const S& s, const spaces::point_t<S>& p) { s.distance(p, p); };
}  // namespace

TEST_CASE("rsc graph: a space with no member distance gets only the derived route", "[rsc]") {
    // observed: the type really has no distance of its own
    const auto chart = spaces::metric_chart<double, 2>(SphereMetric{});
    using Chart = std::remove_const_t<decltype(chart)>;
    STATIC_REQUIRE(!HasMemberDistance<Chart>);
    STATIC_REQUIRE(spaces::Distanced<Chart>);          // yet the customization point answers

    // declared: the graph offers the derived route and nothing else
    const rsc::Graph& g = graph();
    const std::vector<rsc::Route> routes = g.routes("S2Derived", "distance", 6);
    REQUIRE(routes.size() == 1);
    CHECK(edge_with(g, routes[0], "distance", "derived") != nullptr);
    CHECK(edge_with(g, routes[0], "distance", "member") == nullptr);
    CHECK(edge_with(g, routes[0], "distance", "adl") == nullptr);
}

TEST_CASE("rsc graph: routes from a node with nothing to offer are empty", "[rsc]") {
    const rsc::Graph& g = graph();
    CHECK(g.routes("double", "distance", 6).empty());       // a scalar has no space operations
    CHECK(g.routes("Triangle<3>", "distance", 6).empty());
    CHECK(g.routes("no such node", "distance", 6).empty());
    CHECK(g.routes("Sphere<2>", "no such operation", 6).empty());
    CHECK(g.routes("Sphere<2>", "distance", 0).empty());     // depth 0 reaches nothing
}

TEST_CASE("rsc graph: a graph that contradicts itself or is not JSON is refused", "[rsc]") {
    // every check of the loader gets a case it must fail
    const auto refused = [](std::string_view text, std::string_view why) {
        const auto g = rsc::Graph::parse(text);
        INFO("expected a refusal: " << why);
        REQUIRE(!g.has_value());
        return g.error();
    };
    const std::string head = R"({"version": 1, "nodes": [
        {"id": "T", "kind": "type"}, {"id": "C", "kind": "concept"}, {"id": "D", "kind": "concept"}],)";

    CHECK(rsc::Graph::parse(head + R"( "edges": []})").has_value());
    refused("{", "truncated");
    refused("[1, 2]", "not an object");
    refused(R"({"version": 2, "nodes": [], "edges": []})", "unknown version");
    refused(head + R"( "edges": [{"kind": "satisfies", "from": "T", "to": "Nowhere"}]})", "edge to no node");
    refused(head + R"( "edges": [{"kind": "satisfies", "from": "C", "to": "D"}]})", "a concept satisfies");
    refused(head + R"( "edges": [{"kind": "refines", "from": "T", "to": "C"}]})", "a type refines");
    refused(R"({"version": 1, "nodes": [{"id": "A", "kind": "type"}, {"id": "A", "kind": "type"}], "edges": []})",
            "one name twice");
    refused(head + R"( "edges": [{"kind": "operation", "op": "f", "variant": "v", "output": "x",
                     "inputs": ["T", "C"], "needs": [], "level": 1}]})", "an operation needs a concept the type does not satisfy");
    refused(head + R"( "edges": [{"kind": "operation", "op": "f", "variant": "v", "output": "x",
                     "inputs": ["T"], "needs": ["x"], "level": 1}]})", "an operation needs its own output");
    refused(head + R"( "edges": [{"kind": "operation", "op": "f", "variant": "v", "output": "x",
                     "inputs": ["T"], "needs": [], "level": 4}]})", "a level out of range");
}

TEST_CASE("rsc graph: routes close over their needs and are found shortest first", "[rsc]") {
    // x is yielded two ways, y needs x, z needs x and y: a diamond in which x has two routes
    const auto g = rsc::Graph::parse(R"({"version": 1,
        "nodes": [{"id": "T", "kind": "type"}],
        "edges": [
          {"kind": "operation", "op": "x", "variant": "a", "output": "x", "inputs": ["T"], "needs": [], "level": 3},
          {"kind": "operation", "op": "x", "variant": "b", "output": "x", "inputs": ["T"], "needs": [], "level": 1},
          {"kind": "operation", "op": "y", "variant": "a", "output": "y", "inputs": ["T"], "needs": ["x"], "level": 2},
          {"kind": "operation", "op": "z", "variant": "a", "output": "z", "inputs": ["T"], "needs": ["x", "y"], "level": 3}]})");
    REQUIRE(g.has_value());
    const auto routes = g->routes("T", "z", 5);
    REQUIRE(routes.size() == 2);                         // x by a, or x by b; y and z once each
    for (const rsc::Route& r : routes) {
        REQUIRE(r.size() == 3);
        CHECK(g->edges()[r[0]].output == "x");            // dependencies first
        CHECK(g->edges()[r[1]].output == "y");
        CHECK(g->edges()[r[2]].output == "z");
    }
    CHECK(g->route_level(routes[0]) + g->route_level(routes[1]) == 2 + 1);   // min of 3,2,3 and of 1,2,3
    CHECK(g->routes("T", "z", 2).empty());               // three edges do not fit in a depth of two
    CHECK(g->routes("T", "x", 1).size() == 2);
}
