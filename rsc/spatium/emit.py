"""A route of the typed graph as C++: the code that runs it, written by a program.

A route is an ordered list of edges of rsc/spatium/graph.json, each after the edges that yield its
needs. Every edge becomes one function `e<k>` that does what the edge names, calling the functions
of the edges it needs; the last is the route. Nothing is written by hand per route.

  member      `s.exp_map(p, v, t)`                       the form the type states itself
  adl         `exp_map(s, p, v, t)`                       a free function found by argument-dependent lookup
  derived     `sqrt(metric(s, p, log_map(s, p, q), ...))`  access.hpp's third step, from the needs of the edge
  generic     geodesic = exp_p(t log_p q), midpoint = geodesic at 1/2
  components  a product from its factors' functions (factor spaces `s.space1`, `s.space2`)
  group       a Lie group manifold from its group (`s.group`)
  lift:<S>    the same quantity computed on the space over scalar S, the sample converted up and the
              answer back to double (tests/symmetry/rsc_pair_support.hpp: up_pt, up_tan, up_sc, down)

An edge counts one call (rsc_pair::calls()) when it runs; a lift counts none, it only converts. That
is the cost of a route: the operation calls one evaluation makes.

The independence rule (for the first version, as the plan says): two routes are independent when they
share no node. A node is an edge seen without its scalar: (type before '@', operation, variant). Scalar
lifts are conversions of arithmetic, which is trusted, so they are not nodes; what is computed on the
lifted space is, under the type's own name. Arithmetic, sqrt, exp, acos and cos are the trusted base
because every route ends in them. What the graph cannot see is what a node calls inside the library
(a product's member distance calls its factors' distance); the rule is structural, not semantic.
"""
import json
import pathlib
import re

# the arguments of each quantity: (name, kind) with kind pt (a point), tan (a tangent), sc (a scalar)
SIG = {
    "exp_map": [("p", "pt"), ("v", "tan"), ("t", "sc")],
    "log_map": [("p", "pt"), ("q", "pt")],
    "metric": [("p", "pt"), ("u", "tan"), ("v", "tan")],
    "distance": [("p", "pt"), ("q", "pt")],
    "norm_at": [("p", "pt"), ("v", "tan")],
    "project": [("p", "pt")],
    "normal": [("p", "pt")],
    "geodesic": [("p", "pt"), ("q", "pt"), ("t", "sc")],
    "midpoint": [("p", "pt"), ("q", "pt")],
}
# what a draw d (rsc_pair::Draw) is called with: points d.p, d.q, tangents d.u, d.v at p; the scalar of exp_map is 1, of geodesic 1/4
HARNESS_ARGS = {
    "exp_map": "d.p, d.u, 1.0", "log_map": "d.p, d.q", "metric": "d.p, d.u, d.v", "distance": "d.p, d.q",
    "norm_at": "d.p, d.v", "project": "d.p", "normal": "d.p", "geodesic": "d.p, d.q, 0.25", "midpoint": "d.p, d.q",
}
# what each path promises: see rsc_pair_support.hpp (allowance)
TIER = {"exp_map": "Smooth", "metric": "Smooth", "norm_at": "Smooth", "project": "Smooth", "normal": "Smooth",
        "log_map": "Inverse", "distance": "Inverse", "geodesic": "Inverse", "midpoint": "Inverse"}
MEMBER_NAME = {"exp_map": "exp_map", "log_map": "log_map", "metric": "metric_at", "distance": "distance",
               "project": "project", "normal": "normal", "compose": "compose", "inverse": "inverse",
               "lie_exp": "exp", "lie_log": "log"}
PILLS = "void exp_map() = delete; void log_map() = delete; void metric_at() = delete; void distance() = delete; void project() = delete; void normal() = delete;"

LIFT_CTYPE = {"Real50": "spatium::Real50", "Dual": "spatium::Dual<double>"}


def strip_scalar(type_id):
    return type_id.split("@")[0]


def file_stem(type_id):
    return re.sub(r"\W+", "_", type_id).strip("_")


class Graph:
    """graph.json as the emitter reads it."""

    def __init__(self, path):
        self.doc = json.loads(pathlib.Path(path).read_text())
        self.nodes = {n["id"]: n for n in self.doc["nodes"]}
        self.edges = self.doc["edges"]

    def edge(self, i):
        return self.edges[i]

    def subject(self, i):
        return self.edges[i]["inputs"][0]

    def is_lift(self, i):
        return self.edges[i]["variant"].startswith("lift")

    def want(self, i, k):
        """The need k of edge i as (type, quantity)."""
        w = self.edges[i]["needs"][k]
        if "#" in w:
            t, q = w.split("#", 1)
            return t, q
        return self.subject(i), w

    def nodes_of(self, route):
        """The nodes of a route (see the module note on independence): lifts are not nodes."""
        return frozenset((strip_scalar(self.subject(i)), self.edges[i]["op"], self.edges[i]["variant"])
                         for i in route if not self.is_lift(i))

    def describe(self, route, start):
        """The route in one line, dependencies first: `log_map[member] > metric[member] > distance[derived]`."""
        parts = []
        for i in route:
            e = self.edges[i]
            where = "" if self.subject(i) == start else f" on {self.subject(i)}"
            parts.append(f"{e['op']}[{e['variant']}]{where}")
        return " > ".join(parts)

    def short(self, route, start):
        """A name for a pair: the last edge's variant, and what it is built from when that is not obvious."""
        last = self.edges[route[-1]]
        if self.is_lift(route[-1]):
            inner = self.short(route[:-1], self.subject(route[-2]) if len(route) > 1 else start)
            return f"{last['variant']} of {inner}"
        if len(route) == 1:
            return last["variant"]
        subs = [self.edges[i]["variant"] + ("" if self.subject(i) == self.subject(route[-1]) else f" on {strip_scalar(self.subject(i))}")
                for i in route[:-1]]
        return last["variant"] + " (" + ", ".join(f"{self.edges[i]['op']} {s}" for i, s in zip(route[:-1], subs)) + ")"

    def cell_of(self, type_id):
        return self.nodes[type_id].get("cell")


# ── One edge as a function ─────────────────────────────────────

def _args(op):
    return [a for a, _ in SIG[op]]


def emit_edge(g, route, k, fn):
    """The C++ function e<k> for the k-th edge of the route. `fn` maps (type, quantity) -> index in the route."""
    i = route[k]
    e = g.edge(i)
    subj, op, variant = g.subject(i), e["op"], e["variant"]
    name = f"e{k}"
    needs = [f"e{fn[g.want(i, n)]}" for n in range(len(e["needs"]))]
    count = "++rsc_pair::calls();"

    if variant == "member" or variant == "adl":
        call = f"s.{MEMBER_NAME[op]}(a...)" if variant == "member" else f"{MEMBER_NAME[op]}(s, a...)"
        return f"inline auto {name}(const auto& s, const auto&... a) {{ {count} return rsc_pair::val({call}); }}"

    sig = SIG.get(op)
    params = ", ".join(f"const auto& {a}" for a in _args(op)) if sig else ""
    head = f"inline auto {name}(const auto& s{', ' + params if params else ''})"

    if variant == "derived" and op == "distance":
        log, metric = needs
        return (f"{head} {{ {count} using std::sqrt; const auto v = {log}(s, p, q); return rsc_pair::val(sqrt({metric}(s, p, v, v))); }}")
    if variant == "derived" and op == "norm_at":
        (metric,) = needs
        return f"{head} {{ {count} using std::sqrt; return rsc_pair::val(sqrt({metric}(s, p, v, v))); }}"
    if variant == "generic" and op == "geodesic":
        exp, log = needs
        return (f"{head} {{ {count} using Tv = spatium::spaces::tangent_t<std::remove_cvref_t<decltype(s)>>; "
                f"return rsc_pair::val({exp}(s, p, Tv{{{log}(s, p, q)}}, t)); }}")
    if variant == "generic" and op == "midpoint":
        (geo,) = needs
        return (f"{head} {{ {count} using T = spatium::spaces::scalar_t<std::remove_cvref_t<decltype(s)>>; "
                f"return rsc_pair::val({geo}(s, p, q, T{{1}} / T{{2}})); }}")

    if variant == "components":
        a, b = needs[0], needs[-1]   # one need when the two factors are the same type
        pre = "using S = std::remove_cvref_t<decltype(s)>; "
        f1 = lambda part, x: f"S::template first_of<typename S::{part}>({x})"
        f2 = lambda part, x: f"S::template second_of<typename S::{part}>({x})"
        if op == "distance":
            body = (f"using std::sqrt; const auto d1 = {a}(s.space1, {f1('P1', 'p')}, {f1('P1', 'q')}); "
                    f"const auto d2 = {b}(s.space2, {f2('P2', 'p')}, {f2('P2', 'q')}); return rsc_pair::val(sqrt(d1 * d1 + d2 * d2));")
        elif op == "exp_map":
            body = (f"return S::template join_as<typename S::PointType>("
                    f"typename S::P1({a}(s.space1, {f1('P1', 'p')}, {f1('V1', 'v')}, t)), "
                    f"typename S::P2({b}(s.space2, {f2('P2', 'p')}, {f2('V2', 'v')}, t)));")
        elif op == "log_map":
            body = (f"return S::template join_as<typename S::TangentVector>("
                    f"typename S::V1({a}(s.space1, {f1('P1', 'p')}, {f1('P1', 'q')})), "
                    f"typename S::V2({b}(s.space2, {f2('P2', 'p')}, {f2('P2', 'q')})));")
        elif op == "metric":
            body = (f"return rsc_pair::val({a}(s.space1, {f1('P1', 'p')}, {f1('V1', 'u')}, {f1('V1', 'v')}) + "
                    f"{b}(s.space2, {f2('P2', 'p')}, {f2('V2', 'u')}, {f2('V2', 'v')}));")
        else:
            raise ValueError(f"no components emission for {op}")
        return f"{head} {{ {count} {pre}{body} }}"

    if variant == "group":
        pre = "using S = std::remove_cvref_t<decltype(s)>; "
        if op == "exp_map":
            compose, lie_exp = needs
            body = f"return {compose}(s.group, p, {lie_exp}(s.group, typename S::TangentVector{{v * t}}));"
        elif op == "log_map":
            compose, inverse, lie_log = needs
            body = f"return typename S::TangentVector{{{lie_log}(s.group, {compose}(s.group, {inverse}(s.group, p), q))}};"
        else:
            raise ValueError(f"no group emission for {op}")
        return f"{head} {{ {count} {pre}{body} }}"

    if variant.startswith("lift:"):
        suffix = variant.split(":", 1)[1]
        (inner,) = needs
        lifted_type = g.want(i, 0)[0]
        cell = g.nodes[lifted_type]["cell"]
        conv = []
        for a, kind in SIG[op]:
            conv.append({"pt": f"rsc_pair::up_pt<R>(sl, {a})", "tan": f"rsc_pair::up_tan<R>({a})", "sc": f"rsc_pair::up_sc<R>({a})"}[kind])
        return (f"{head} {{ using R = {LIFT_CTYPE[suffix]}; static const auto sl = {cell}::make<R>(); "
                f"return rsc_pair::down({inner}(sl, {', '.join(conv)})); }}")

    raise ValueError(f"no emission for {op} [{variant}]")


def emit_route(g, route):
    """The functions of a route, to be placed in a namespace of their own. The route is e<len-1>."""
    fn = {}
    for k, i in enumerate(route):
        fn[(g.subject(i), g.edge(i)["output"])] = k
    lines = [PILLS]
    lines += [emit_edge(g, route, k, fn) for k in range(len(route))]
    return "\n".join(lines)


def entry_name(route):
    return f"e{len(route) - 1}"


def call(ns, route, op):
    """The call of the route on a draw `d` in a context with the space `space`."""
    return f"{ns}::{entry_name(route)}(space, {HARNESS_ARGS[op]})"


def bound_expr(op, cell):
    return f"rsc_pair::bound<{cell}>(rsc_pair::Tier::{TIER[op]}, d)"
