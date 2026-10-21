#!/usr/bin/env python3
"""The library as a typed graph: rsc/spatium/graph.json.

Nodes
  concept   every `concept` declared under include/spatium/ (parsed from the headers)
  type      the library's spaces (every space of the connectivity matrix, Minkowski), the
            geometry primitives, SO3/SE3 and LieGroupManifold over them, FiniteSet
  scalar    double, float, Real50, Dual<double>

Edges
  refines    concept -> concept   parsed: `A<S>` in the body of B means B is built from A
  satisfies  type -> concept      COMPILED, not parsed: one translation unit per type prints
                                  `Concept<T>` for every concept that takes one type; nothing is
                                  compiled to be expected to fail, a false is a printed 0
  operation  one way to get a quantity out of a type, by what the type provides
               exp_map log_map metric project normal distance   member, ADL (a free function
                                  found by argument-dependent lookup, access.hpp's second step)
               distance, norm_at   derived: |log_p q|_g, sqrt(g(v, v))   (access.hpp's third step)
               geodesic midpoint frechet_mean verify_exp_log   the generic algorithms, which ask
                                  for the capabilities of core/access.hpp
             The variants are separate edges; each carries the level the connectivity matrix
             (tests/connectivity/expected.hpp) reached for it on the type: 0 declared and not
             exercised by any cell, 1 runs and is finite, 2 its axiom or symmetry holds, 3 its
             numbers agree across scalars. An edge is as verified as the least verified probe
             that exercises it.

Edges between types (what lets a route cross them; each is verified by compilation, in the translation
unit of the type it is on, or by the finder that compiles and runs it):
  scalar lift   `T` -> `T@Real50`, `T@Dual`: the same space over another scalar is a node of its own,
                with its own satisfies and operation edges. An `entry` edge `lift:<scalar>` on T needs the
                quantity on the lifted node ("compute on Real50 and convert"); it may close a route, never
                be a step inside one. Real50 is left out for the spaces whose Real50 cells cost minutes
                (gen_connectivity.HEAVY_SPACES), the erased surfaces (SLOW_REAL50) and the typed surfaces the
                matrix holds against long double.
  components    a ProductSpace from its factors: `components` edges for exp_map, log_map, metric and
                distance need the factor's quantity on each factor node. The factor table is checked
                against the type by a static_assert.
  group         a LieGroupManifold from its group: `group` edges for exp_map and log_map need the group's
                compose / inverse / lie_exp / lie_log, which are operations of the group node (member, detected).
  The matrix has no cell for these edges, so they sit at level 0 until a finder run witnesses them.

Run inside `nix develop` from the repository root:
    python3 rsc/spatium/gen_graph.py [-j 6] [--check] [--stats]
--check regenerates into a temporary directory and fails if rsc/spatium/graph.json would change.
"""
import argparse
import concurrent.futures as cf
import json
import pathlib
import re
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "rsc/spatium/graph.json"
sys.dont_write_bytecode = True   # importing the matrix script must not leave a __pycache__ in scripts/
sys.path.insert(0, str(ROOT / "scripts"))
import gen_connectivity as gc   # noqa: E402  (what the matrix ran: spaces, scalars, their pairings)

# ── The types ────────────────────────────────────────────────────
# (id, aliases, the C++ type, the matrix name if it is a cell of the matrix)
# `cpp` is what the library spells; for a matrix space it is checked against the matrix's own
# `make<double>()` by a static_assert in the translation unit, so the table cannot lie about it.
def _m(name):
    return f"decltype(connectivity::{name}::make<double>())"

# a spelling and, for the matrix spaces, the cell name
MATRIX_TYPES = [
    ("Euclidean<3>", "E3", "spatium::Euclidean<3>"),
    ("Sphere<2>", "S2", "spatium::Sphere<2>"),
    ("Hyperbolic<2>", "H2", "spatium::Hyperbolic<2>"),
    ("S2xE1", "S2xE1", "spatium::ProductSpace<spatium::Sphere<2>, spatium::Euclidean<1>>"),
    ("SPD<2>", "SPDLogE", "spatium::SPD<2>"),
    ("SPDAffineInvariant<2>", "SPDAff", "spatium::SPDAffineInvariant<2>"),
    ("TorusChart", "TorusChart", "spatium::ParametricSurface<double>"),
    ("SphereLevelSet", "SphereLevelSet", "spatium::ImplicitSurface<double>"),
    ("CylinderADL", "CylinderADL", "third_party::Cyl<double>"),
    ("S2xH2", "S2xH2", None),
    ("H2xE3", "H2xE3", None),
    ("S2xS2", "S2xS2", None),
    ("SPDLogExS2", "SPDLogExS2", None),
    ("SPDAffxS2", "SPDAffxS2", None),
    ("CylxE3", "CylxE3", None),
    ("S2Derived", "S2Derived", None),
    ("H2Derived", "H2Derived", None),
    ("H2DerivedxS2", "H2DerivedxS2", None),
    ("Euclidean<1>", "E1", "spatium::Euclidean<1>"),
    ("Euclidean<4>", "E4", "spatium::Euclidean<4>"),
    ("Euclidean<8>", "E8", "spatium::Euclidean<8>"),
    ("Sphere<1>", "S1", "spatium::Sphere<1>"),
    ("Sphere<3>", "S3", "spatium::Sphere<3>"),
    ("Sphere<4>", "S4", "spatium::Sphere<4>"),
    ("Sphere<7>", "S7", "spatium::Sphere<7>"),
    ("Hyperbolic<1>", "H1", "spatium::Hyperbolic<1>"),
    ("Hyperbolic<4>", "H4", "spatium::Hyperbolic<4>"),
    ("Hyperbolic<8>", "H8", "spatium::Hyperbolic<8>"),
    ("TorusTyped", "TorusTyped", None),
    ("SphereLevelSetTyped", "SphereLevelSetTyped", None),
]
# the readable ids above are the primary names; the cell names are aliases (an id equal to its cell
# name has none). `cpp` for a type with no short spelling says which cell it is.
assert [c for _, c, _ in MATRIX_TYPES] == gc.SPACES, "the matrix has spaces this table does not know"

# (id, the C++ type, the header it needs)
OTHER_TYPES = [
    ("Minkowski<4>", "spatium::Minkowski<4>", "spaces/minkowski.hpp"),
    ("Line<3>", "spatium::geometry::Line<3>", "geometry/geometry.hpp"),
    ("Ray<3>", "spatium::geometry::Ray<3>", "geometry/geometry.hpp"),
    ("Segment<3>", "spatium::geometry::Segment<3>", "geometry/geometry.hpp"),
    ("Hyperplane<3>", "spatium::geometry::Hyperplane<3>", "geometry/geometry.hpp"),
    ("Triangle<3>", "spatium::geometry::Triangle<3>", "geometry/geometry.hpp"),
    ("Polygon<3>", "spatium::geometry::Polygon<3>", "geometry/geometry.hpp"),
    ("Circle<3>", "spatium::geometry::Circle<3>", "geometry/geometry.hpp"),
    ("Disk<3>", "spatium::geometry::Disk<3>", "geometry/geometry.hpp"),
    ("Box<3>", "spatium::geometry::Box<3>", "geometry/geometry.hpp"),
    ("Simplex<3,2>", "spatium::geometry::Simplex<3, 2>", "geometry/geometry.hpp"),
    ("Quadric<double>", "spatium::geometry::Quadric<double>", "geometry/geometry.hpp"),
    ("SO3", "spatium::algebra::SO3<>", "algebra/groups/so3.hpp"),
    ("SE3", "spatium::algebra::SE3<>", "algebra/groups/se3.hpp"),
    ("LieGroupManifold<SO3>", "spatium::LieGroupManifold<spatium::algebra::SO3<>>", "spaces/lie_group.hpp"),
    ("LieGroupManifold<SE3>", "spatium::LieGroupManifold<spatium::algebra::SE3<>>", "spaces/lie_group.hpp"),
    ("FiniteSet<int>", "spatium::discrete::FiniteSet<int>", "discrete/finite_set.hpp"),
]

# The factors of the products in the matrix (checked against the type by a static_assert), the group of
# a Lie group manifold, and the scalars a space is lifted to. Lift names are the node suffixes.
FACTORS = {
    "S2xE1": ("Sphere<2>", "Euclidean<1>"), "S2xH2": ("Sphere<2>", "Hyperbolic<2>"),
    "H2xE3": ("Hyperbolic<2>", "Euclidean<3>"), "S2xS2": ("Sphere<2>", "Sphere<2>"),
    "SPDLogExS2": ("SPD<2>", "Sphere<2>"), "SPDAffxS2": ("SPDAffineInvariant<2>", "Sphere<2>"),
    "CylxE3": ("CylinderADL", "Euclidean<3>"), "H2DerivedxS2": ("H2Derived", "Sphere<2>"),
}
GROUPS = {"LieGroupManifold<SO3>": "SO3", "LieGroupManifold<SE3>": "SE3"}
LIFTS = [("Real50", "spatium::Real50", True, "Real50"), ("Dual", "spatium::Dual<double>", False, "Dual<double>")]   # suffix, C++, needs Boost, scalar node
# the samplers of the types that are not cells of the matrix (tests/symmetry/rsc_pair_support.hpp)
OTHER_CELLS = {"LieGroupManifold<SO3>": "rsc_pair::LieSO3", "LieGroupManifold<SE3>": "rsc_pair::LieSE3"}
# the quantities a route can be emitted for, and so lifted; frechet_mean and verify_exp_log are algorithms over a set of points
EMITTABLE = ["exp_map", "log_map", "metric", "project", "normal", "distance", "norm_at", "geodesic", "midpoint"]
COMPONENT_OPS = {"exp_map": "Exponential", "log_map": "Logarithmic", "metric": "Metrized", "distance": "Distanced"}
GROUP_OPS = ["compose", "inverse", "lie_exp", "lie_log"]   # operations of a group node, each the member of that name (exp, log for the Lie ones)

# (id, the C++ type, the matrix's name for it, needs Boost)
SCALARS = [
    ("double", "double", "double", False),
    ("float", "float", "float", False),
    ("Real50", "spatium::Real50", "Real50", True),
    ("Dual<double>", "spatium::Dual<double>", "Dual", False),
]

# ── The operations ───────────────────────────────────────────────
# quantity -> the matrix probes that exercise an edge yielding it (by variant: None = any)
PROBES = {
    ("exp_map", None): ["ExpLog"],
    ("log_map", None): ["ExpLog"],
    ("metric", None): ["DerivedDistance"],
    ("distance", "member"): ["MetricAxioms"],
    ("distance", "adl"): ["MetricAxioms"],
    ("distance", "derived"): ["DerivedDistance"],
    ("norm_at", None): ["DerivedDistance"],
    ("geodesic", None): ["Midpoint"],
    ("midpoint", None): ["Midpoint"],
    ("frechet_mean", None): ["FrechetMean"],
    ("verify_exp_log", None): ["VerifyExpLog"],
    # project and normal: no cell exercises them, so they stay at level 0
}
# quantities a type provides itself (or by ADL): the member / ADL detection in the translation unit
PROVIDED = ["exp_map", "log_map", "metric", "project", "normal", "distance"]
# derived and generic: quantity, variant, concepts that make it available, quantities it is built from
BUILT = [
    ("distance", "derived", ["Logarithmic", "Metrized"], ["log_map", "metric"]),
    ("norm_at", "derived", ["Metrized"], ["metric"]),
    ("geodesic", "generic", ["Exponential", "Logarithmic"], ["exp_map", "log_map"]),
    ("midpoint", "generic", ["Exponential", "Logarithmic"], ["geodesic"]),
    ("frechet_mean", "generic", ["Riemannian"], ["exp_map", "log_map", "metric"]),
    ("verify_exp_log", "generic", ["Exponential", "Logarithmic"], ["exp_map", "log_map"]),
]
OP_ORDER = ["exp_map", "log_map", "metric", "project", "normal", "distance", "norm_at", "geodesic", "midpoint",
            "frechet_mean", "verify_exp_log"]
VARIANT_ORDER = ["member", "adl", "derived", "generic"]

# ── The concepts (parsed) ────────────────────────────────────────

def strip_source(text):
    text = re.sub(r"/\*.*?\*/", lambda m: "\n" * m.group(0).count("\n"), text, flags=re.S)
    text = re.sub(r"//[^\n]*", "", text)
    text = re.sub(r'"(?:\\.|[^"\\\n])*"', '""', text)
    return text


def balanced(text, start, open_ch, close_ch):
    """Index just past the bracket that closes the one at `start`."""
    depth = 0
    for i in range(start, len(text)):
        if text[i] == open_ch:
            depth += 1
        elif text[i] == close_ch:
            depth -= 1
            if depth == 0:
                return i + 1
    raise ValueError("unbalanced")


def split_top(s):
    parts, depth, cur = [], 0, ""
    for ch in s:
        if ch in "<([{":
            depth += 1
        elif ch in ">)]}":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append(cur)
            cur = ""
        else:
            cur += ch
    if cur.strip():
        parts.append(cur)
    return [p.strip() for p in parts]


def parse_concepts():
    """[{name, ns, header, params, first, body}] for every `concept` under include/spatium/."""
    found = {}
    for path in sorted((ROOT / "include/spatium").rglob("*.hpp")):
        rel = path.relative_to(ROOT / "include/spatium").as_posix()
        if rel.startswith("vendor/"):
            continue
        text = strip_source(path.read_text())
        stack = []
        for m in re.finditer(r"\bnamespace\s+((?:inline\s+)?\w+(?:\s*::\s*(?:inline\s+)?\w+)*)\s*\{|\{|\}|\bconcept\s+(\w+)\s*=", text):
            tok = m.group(0)
            if m.group(1):
                stack.append(("ns", re.sub(r"\binline\s+|\s+", "", m.group(1))))
            elif tok == "{":
                stack.append(("blk", ""))
            elif tok == "}":
                if stack:
                    stack.pop()
            else:
                name = m.group(2)
                ns = "::".join(n for k, n in stack if k == "ns")
                head = text.rfind("template", 0, m.start())
                lt = text.find("<", head)
                gt = balanced(text, lt, "<", ">")
                if text[gt:m.start()].strip():
                    raise ValueError(f"{rel}: concept {name}: something between the template header and the name")
                params = split_top(text[lt + 1:gt - 1])
                end = m.end()
                depth, i = 0, end
                while i < len(text) and not (text[i] == ";" and depth == 0):
                    depth += text[i] in "({[" and 1 or (text[i] in ")}]" and -1 or 0)
                    i += 1
                key = f"{ns}::{name}"
                found.setdefault(key, {"name": name, "ns": ns, "header": rel, "params": params, "body": text[end:i]})
    return list(found.values())


def param_name(decl):
    """The identifier of a template parameter declaration, and the constraint if it is `Concept T`."""
    decl = decl.split("=")[0].strip()
    m = re.match(r"^(?:(?:class|typename)\s+)?(?:([\w:]+)\s+)?(\w+)$", decl)
    if m and m.group(1) in (None, "class", "typename"):
        return m.group(2), None
    if m:
        return m.group(2), m.group(1)
    return None, None


def split_at(expr, op):
    """expr cut at the operator `op` that is not inside any bracket."""
    parts, depth, cur, i = [], 0, "", 0
    while i < len(expr):
        ch = expr[i]
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        if depth == 0 and expr.startswith(op, i):
            parts.append(cur)
            cur = ""
            i += len(op)
            continue
        cur += ch
        i += 1
    parts.append(cur)
    return [p.strip() for p in parts]


def conjunct_refs(expr, pname):
    """(written namespace, name) of the concepts applied to `pname` that the expression REQUIRES: the
    terms of its top-level conjunction, through parentheses that hold no disjunction. `A<S> || B<S>`
    refines neither."""
    expr = expr.strip()
    if split_at(expr, "||") != [expr]:
        return []
    refs = []
    for term in split_at(expr, "&&"):
        if term.startswith("(") and balanced(term, 0, "(", ")") == len(term):
            refs += conjunct_refs(term[1:-1], pname)
            continue
        m = re.fullmatch(r"((?:\w+::)*)(\w+)\s*<\s*" + re.escape(pname) + r"\s*(?:,[^;]*)?>", term)
        if m:
            refs.append((m.group(1), m.group(2)))
            continue
        for n in re.finditer(r"\brequires\s+((?:\w+::)*)(\w+)\s*<\s*" + re.escape(pname) + r"\s*(?:,[^;]*)?>\s*;", term):
            refs.append((n.group(1), n.group(2)))
    return refs


def concept_graph(concepts):
    """Node ids, which concepts are unary, and the refinement edges."""
    by_name = {}
    for c in concepts:
        by_name.setdefault(c["name"], []).append(c)
    for c in concepts:
        dup = len(by_name[c["name"]]) > 1
        short = c["ns"][len("spatium::"):] if c["ns"].startswith("spatium::") else ""
        c["id"] = f"{short}::{c['name']}" if dup and short else c["name"]
        c["cpp"] = f"{c['ns']}::{c['name']}"
        first = c["params"][0] if c["params"] else ""
        starts_type = re.match(r"^(class|typename)\b", first) is not None
        constrained = param_name(first)[1] if first else None
        rest_default = all("=" in p for p in c["params"][1:])
        c["unary"] = bool(c["params"]) and (starts_type or constrained is not None) and rest_default
    ids = {c["cpp"]: c for c in concepts}

    def resolve(name, scope_ns, written_ns):
        cands = by_name.get(name, [])
        if len(cands) <= 1:
            return cands[0] if cands else None
        # several of that name: the one in the same namespace, else the nearest enclosing
        for c in cands:
            if c["ns"] == scope_ns:
                return c
        parts = scope_ns.split("::")
        while parts:
            parts.pop()
            for c in cands:
                if c["ns"] == "::".join(parts):
                    return c
        return cands[0]

    refines = set()
    for c in concepts:
        if not c["params"]:
            continue
        pname, constraint = param_name(c["params"][0])
        if pname is None:
            continue
        if constraint:
            t = resolve(constraint.split("::")[-1], c["ns"], constraint)
            if t and t is not c:
                refines.add((c["id"], t["id"]))
        for written, name in conjunct_refs(c["body"], pname):
            t = resolve(name, c["ns"], written)
            if t and t is not c:
                refines.add((c["id"], t["id"]))
    return sorted(refines)


# ── Compiling the facts ──────────────────────────────────────────

def compile_flags(boost):
    return ["-std=c++23", "-O0", "-w", f"-I{ROOT}/include", f"-I{ROOT}/tests/connectivity",
            f"-DSPATIUM_HAS_BOOST_MULTIPRECISION={1 if boost else 0}"]


PROBE_PRELUDE = r"""
namespace rsc_probe {
// The poison pills of core/access.hpp: an unqualified call finds a user's overload by ADL and nothing else.
void exp_map() = delete;
void log_map() = delete;
void metric_at() = delete;
void distance() = delete;
void project() = delete;
void normal() = delete;

template<class S> concept IsSpace = requires {
    typename spatium::spaces::point_t<S>;
    typename spatium::spaces::tangent_t<S>;
    typename spatium::spaces::scalar_t<S>;
};
template<class S> using P = spatium::spaces::point_t<S>;
template<class S> using V = spatium::spaces::tangent_t<S>;
template<class S> using T = spatium::spaces::scalar_t<S>;

template<class G> concept IsGroup = requires { typename G::ElementType; typename G::AlgebraType; };
template<class G> concept GroupCompose = IsGroup<G> && requires(const G& g, const typename G::ElementType& a) { g.compose(a, a); };
template<class G> concept GroupInverse = IsGroup<G> && requires(const G& g, const typename G::ElementType& a) { g.inverse(a); };
template<class G> concept GroupExp = IsGroup<G> && requires(const G& g, const typename G::AlgebraType& v) { g.exp(v); };
template<class G> concept GroupLog = IsGroup<G> && requires(const G& g, const typename G::ElementType& a) { g.log(a); };

template<class S> concept MemExp = IsSpace<S> && requires(const S& s, const P<S>& p, const V<S>& v, T<S> t) { s.exp_map(p, v, t); };
template<class S> concept AdlExp = IsSpace<S> && requires(const S& s, const P<S>& p, const V<S>& v, T<S> t) { exp_map(s, p, v, t); };
template<class S> concept MemLog = IsSpace<S> && requires(const S& s, const P<S>& p, const P<S>& q) { s.log_map(p, q); };
template<class S> concept AdlLog = IsSpace<S> && requires(const S& s, const P<S>& p, const P<S>& q) { log_map(s, p, q); };
template<class S> concept MemMetric = IsSpace<S> && requires(const S& s, const P<S>& p, const V<S>& u, const V<S>& v) { s.metric_at(p, u, v); };
template<class S> concept AdlMetric = IsSpace<S> && requires(const S& s, const P<S>& p, const V<S>& u, const V<S>& v) { metric_at(s, p, u, v); };
template<class S> concept MemProject = IsSpace<S> && requires(const S& s, const P<S>& p) { s.project(p); };
template<class S> concept AdlProject = IsSpace<S> && requires(const S& s, const P<S>& p) { project(s, p); };
template<class S> concept MemNormal = IsSpace<S> && requires(const S& s, const P<S>& p) { s.normal(p); };
template<class S> concept AdlNormal = IsSpace<S> && requires(const S& s, const P<S>& p) { normal(s, p); };
template<class S> concept MemDistance = IsSpace<S> && requires(const S& s, const P<S>& p, const P<S>& q) { s.distance(p, q); };
template<class S> concept AdlDistance = IsSpace<S> && requires(const S& s, const P<S>& p, const P<S>& q) { distance(s, p, q); };
}
"""

MEMBER_ADL = {   # quantity -> (member concept, ADL concept)
    "exp_map": ("MemExp", "AdlExp"), "log_map": ("MemLog", "AdlLog"), "metric": ("MemMetric", "AdlMetric"),
    "project": ("MemProject", "AdlProject"), "normal": ("MemNormal", "AdlNormal"),
    "distance": ("MemDistance", "AdlDistance"),
}


def all_header_source(concepts, boost):
    heads = sorted({c["header"] for c in concepts if c["unary"]} | {t[2] for t in OTHER_TYPES}
                   | {"core/access.hpp", "algebra/dual.hpp", "spaces/product.hpp", "algebra/calculus.hpp",
                      "core/verify.hpp", "spaces/euclidean.hpp", "spaces/sphere.hpp", "spaces/hyperbolic.hpp"})
    lines = ["#pragma once"]
    if boost:
        lines.append("#include <spatium/core/precision.hpp>")
    lines += [f"#include <spatium/{h}>" for h in heads]
    lines += ['#include "probes.hpp"', "#include <cstdio>", "#include <type_traits>"]
    return "\n".join(lines) + "\n"


def tu_source(tid, cpp_expr, check_cpp, concepts, asserts=()):
    out = ['#include "rsc_all.hpp"', PROBE_PRELUDE, f"using Self = {cpp_expr};"]
    if check_cpp:
        out.append(f"static_assert(std::is_same_v<Self, {check_cpp}>, \"the type table disagrees with the matrix cell\");")
    out += list(asserts)
    out.append("int main() {")
    out.append("    using namespace rsc_probe;")
    for c in concepts:
        if c["unary"]:
            out.append(f'    std::printf("sat|{c["id"]}|%d\\n", int({c["cpp"]}<Self>));')
    for q, (mem, adl) in MEMBER_ADL.items():
        out.append(f'    std::printf("has|{q}|member|%d\\n", int({mem}<Self>));')
        out.append(f'    std::printf("has|{q}|adl|%d\\n", int({adl}<Self>));')
    for q, concept in (("compose", "GroupCompose"), ("inverse", "GroupInverse"), ("lie_exp", "GroupExp"), ("lie_log", "GroupLog")):
        out.append(f'    std::printf("has|{q}|member|%d\\n", int({concept}<Self>));')
    out.append("}")
    return "\n".join(out) + "\n"


def tu_run(work, tid, cpp_expr, check_cpp, concepts, boost, asserts=()):
    stem = re.sub(r"\W", "_", tid)
    src, exe = work / f"{stem}.cpp", work / stem
    src.write_text(tu_source(tid, cpp_expr, check_cpp, concepts, asserts))
    c = subprocess.run(["g++", *compile_flags(boost), str(src), "-o", str(exe)], capture_output=True, text=True,
                       cwd=work)
    if c.returncode != 0:
        errs = [l for l in c.stderr.splitlines() if "error" in l][:3]
        raise SystemExit(f"gen_graph: {tid}: the translation unit does not compile\n  " + "\n  ".join(errs))
    r = subprocess.run([str(exe)], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"gen_graph: {tid}: the translation unit failed to run ({r.returncode})")
    sat, has = {}, {}
    for line in r.stdout.splitlines():
        p = line.split("|")
        if p[0] == "sat":
            sat[p[1]] = p[2] == "1"
        elif p[0] == "has":
            has[(p[1], p[2])] = p[3] == "1"
    return tid, sat, has


def precompile(work, concepts, boost):
    (work / "rsc_all.hpp").write_text(all_header_source(concepts, boost))
    (work / "probes.hpp").write_text((ROOT / "tests/connectivity/probes.hpp").read_text())
    cmd = ["g++", *compile_flags(boost), "-x", "c++-header", str(work / "rsc_all.hpp"), "-o", str(work / "rsc_all.hpp.gch")]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=work)
    if r.returncode != 0:
        errs = [l for l in r.stderr.splitlines() if "error" in l][:5]
        raise SystemExit("gen_graph: the headers do not compile together\n  " + "\n  ".join(errs))


# ── The matrix's levels ──────────────────────────────────────────

def matrix_cells():
    text = (ROOT / "tests/connectivity/expected.hpp").read_text()
    cells = {}
    for m in re.finditer(r'X\(([^,]+(?:<[^)]*>)?),\s*"([^"]+)",\s*(\w+),\s*(\w+),\s*(\d)\)', text):
        cells[(m.group(2), m.group(3), m.group(4))] = int(m.group(5))
    if not cells:
        raise SystemExit("gen_graph: no cells read from tests/connectivity/expected.hpp")
    return cells


LEVEL_KEY = {"double": "double", "float": "float", "Real50": "Real50", "Dual": "Dual<double>"}   # the matrix's name -> the scalar node


def levels_for(cells, matrix_name, quantity, variant, boost):
    """(level on double, {scalar: level}) of an edge on a matrix space; level 0 with no scalars if no probe."""
    probes = PROBES.get((quantity, variant)) or PROBES.get((quantity, None))
    if not probes or matrix_name is None:
        return 0, {}
    ran = gc.SPACE_SCALARS.get(matrix_name)
    per = {}
    for s in gc.SCALARS:
        if s[0] not in LEVEL_KEY or (s[3] and not boost):
            continue
        if ran is not None and s[0] not in ran:
            continue
        per[LEVEL_KEY[s[0]]] = min(cells.get((s[0], matrix_name, p), 0) for p in probes)
    return per.get("double", 0), per


# ── The graph ────────────────────────────────────────────────────

# The erased surfaces' exp, log and distance integrate a geodesic flow: over Real50 one call took 0.24 s
# (measured 2026-10-20, rsc_find on TorusChart distance), so a route times the inputs of a pair times the
# operations is minutes per type. Their Real50 lifts are left out; Dual stays.
SLOW_REAL50 = {"TorusChart", "SphereLevelSet"}


def lifted(cell, scalar_suffix, boost):
    """Is the space of this cell lifted to this scalar? Where the matrix ran it, and not where it costs minutes."""
    if cell is None:
        return False
    if scalar_suffix == "Real50":
        return boost and cell not in gc.HEAVY_SPACES | SLOW_REAL50 and "Real50" in (gc.SPACE_SCALARS.get(cell) or {"Real50"})
    return "Dual" in (gc.SPACE_SCALARS.get(cell) or {"Dual"})


def build(work, jobs, boost):
    concepts = parse_concepts()
    refines = concept_graph(concepts)
    precompile(work, concepts, boost)
    cells = matrix_cells()
    cell_of = {tid: cell for tid, cell, _ in MATRIX_TYPES}

    # (id, aliases, cpp for the TU, cpp to show, check, matrix name, kind, attributes, asserts)
    types = []
    for tid, cell, spelling in MATRIX_TYPES:
        shown = spelling or f"the matrix cell {cell}"
        attrs = {"cell": f"connectivity::{cell}"}
        asserts = []
        if cell in FACTORS:
            attrs["factors"] = list(FACTORS[cell])
            for k, f in enumerate(FACTORS[cell], 1):
                asserts.append(f"static_assert(std::is_same_v<std::remove_cvref_t<decltype(std::declval<Self&>().space{k})>, "
                               f"std::remove_cvref_t<decltype(connectivity::{cell_of[f]}::make<double>())>>, \"the factor table disagrees with the type\");")
        types.append((tid, [] if tid == cell else [cell], _m(cell), shown, spelling, cell, "type", attrs, asserts))
    for tid, cpp, _ in OTHER_TYPES:
        attrs = {}
        asserts = []
        if tid in OTHER_CELLS:
            attrs["cell"] = OTHER_CELLS[tid]
        if tid in GROUPS:
            attrs["group"] = GROUPS[tid]
            gcpp = {t[0]: t[1] for t in OTHER_TYPES}[GROUPS[tid]]
            asserts.append(f"static_assert(std::is_same_v<std::remove_cvref_t<decltype(std::declval<Self&>().group)>, {gcpp}>, \"the group table disagrees with the type\");")
        types.append((tid, [], cpp, cpp, None, None, "type", attrs, asserts))
    for tid, cpp, cell, needs_boost in SCALARS:
        if needs_boost and not boost:
            continue
        types.append((tid, [], cpp, cpp, None, None, "scalar", {}, []))
    base_count = len(types)
    # the same spaces over another scalar
    for tid, cell, spelling in MATRIX_TYPES:
        for suffix, ctype, needs_boost, _ in LIFTS:
            if not lifted(cell, suffix, boost):
                continue
            attrs = {"cell": f"connectivity::{cell}", "lift_of": tid, "scalar": suffix}
            if cell in FACTORS:
                attrs["factors"] = [f"{f}@{suffix}" for f in FACTORS[cell]]
            asserts = []
            types.append((f"{tid}@{suffix}", [], f"decltype(connectivity::{cell}::make<{ctype}>())",
                          f"{spelling or 'the matrix cell ' + cell} over {ctype}", None, cell, "type", attrs, asserts))
    # a product is lifted only with its factors lifted, or its components edges would name nodes that are not there
    ids = {t[0] for t in types}
    types = [t for t in types if "factors" not in t[7] or all(f in ids for f in t[7]["factors"])]

    with cf.ThreadPoolExecutor(jobs) as ex:
        facts = {t[0]: (sat, has) for t, (tid, sat, has) in zip(types, ex.map(
            lambda t: tu_run(work, t[0], t[2], t[4], concepts, boost, t[8]), types))}

    nodes = []
    for c in sorted(concepts, key=lambda c: c["id"]):
        n = {"id": c["id"], "kind": "concept", "cpp": c["cpp"], "header": c["header"], "arity": len(c["params"])}
        if not c["unary"]:
            n["unary"] = False
        nodes.append(n)
    for tid, aliases, _, shown, _, _, kind, attrs, _ in types:
        n = {"id": tid, "kind": kind, "cpp": shown}
        if aliases:
            n["aliases"] = aliases
        n.update(attrs)
        nodes.append(n)

    edges = [{"kind": "refines", "from": a, "to": b} for a, b in refines]
    unary = [c["id"] for c in sorted(concepts, key=lambda c: c["id"]) if c["unary"]]
    ops_of = {}   # type id -> {quantity: [edge]}, the operation edges that stay on the type
    for tid, _, _, _, _, cell, kind, attrs, _ in types:
        sat, has = facts[tid]
        # the same fact two ways: the customization-point concepts of access.hpp against the member / ADL detection
        for q, concept in (("exp_map", "Exponential"), ("log_map", "Logarithmic"), ("metric", "Metrized")):
            if sat[concept] != (has[(q, "member")] or has[(q, "adl")]):
                raise SystemExit(f"gen_graph: {tid}: {concept} says {sat[concept]} but member/ADL detection of {q} disagrees")
        # core/access.hpp withholds the derived distance from a space with an interval and proper time
        derived_ok = sat["Logarithmic"] and sat["Metrized"] and not sat.get("HasInterval", False)
        if sat["Distanced"] != (has[("distance", "member")] or has[("distance", "adl")] or derived_ok):
            raise SystemExit(f"gen_graph: {tid}: Distanced disagrees with member/ADL/derived detection of distance")
        for cid in unary:
            if sat[cid]:
                edges.append({"kind": "satisfies", "from": tid, "to": cid})
        ops = []
        for q in PROVIDED:
            for variant in ("member", "adl"):
                if has[(q, variant)]:
                    ops.append((q, variant, [tid], []))
        for q, variant, needs_concepts, needs in BUILT:
            if q == "distance" and variant == "derived" and sat.get("HasInterval", False):
                continue   # withheld from a space with an interval and proper time (core/access.hpp)
            if all(sat[c] for c in needs_concepts):
                ops.append((q, variant, [tid] + needs_concepts, needs))
        # a group's own operations (the members compose, inverse, exp, log)
        if kind == "type" and tid in GROUPS.values():
            for q in GROUP_OPS:
                if has[(q, "member")]:
                    ops.append((q, "member", [tid], []))
        ops.sort(key=lambda o: (OP_ORDER.index(o[0]) if o[0] in OP_ORDER else len(OP_ORDER) + GROUP_OPS.index(o[0]),
                                VARIANT_ORDER.index(o[1])))
        lift = attrs.get("lift_of")
        scalar_node = next((l[3] for l in LIFTS if l[0] == attrs.get("scalar")), None)
        for q, variant, inputs, needs in ops:
            level, per = levels_for(cells, cell, q, variant, boost)
            if lift:   # the level the matrix reached on this scalar
                level = per.get(scalar_node, 0)
                per = {scalar_node: level} if scalar_node in per else {}
            e = {"kind": "operation", "op": q, "variant": variant, "output": q, "inputs": inputs, "needs": needs,
                 "level": level}
            if per:
                e["levels"] = per
            edges.append(e)
            ops_of.setdefault(tid, {}).setdefault(q, []).append(e)

    # edges between types
    for tid, _, _, _, _, cell, kind, attrs, _ in types:
        sat, has = facts[tid]
        # a space over another scalar: closes a route that is computed there and converted back
        for suffix, _, _, _ in LIFTS:
            lid = f"{tid}@{suffix}"
            for q in EMITTABLE:
                if lid in ops_of and q in ops_of[lid] and tid in ops_of and q in ops_of[tid]:
                    edges.append({"kind": "operation", "op": q, "variant": f"lift:{suffix}", "output": q, "inputs": [tid],
                                  "needs": [f"{lid}#{q}"], "level": max(e["level"] for e in ops_of[lid][q]), "entry": True})
        # a product from its factors
        if "factors" in attrs:
            f1, f2 = attrs["factors"]
            for q, concept in COMPONENT_OPS.items():
                if sat[concept] and q in ops_of.get(f1, {}) and q in ops_of.get(f2, {}):
                    needs = [f"{f1}#{q}"] + ([f"{f2}#{q}"] if f2 != f1 else [])
                    edges.append({"kind": "operation", "op": q, "variant": "components", "output": q, "inputs": [tid],
                                  "needs": needs, "level": 0})
        # a Lie group manifold from its group
        if "group" in attrs:
            g = attrs["group"]
            have = lambda *qs: all(q in ops_of.get(g, {}) for q in qs)
            if have("compose", "lie_exp"):
                edges.append({"kind": "operation", "op": "exp_map", "variant": "group", "output": "exp_map", "inputs": [tid],
                              "needs": [f"{g}#compose", f"{g}#lie_exp"], "level": 0})
            if have("compose", "inverse", "lie_log"):
                edges.append({"kind": "operation", "op": "log_map", "variant": "group", "output": "log_map", "inputs": [tid],
                              "needs": [f"{g}#compose", f"{g}#inverse", f"{g}#lie_log"], "level": 0})
    return {"version": 1, "generator": "rsc/spatium/gen_graph.py", "nodes": nodes, "edges": edges}


def render(graph):
    """One node or edge per line, so a change in the library is a few lines of diff."""
    def block(name, items):
        body = ",\n".join("    " + json.dumps(i, ensure_ascii=False) for i in items)
        return f'  "{name}": [\n{body}\n  ]'
    return ("{\n" + f'  "version": {graph["version"]},\n  "generator": {json.dumps(graph["generator"])},\n'
            + block("nodes", graph["nodes"]) + ",\n" + block("edges", graph["edges"]) + "\n}\n")


def stats(graph):
    nodes, edges = graph["nodes"], graph["edges"]
    count = lambda xs, k, v: sum(1 for x in xs if x[k] == v)
    return (f"concepts {count(nodes, 'kind', 'concept')} "
            f"(unary {sum(1 for n in nodes if n['kind'] == 'concept' and n.get('unary', True))}), "
            f"types {count(nodes, 'kind', 'type')}, scalars {count(nodes, 'kind', 'scalar')}, "
            f"refines {count(edges, 'kind', 'refines')}, satisfies {count(edges, 'kind', 'satisfies')}, "
            f"operations {count(edges, 'kind', 'operation')}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-j", type=int, default=6, help="translation units compiled at once (at most 6)")
    ap.add_argument("--check", action="store_true", help="write nothing; fail if graph.json would change")
    ap.add_argument("--stats", action="store_true", help="print the counts of the committed graph and exit")
    ap.add_argument("--no-boost", action="store_true")
    a = ap.parse_args()
    if a.stats:
        print(stats(json.loads(OUT.read_text())))
        return 0
    with tempfile.TemporaryDirectory() as tmp:
        graph = build(pathlib.Path(tmp), max(1, min(a.j, 6)), not a.no_boost)
    text = render(graph)
    print(stats(graph), file=sys.stderr)
    if a.check:
        if not OUT.exists() or OUT.read_text() != text:
            print(f"stale: {OUT.relative_to(ROOT)} -- rerun rsc/spatium/gen_graph.py", file=sys.stderr)
            if OUT.exists():
                import difflib
                diff = list(difflib.unified_diff(OUT.read_text().splitlines(), text.splitlines(), "committed", "regenerated", lineterm="", n=0))
                print("\n".join(diff[:40]), file=sys.stderr)
            return 1
        return 0
    OUT.write_text(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
