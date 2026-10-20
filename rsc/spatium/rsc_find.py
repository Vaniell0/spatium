#!/usr/bin/env python3
"""rsc_find: synthesis and certification over the typed graph of the library.

    python3 rsc/spatium/rsc_find.py <type> <operation> [--samples N]    one (type, operation), a table
    python3 rsc/spatium/rsc_find.py --all [--check]                     every (type, operation) with >= 2 routes

Routes are found in rsc/spatium/graph.json by rsc_route (the engine in rsc/core), turned into C++ by
emit.py (one function per edge, nothing written by hand per route), compiled, run against a reference
on inputs drawn from the connectivity matrix's sampler for the type, and the outcome of each route is
recorded.

The reference of a (type, operation) is the route that is one member (or ADL) closed form the matrix
verified at L3; where the type has none, the Real50 lift of its shortest route; failing that, the
shortest of the most verified. Each route is held to a bound that is the sum of what the two paths promise
(tests/symmetry/rsc_pair_support.hpp). Outcomes, for a route that is not the reference:

  does not compile    a red cell: the error is listed
  disagrees           a finding: the worst input, the error and the bound
  agrees              the chain is the same as another route at another scalar (a lift), or it shares a
                      node with the reference, or it is independent of it (emit.py: the rule);
                      an independent route that brings a node no earlier witness covered is written as a
                      pair in tests/symmetry/pairs/rsc_<type>_<operation>.cpp
  cheaper             any agreeing route that makes fewer calls than the reference is listed as well

--all writes rsc/spatium/find-report.md and the witness files; with --check it writes nothing and fails
if either would change. The report carries calls (the cost), not times: a time differs between runs and
would make --check fail. The table of a single (type, operation) shows both.

Run inside `nix develop` from the repository root. At most 6 translation units are compiled at once.
"""
import argparse
import concurrent.futures as cf
import difflib
import json
import pathlib
import re
import subprocess
import sys
import tempfile
import time

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
import emit  # noqa: E402

GRAPH = HERE / "graph.json"
REPORT = HERE / "find-report.md"
PAIRS = ROOT / "tests/symmetry/pairs"
COMMAND = "python3 rsc/spatium/rsc_find.py --all"
OPS = list(emit.SIG)          # the order operations are reported in
DEFAULT_SAMPLES = 300         # the rule of three then bounds the unmatched measure by 1 %
DEPTH = 9                     # a product of two derived factors, lifted, is 8 edges
GENERATED = re.compile(r"^rsc_.*\.cpp$")


# ── Routes and the reference ───────────────────────────────────

def build_route_tool(work):
    exe = work / "rsc_route"
    cmd = ["g++", "-std=c++23", "-O1", "-w", f"-I{ROOT}", f'-DRSC_GRAPH_JSON="{GRAPH}"', str(HERE / "rsc_route.cpp"), "-o", str(exe)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("rsc_find: rsc_route does not build\n" + r.stderr[:2000])
    return exe


def all_routes(work):
    exe = build_route_tool(work)
    r = subprocess.run([str(exe), "--all", "--json", "--depth", str(DEPTH)], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("rsc_find: rsc_route failed\n" + r.stderr)
    return json.loads(r.stdout)["routes"]


def choose_reference(g, routes):
    """(index, why). A closed form of the type's own at L3; else the shortest Real50 lift; else the most verified shortest."""
    for k, r in enumerate(routes):
        e = g.edge(r[0])
        if len(r) == 1 and e["variant"] in ("member", "adl") and e["level"] == 3:
            return k, f"{e['variant']} closed form, L3"
    for k, r in enumerate(routes):
        if g.edge(r[-1])["variant"] == "lift:Real50":
            return k, "Real50 lift"
    def key(k):
        return (-min(g.edge(i)["level"] for i in routes[k]), len(routes[k]), k)
    k = min(range(len(routes)), key=key)
    return k, "shortest of the most verified (no L3 closed form, no Real50 lift)"


class Target:
    def __init__(self, g, tid, op, routes):
        self.tid, self.op, self.routes = tid, op, routes
        self.ref, self.why = choose_reference(g, routes)
        self.rows = {}        # route index -> measured row
        self.red = {}         # route index -> first compiler error
        self.outcome = {}     # route index -> (kind, text)
        self.cheaper = set()


def targets_of(g, listing, only=None):
    out, skipped = [], {"not emitted": 0, "no sampler": 0, "single": 0}
    for entry in listing:
        tid, op, routes = entry["type"], entry["op"], entry["routes"]
        if "@" in tid:
            continue
        if only and (tid, op) != only:
            continue
        if op not in emit.SIG:
            if len(routes) >= 2:
                skipped["not emitted"] += 1
            continue
        if len(routes) < 2:
            skipped["single"] += 1
            continue
        if not g.cell_of(tid):
            skipped["no sampler"] += 1
            continue
        out.append(Target(g, tid, op, routes))
    order = {n["id"]: k for k, n in enumerate(g.doc["nodes"])}
    # factors before products, so a product's factors are witnessed first
    out.sort(key=lambda t: ("factors" in g.nodes[t.tid], order[t.tid], OPS.index(t.op)))
    return out, skipped


# ── The translation unit of a type ─────────────────────────────

def tu_source(g, tid, items, samples):
    """items: [(op, routes, ref, [route indices to run])] for one type."""
    cell = g.cell_of(tid)
    out = ['#include "rsc_find_all.hpp"', "#include <cstdio>", "namespace {", f"using Cell = {cell};",
           "const auto space = Cell::make<double>();"]
    for oi, (op, routes, ref, run) in enumerate(items):
        for j in sorted(set(run) | {ref}):
            out.append(f"namespace o{oi}_r{j} {{\n{emit.emit_route(g, routes[j])}\n}}")
    out.append('''void print_row(const char* type, const char* op, int route, const rsc_pair::Row& r) {
    std::printf("ROW\\t%s\\t%s\\t%d\\t%zu\\t%zu\\t%zu\\t%.9g\\t%.9g\\t%.9g\\t%zu\\t%zu\\t%zu\\t%.1f\\t%.1f\\t%s\\n", type, op, route, r.samples,
                r.failures, r.first_failure, r.worst_ratio, r.worst_dist, r.worst_bound, r.worst_index, r.calls, r.ref_calls,
                r.ns, r.ref_ns, r.worst_input.c_str());
    std::fflush(stdout);
}
}  // namespace
int main() {
    const std::size_t n = @N@;'''.replace("@N@", str(samples)))
    for oi, (op, routes, ref, run) in enumerate(items):
        refcall = f"o{oi}_r{ref}::{emit.entry_name(routes[ref])}(space, {emit.HARNESS_ARGS[op]})"
        for j in run:
            rcall = f"o{oi}_r{j}::{emit.entry_name(routes[j])}(space, {emit.HARNESS_ARGS[op]})"
            out.append(f'''    {{
        const auto row = rsc_pair::compare<Cell>(n, [&](const auto& d) {{ return {refcall}; }},
            [&](const auto& d) {{ return {rcall}; }},
            [&](const auto& d) {{ return rsc_pair::bound<Cell>(rsc_pair::Tier::{emit.TIER[op]}, d); }});
        print_row("{tid}", "{op}", {j}, row);
    }}''')
    out.append("}")
    return "\n".join(out) + "\n"


FLAGS = ["-std=c++23", "-O0", "-w", f"-I{ROOT}/include", f"-I{ROOT}/tests/connectivity", f"-I{ROOT}/tests/symmetry",
         "-DSPATIUM_HAS_BOOST_MULTIPRECISION=1"]


def precompile(work):
    (work / "rsc_find_all.hpp").write_text('#pragma once\n#include "rsc_pair_support.hpp"\n')
    r = subprocess.run(["g++", *FLAGS, "-x", "c++-header", str(work / "rsc_find_all.hpp"), "-o", str(work / "rsc_find_all.hpp.gch")],
                       capture_output=True, text=True, cwd=work)
    if r.returncode != 0:
        raise SystemExit("rsc_find: the support header does not compile\n" + "\n".join(l for l in r.stderr.splitlines() if "error" in l)[:3000])


def build_and_run(work, stem, src):
    """(ok, rows | first error line)"""
    path = work / f"{stem}.cpp"
    path.write_text(src)
    c = subprocess.run(["g++", *FLAGS, str(path), "-o", str(work / stem)], capture_output=True, text=True, cwd=work)
    if c.returncode != 0:
        errs = [l for l in c.stderr.splitlines() if "error" in l]
        return False, (errs[0] if errs else c.stderr[:300])[:400]
    r = subprocess.run([str(work / stem)], capture_output=True, text=True, timeout=3600)
    if r.returncode != 0:
        return False, f"crashed ({r.returncode})"
    return True, [l.split("\t") for l in r.stdout.splitlines() if l.startswith("ROW\t")]


def parse_row(f):
    return {"samples": int(f[4]), "failures": int(f[5]), "first": int(f[6]), "ratio": float(f[7]), "dist": float(f[8]),
            "bound": float(f[9]), "index": int(f[10]), "calls": int(f[11]), "ref_calls": int(f[12]), "ns": float(f[13]),
            "ref_ns": float(f[14]), "input": f[15] if len(f) > 15 else ""}


def run_type(work, g, tid, tgts, samples):
    """Fill rows and red cells for the targets of one type; compile the whole type, then by operation, then by route."""
    stem = emit.file_stem(tid)
    started = time.monotonic()

    def items(sel):
        return [(t.op, t.routes, t.ref, [j for j in range(len(t.routes))]) for t in sel]

    def take(sel, rows):
        for f in rows:
            op, j = f[2], int(f[3])
            for t in sel:
                if t.op == op:
                    t.rows[j] = parse_row(f)

    def done(how):
        print(f"rsc_find: {tid}: {sum(len(t.routes) for t in tgts)} routes, {how}, {time.monotonic() - started:.0f} s", file=sys.stderr, flush=True)

    ok, res = build_and_run(work, stem, tu_source(g, tid, items(tgts), samples))
    if ok:
        take(tgts, res)
        done("compiled as one")
        return
    for t in tgts:   # by operation
        ok, res = build_and_run(work, f"{stem}_{t.op}", tu_source(g, tid, items([t]), samples))
        if ok:
            take([t], res)
            continue
        for j in range(len(t.routes)):   # by route
            if j == t.ref:
                continue
            sel = [(t.op, t.routes, t.ref, [j])]
            ok, res = build_and_run(work, f"{stem}_{t.op}_{j}", tu_source(g, tid, sel, samples))
            if ok:
                take([t], res)
            else:
                t.red[j] = res
        ok, res = build_and_run(work, f"{stem}_{t.op}_ref", tu_source(g, tid, [(t.op, t.routes, t.ref, [t.ref])], samples))
        if ok:
            take([t], res)
        else:
            t.red[t.ref] = res
    done("compiled apart")


# ── Outcomes ───────────────────────────────────────────────────

def classify(g, targets):
    """Fill outcome and cheaper for every route; return the witnesses: [(target, route index)]."""
    covered = set()
    witnesses = []
    for t in targets:
        ref_nodes = g.nodes_of(t.routes[t.ref])
        chain = {ref_nodes: t.ref}
        ref_row = t.rows.get(t.ref)
        for j, r in enumerate(t.routes):
            nodes = g.nodes_of(r)
            if j == t.ref:
                t.outcome[j] = ("reference", f"reference: {t.why}")
                continue
            if j in t.red:
                t.outcome[j] = ("red", ("does not run: " if t.red[j].startswith("crashed") else "does not compile: ") + t.red[j])
                continue
            row = t.rows.get(j)
            if row is None:
                t.outcome[j] = ("red", "not run")
                continue
            if row["failures"]:
                t.outcome[j] = ("disagrees", f"DISAGREES on {row['failures']} of {row['samples']} inputs")
                chain.setdefault(nodes, j)
                continue
            if ref_row and row["calls"] < row["ref_calls"]:
                t.cheaper.add(j)
            if nodes in chain:
                other = chain[nodes]
                t.outcome[j] = ("agrees", "agrees; the chain of the reference at another scalar" if other == t.ref
                                else f"agrees; the chain of route {other + 1} at another scalar")
                continue
            chain[nodes] = j
            shared = sorted(ref_nodes & nodes)
            if shared:
                t.outcome[j] = ("agrees", "agrees; shares " + ", ".join(f"{a}.{b}[{c}]" for a, b, c in shared[:3]) + (" ..." if len(shared) > 3 else "") + " with the reference")
                continue
            new = nodes - covered
            if new:
                t.outcome[j] = ("witness", "agrees and independent: witness")
                witnesses.append((t, j))
                covered |= nodes | ref_nodes
            else:
                t.outcome[j] = ("independent", "agrees and independent; every node already witnessed")
    return witnesses


# ── Witness files ──────────────────────────────────────────────

def witness_file(g, tid, op, entries, samples):
    """entries: [(Target, route index)] of one (type, op)."""
    t0 = entries[0][0]
    stem = f"{emit.file_stem(tid)}_{op}"
    ns = f"rsc_{stem}"
    cell = g.cell_of(tid)
    ref = t0.routes[t0.ref]
    out = [f"// generated by rsc_find, do not edit; regenerate: {COMMAND}",
           f"// {tid} {op}. The reference is {g.describe(ref, tid)} ({t0.why}).",
           f"// Each pair is a route that agreed with it on {samples} inputs of the matrix's sampler and shares no node of",
           "// the typed graph with it; the bound is the sum of what the two paths promise (rsc_pair_support.hpp).",
           '#include "../pair.hpp"', '#include "../rsc_pair_support.hpp"', "",
           "namespace symmetry {", "namespace {", f"namespace {ns} {{", "",
           f"using Cell = {cell};", "const auto space = Cell::make<double>();", "",
           f"namespace ref {{\n{emit.emit_route(g, ref)}\n}}  // namespace ref"]
    regs = []
    for t, j in entries:
        r = t.routes[j]
        label = f"rsc: {tid} {op}: {g.short(r, tid)} vs {g.short(ref, tid)}"
        out.append(f"namespace r{j + 1} {{\n{emit.emit_route(g, r)}\n}}  // namespace r{j + 1}")
        out.append(f'''inline Entry entry_r{j + 1}(std::size_t n) {{
    const std::string name = "{label}";
    return {{name, n, [n, name] {{
        return check<rsc_pair::Draw<Cell>>(
            name, Kind::Function, n,
            [](std::size_t i) {{ return rsc_pair::draw<Cell>(i); }},
            [](const rsc_pair::Draw<Cell>& d) {{ return {emit.call("ref", ref, op)}; }},
            [](const rsc_pair::Draw<Cell>& d) {{ return {emit.call(f"r{j + 1}", r, op)}; }},
            [](const auto& a, const auto& b) {{ return rsc_pair::dist(a, b); }},
            [](const rsc_pair::Draw<Cell>& d) {{ return rsc_pair::bound<Cell>(rsc_pair::Tier::{emit.TIER[op]}, d); }});
    }}}};
}}''')
        regs.append(f'const Registrar registered_{ns}_r{j + 1}{{"{ns}_r{j + 1}", {samples}, [](std::size_t n) {{ return {ns}::entry_r{j + 1}(n); }}}};')
    out += ["", f"}}  // namespace {ns}", "}  // namespace", ""] + regs + ["", "}  // namespace symmetry", ""]
    return "\n".join(out)


def witness_files(g, witnesses, samples):
    groups = {}
    for t, j in witnesses:
        groups.setdefault((t.tid, t.op), []).append((t, j))
    return {f"rsc_{emit.file_stem(tid)}_{op}.cpp": witness_file(g, tid, op, es, samples) for (tid, op), es in groups.items()}


# ── The report ─────────────────────────────────────────────────

def fmt(x):
    return f"{x:.2g}"


def report(g, targets, skipped, witnesses, samples, listing):
    total_routes = sum(len(t.routes) for t in targets)
    findings = [(t, j) for t in targets for j in sorted(t.outcome) if t.outcome[j][0] == "disagrees"]
    reds = [(t, j) for t in targets for j in sorted(t.outcome) if t.outcome[j][0] == "red"]
    cheaper = [(t, j) for t in targets for j in sorted(t.cheaper)]
    agreed = sum(1 for t in targets for j in t.outcome if t.outcome[j][0] in ("agrees", "witness", "independent"))
    refs = len(targets)
    pairs = sum(1 for t in targets for j in t.rows)
    clean = sum(1 for t in targets for j, r in t.rows.items() if r["failures"] == 0)
    L = ["# rsc_find: routes, run against a reference, and what came of each", "",
         f"Generated by `{COMMAND}`; do not edit. `--check` regenerates this file and the witness pairs and fails if either changes.", "",
         "A route of the typed graph (`rsc/spatium/graph.json`) is turned into C++ by `rsc/spatium/emit.py`, compiled, and run against a",
         f"reference on {samples} inputs drawn from the connectivity matrix's sampler for the type (seeded from the Halton sequence).",
         "The bound of a pair is the sum of what its two paths promise (`tests/symmetry/rsc_pair_support.hpp`). The cost of a route is the",
         "number of operation calls one evaluation makes (a scalar lift counts none: it converts). Times are not in this file, they differ",
         "from run to run; `rsc_find.py <type> <operation>` prints them.", "",
         "The reference is the type's own closed form (member or ADL) that the matrix verified at L3; where there is none, the Real50 lift",
         "of its shortest route. Two routes are independent when they share no node, a node being (type, operation, variant) with the",
         "scalar of a lift ignored; arithmetic, sqrt, exp, acos and cos are the trusted base (`emit.py` states the rule and what it cannot see).", "",
         "## Totals", "",
         "| | |", "|---|---|",
         f"| (type, operation) pairs with at least two routes | {len(targets)} |",
         f"| routes found in them | {total_routes} |",
         f"| routes run against a reference | {pairs} |",
         f"| agree (including the references' own chains at another scalar) | {agreed} |",
         f"| witnesses added (pairs in `tests/symmetry/pairs/`) | {len(witnesses)} in {len({(t.tid, t.op) for t, _ in witnesses})} files |",
         f"| agree and cheaper than the reference | {len(cheaper)} |",
         f"| disagreements (findings) | {len(findings)} |",
         f"| red cells (do not compile or run) | {len(reds)} |",
         f"| unmatched measure | {clean} pairs with no failure: below 3/{samples} = {3 / samples:.3f} each (95 %); {pairs - clean} with failures |",
         f"| not run: operations with no emitter (`frechet_mean`, `verify_exp_log`: algorithms over a set of points) | {skipped['not emitted']} pairs |",
         f"| not run: no sampler for the type | {skipped['no sampler']} pairs |",
         f"| operations with one route (nothing to compare) | {skipped['single']} |", ""]
    L += ["## Disagreements", ""]
    if findings:
        L += ["| type | operation | route | worst error | bound | ratio | inputs failed | worst input |", "|---|---|---|---|---|---|---|---|"]
        for t, j in findings:
            r = t.rows[j]
            L.append(f"| {t.tid} | {t.op} | {g.describe(t.routes[j], t.tid)} | {r['dist']:.3g} | {r['bound']:.3g} | {fmt(r['ratio'])} | "
                     f"{r['failures']} of {r['samples']} | input {r['index']}: `{r['input']}` |")
    else:
        L.append("None.")
    L += ["", "## Red cells", ""]
    if reds:
        L += ["| type | operation | route | error |", "|---|---|---|---|"]
        for t, j in reds:
            L.append(f"| {t.tid} | {t.op} | {g.describe(t.routes[j], t.tid)} | `{t.red.get(j, 'not run').replace('|', '/')}` |")
    else:
        L.append("None.")
    L += ["", "## Cheaper than the reference", ""]
    if cheaper:
        L += ["| type | operation | route | calls | reference calls |", "|---|---|---|---|---|"]
        for t, j in cheaper:
            r = t.rows[j]
            L.append(f"| {t.tid} | {t.op} | {g.describe(t.routes[j], t.tid)} | {r['calls']} | {r['ref_calls']} |")
    else:
        L.append("None.")
    L += ["", "## Witnesses added", ""]
    if witnesses:
        L += ["| file | route | against |", "|---|---|---|"]
        for t, j in witnesses:
            L.append(f"| `tests/symmetry/pairs/rsc_{emit.file_stem(t.tid)}_{t.op}.cpp` | {g.describe(t.routes[j], t.tid)} | {g.describe(t.routes[t.ref], t.tid)} |")
    else:
        L.append("None.")
    L += ["", "## By type", ""]
    current = None
    for t in sorted(targets, key=lambda t: [n["id"] for n in g.doc["nodes"]].index(t.tid) * 100 + OPS.index(t.op)):
        if t.tid != current:
            current = t.tid
            L += [f"### {t.tid}", ""]
        L += [f"#### {t.op}  (reference: route {t.ref + 1}, {t.why})", "",
              "| # | route | level | calls | worst ratio | outcome |", "|---|---|---|---|---|---|"]
        for j, r in enumerate(t.routes):
            row = t.rows.get(j)
            lvl = min(g.edge(i)["level"] for i in r)
            ratio = fmt(row["ratio"]) if row else "-"
            calls = str(row["calls"]) if row else "-"
            kind, text = t.outcome[j]
            if j in t.cheaper:
                text += "; cheaper"
            L.append(f"| {j + 1} | {g.describe(r, t.tid)} | L{lvl} | {calls} | {ratio} | {text} |")
        L.append("")
    return "\n".join(L).rstrip("\n") + "\n"


# ── The single table ───────────────────────────────────────────

def table(g, t):
    print(f"{t.tid} {t.op}: {len(t.routes)} routes, reference route {t.ref + 1} ({t.why})")
    print(f"{'#':>3} {'level':>5} {'calls':>5} {'ns/eval':>9} {'ref ns':>9} {'worst':>9}  outcome / route")
    for j, r in enumerate(t.routes):
        row = t.rows.get(j)
        lvl = min(g.edge(i)["level"] for i in r)
        kind, text = t.outcome[j]
        extra = "; cheaper" if j in t.cheaper else ""
        if row:
            print(f"{j + 1:>3} {'L' + str(lvl):>5} {row['calls']:>5} {row['ns']:>9.0f} {row['ref_ns']:>9.0f} {fmt(row['ratio']):>9}  {text}{extra}")
        else:
            print(f"{j + 1:>3} {'L' + str(lvl):>5} {'-':>5} {'-':>9} {'-':>9} {'-':>9}  {text}")
        print(f"{'':>45}{g.describe(r, t.tid)}")
        if kind == "disagrees":
            row = t.rows[j]
            print(f"{'':>45}worst error {row['dist']:.3g} against bound {row['bound']:.3g} at input {row['index']}: {row['input']}")


# ── Main ───────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("type", nargs="?")
    ap.add_argument("operation", nargs="?")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--check", action="store_true", help="with --all: write nothing, fail if the report or a witness would change")
    ap.add_argument("--samples", type=int, default=DEFAULT_SAMPLES)
    ap.add_argument("-j", type=int, default=6, help="translation units compiled at once (at most 6)")
    a = ap.parse_args()
    if a.all == bool(a.type):
        ap.error("give a type and an operation, or --all")
    if a.check and not a.all:
        ap.error("--check goes with --all")
    if a.all and a.samples != DEFAULT_SAMPLES:
        ap.error("--all uses the default number of samples, which the witness files record")
    if a.type and not a.operation:
        ap.error("an operation is missing")
    jobs = max(1, min(a.j, 6))

    g = emit.Graph(GRAPH)
    only = None
    if a.type:
        node = next((n for n in g.doc["nodes"] if a.type == n["id"] or a.type in n.get("aliases", [])), None)
        if node is None or node["kind"] != "type":
            raise SystemExit(f"rsc_find: no type named {a.type!r} in the graph")
        only = (node["id"], a.operation)

    with tempfile.TemporaryDirectory() as tmp:
        work = pathlib.Path(tmp)
        listing = all_routes(work)
        targets, skipped = targets_of(g, listing, only)
        if only and not targets:
            raise SystemExit(f"rsc_find: {only[0]} {only[1]} has fewer than two routes, or cannot be emitted")
        precompile(work)
        by_type = {}
        for t in targets:
            by_type.setdefault(t.tid, []).append(t)
        with cf.ThreadPoolExecutor(jobs) as ex:
            list(ex.map(lambda kv: run_type(work, g, kv[0], kv[1], a.samples), by_type.items()))

    witnesses = classify(g, targets)
    if only:
        table(g, targets[0])
        return 0

    files = witness_files(g, witnesses, a.samples)
    text = report(g, targets, skipped, witnesses, a.samples, listing)
    stale = []
    existing = {p.name: p for p in PAIRS.glob("rsc_*.cpp")}
    if a.check:
        if not REPORT.exists() or REPORT.read_text() != text:
            stale.append(str(REPORT.relative_to(ROOT)))
            if REPORT.exists():
                diff = list(difflib.unified_diff(REPORT.read_text().splitlines(), text.splitlines(), "committed", "regenerated", lineterm="", n=0))
                print("\n".join(diff[:60]), file=sys.stderr)
        for name, src in files.items():
            if name not in existing or existing[name].read_text() != src:
                stale.append(f"tests/symmetry/pairs/{name}")
        stale += [f"tests/symmetry/pairs/{n} (no longer generated)" for n in existing if n not in files]
        for s in stale:
            print(f"stale: {s} -- rerun {COMMAND}", file=sys.stderr)
        return 1 if stale else 0
    REPORT.write_text(text)
    for name, src in files.items():
        (PAIRS / name).write_text(src)
    for n, p in existing.items():
        if n not in files:
            p.unlink()
    print(f"{len(targets)} (type, operation) pairs, {sum(len(t.routes) for t in targets)} routes, {len(witnesses)} witnesses in {len(files)} files, "
          f"{sum(1 for t in targets for o in t.outcome.values() if o[0] == 'disagrees')} disagreements, "
          f"{sum(1 for t in targets for o in t.outcome.values() if o[0] == 'red')} red cells", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
