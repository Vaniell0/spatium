#!/usr/bin/env python3
"""rsc_laws: the laws of the library's spaces, found once and then checked column by column.

    python3 rsc/spatium/rsc_laws.py discover [--size N]          enumerate; write laws.json, laws-report.md, tasks.json
    python3 rsc/spatium/rsc_laws.py check [SPACE ...]            the laws of laws.json again on those spaces (all if none)
    python3 rsc/spatium/rsc_laws.py check --changed REV          on the spaces whose headers changed since REV
    python3 rsc/spatium/rsc_laws.py missing                      the arithmetic a space's tangent lacks, now

`discover` runs rsc_laws.cpp: terms over the operations of core/access.hpp, typed by sort and by
degree in length, evaluated on every space of the connectivity matrix on inputs a breaker moves,
merged where equal; every law kept says on which spaces it holds and on which it can be asked.
It is the survey: run when a branch is grown, by a person or an agent, not on every change. Its
product is data -- laws.json -- and what reads it:

  laws-report.md   the classes of spaces the laws draw, each against the library's concepts
                   (rsc/spatium/graph.json): the same set, or a class the library has no name for
  tasks.json       what a law says is wrong: a space that breaks laws nearly all its neighbours
                   keep, and arithmetic a space's tangent does not have

`check` is the regression: the laws of laws.json evaluated again, only on the spaces asked, and any
change of where a law holds or is defined fails it. Its cost is one column per space, so a change
pays for the spaces it touches and the snapshot does not grow with every run.

Run inside `nix develop` from the repository root.
"""
import argparse
import json
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
HERE = ROOT / "rsc" / "spatium"
LAWS = HERE / "laws.json"
REPORT = HERE / "laws-report.md"
TASKS = HERE / "tasks.json"
FLAGS = ["-std=c++23", "-O2", "-w", f"-I{ROOT}/include", f"-I{ROOT}/tests/connectivity"]

# The header each space of the matrix lives in; its #include closure is what a change must touch
# for the space's column to be checked again. tests/connectivity/probes.hpp is in every column.
SPACE_HEADERS = {
    "E1": ["spaces/euclidean.hpp"], "E3": ["spaces/euclidean.hpp"], "E4": ["spaces/euclidean.hpp"],
    "S1": ["spaces/sphere.hpp"], "S2": ["spaces/sphere.hpp"], "S3": ["spaces/sphere.hpp"], "S7": ["spaces/sphere.hpp"],
    "H1": ["spaces/hyperbolic.hpp"], "H2": ["spaces/hyperbolic.hpp"], "H4": ["spaces/hyperbolic.hpp"],
    "S2xE1": ["spaces/product.hpp", "spaces/sphere.hpp", "spaces/euclidean.hpp"],
    "SPD2-logE": ["spaces/spd.hpp"], "SPD2-affine": ["spaces/spd.hpp"],
    "Torus-chart": ["spaces/parametric.hpp"], "Torus-typed": ["spaces/parametric.hpp"],
    "Sphere-implicit": ["spaces/implicit.hpp"], "Sphere-implicit-typed": ["spaces/implicit.hpp"],
    "S2-derived": ["spaces/metric_chart.hpp"], "H2-derived": ["spaces/metric_chart.hpp"],
    "Cylinder-ADL": ["core/access.hpp"],
    "S2xH2": ["spaces/product.hpp", "spaces/sphere.hpp", "spaces/hyperbolic.hpp"],
    "H2xE3": ["spaces/product.hpp", "spaces/hyperbolic.hpp", "spaces/euclidean.hpp"],
    "S2xS2": ["spaces/product.hpp", "spaces/sphere.hpp"],
    "SPD2-affinexS2": ["spaces/product.hpp", "spaces/spd.hpp", "spaces/sphere.hpp"],
    "CylinderxE3": ["spaces/product.hpp", "spaces/euclidean.hpp"],
}
# The matrix's spaces that are nodes of graph.json, where the library's concepts can be read.
GRAPH_TYPES = {
    "E1": "Euclidean<1>", "E3": "Euclidean<3>", "E4": "Euclidean<4>", "S1": "Sphere<1>", "S2": "Sphere<2>",
    "S3": "Sphere<3>", "S7": "Sphere<7>", "H1": "Hyperbolic<1>", "H2": "Hyperbolic<2>", "H4": "Hyperbolic<4>",
    "SPD2-logE": "SPD<2>", "SPD2-affine": "SPDAffineInvariant<2>", "Sphere-implicit": "SphereLevelSet",
    "Sphere-implicit-typed": "SphereLevelSetTyped", "SPD2-affinexS2": "SPDAffxS2",
}
# A space breaking a law is a task when the law holds on at least this share of where it is defined.
NEIGHBOURS = 0.8
# Tasks a person or an agent has looked at and found to be no defect: the space and the reason. A
# triaged task is still listed in the report, apart, and drops back into the queue if its space
# starts breaking a different first law.
TRIAGED = HERE / "laws-triaged.json"


def build(work):
    exe = work / "rsc_laws"
    r = subprocess.run(["g++", *FLAGS, str(HERE / "rsc_laws.cpp"), "-o", str(exe)], capture_output=True, text=True)
    if r.returncode:
        sys.exit("rsc_laws.cpp does not compile:\n" + r.stderr[-4000:])
    return exe


def concepts_by_space():
    g = json.loads((HERE / "graph.json").read_text())
    sat = {}
    for e in g["edges"]:
        if e["kind"] == "satisfies":
            sat.setdefault(e["from"], set()).add(e["to"])
    return {s: sat.get(t, set()) for s, t in GRAPH_TYPES.items()}


def members(bits, spaces):
    return [s for s, b in zip(spaces, bits) if b == "1"]


def classes(data):
    """The laws grouped by the spaces they hold on, largest class first."""
    spaces = data["spaces"]
    groups = {}
    for law in data["laws"]:
        hold = "all" if law["hold"] == "1" * len(spaces) else law["hold"]
        groups.setdefault(hold, []).append(law)
    def size(h):
        return len(spaces) if h == "all" else h.count("1")
    return sorted(groups.items(), key=lambda kv: (-size(kv[0]), -len(kv[1]), kv[0]))


def tasks_of(data):
    """A space that breaks laws its neighbours keep, and arithmetic a space's tangent lacks."""
    spaces = data["spaces"]
    by_space = {}
    for law in data["laws"]:
        if law["hold"] == "all":
            continue
        defined = members(law["defined"], spaces)
        held = set(members(law["hold"], spaces))
        if not defined or len(held) < NEIGHBOURS * len(defined):
            continue
        for s in defined:
            if s not in held:
                by_space.setdefault(s, []).append(law)
    triaged = json.loads(TRIAGED.read_text()) if TRIAGED.exists() else []
    tasks = []
    for s in sorted(by_space):
        laws = sorted(by_space[s], key=lambda l: (l["size"], l["lhs"], l["rhs"]))
        tasks.append({
            "kind": "breaks a law of its neighbours",
            "space": s,
            "claim": f"{laws[0]['lhs']} {laws[0]['rel']} {laws[0]['rhs']}",
            "laws": len(laws),
            "examples": [f"{l['lhs']} {l['rel']} {l['rhs']}" for l in laws[:5]],
            "check": f"python3 rsc/spatium/rsc_laws.py check {s}",
        })
        why = [t["reason"] for t in triaged if t["space"] == s and t["claim"] == tasks[-1]["claim"]]
        if why:
            tasks[-1]["triaged"] = why[0]
    for m in data["missing"]:
        tasks.append({"kind": "missing operation", "space": m["space"], "claim": m["operation"],
                      "check": "python3 rsc/spatium/rsc_laws.py discover"})
    return tasks


def report(data, tasks):
    spaces = data["spaces"]
    concepts = concepts_by_space()
    mapped = [s for s in spaces if s in concepts]
    named = sorted(set().union(*concepts.values()))
    out = ["# Laws of the spaces", "",
           "Written by `rsc/spatium/rsc_laws.py discover`; checked by `rsc_laws.py check`. Terms up to size "
           f"{data['max_size']} over log, exp, the metric and the distance, {data['terms']} kept and {data['merged']} merged, "
           f"on {len(spaces)} spaces of the connectivity matrix, {data['tuples']} inputs each.", "",
           f"{len(data['laws'])} laws in {len(classes(data))} classes. A class is the set of spaces a law holds on; "
           "it is set against the library's concepts on the spaces that are nodes of graph.json "
           f"({', '.join(mapped)}).", "",
           "| spaces | laws | the library's concept on those spaces | shortest law |", "|---|---|---|---|"]
    for hold, laws in classes(data):
        held = spaces if hold == "all" else members(hold, spaces)
        on_mapped = {s for s in held if s in concepts}
        same = [c for c in named if {s for s in mapped if c in concepts[s]} == on_mapped] if hold != "all" else []
        first = min(laws, key=lambda l: (l["size"], l["lhs"]))
        who = "every space" if hold == "all" else ", ".join(held)
        out.append(f"| {who} | {len(laws)} | {', '.join(same) or '--'} | `{first['lhs']} {first['rel']} {first['rhs']}` |")
    open_ = [t for t in tasks if "triaged" not in t]
    done = [t for t in tasks if "triaged" in t]
    out += ["", f"## Tasks ({len(open_)} open, {len(done)} triaged)", ""]
    for t in open_:
        extra = f" ({t['laws']} laws)" if "laws" in t else ""
        out.append(f"- {t['space']}: {t['kind']}{extra}: `{t['claim']}`")
    if done:
        out += ["", "Triaged, no defect:", ""]
        for t in done:
            out.append(f"- {t['space']}: `{t['claim']}` -- {t['triaged']}")
    return "\n".join(out) + "\n"


def discover(size):
    with tempfile.TemporaryDirectory() as work:
        exe = build(pathlib.Path(work))
        r = subprocess.run([str(exe), str(size)], capture_output=True, text=True)
        if r.returncode:
            sys.exit("rsc_laws failed:\n" + r.stderr[-4000:])
    data = json.loads(r.stdout)
    LAWS.write_text(json.dumps(data, indent=1) + "\n")
    tasks = tasks_of(data)
    TASKS.write_text(json.dumps(tasks, indent=1) + "\n")
    REPORT.write_text(report(data, tasks))
    print(f"{len(data['laws'])} laws, {len(classes(data))} classes, "
          f"{sum('triaged' not in t for t in tasks)} open tasks", file=sys.stderr)


def closure(header):
    r = subprocess.run(["g++", "-std=c++23", "-MM", f"-I{ROOT}/include", str(ROOT / "include" / "spatium" / header)],
                       capture_output=True, text=True, check=True)
    deps = r.stdout.replace("\\\n", " ").split(":", 1)[1].split()
    return {str(pathlib.Path(d).resolve().relative_to(ROOT)) for d in deps}


def changed_spaces(rev):
    diff = subprocess.run(["git", "diff", "--name-only", rev], capture_output=True, text=True, check=True, cwd=ROOT)
    files = set(diff.stdout.split())
    if {"tests/connectivity/probes.hpp", "rsc/spatium/rsc_laws.cpp"} & files:
        return list(SPACE_HEADERS)
    hit = []
    for s, headers in SPACE_HEADERS.items():
        if any(closure(h) & files for h in headers):
            hit.append(s)
    return hit


def changes_of(data, asked, rows):
    """Where a law's holding or being defined differs from the snapshot, on the spaces asked.

    `rows` are rsc_laws check's lines: a law's index, then a 0/1 per space asked for whether it
    holds there and whether it is defined there."""
    spaces = data["spaces"]
    column = {s: spaces.index(s) for s in asked}
    changes = []
    for row in rows:
        i, hold, defined = row.split("\t")
        law = data["laws"][int(i)]
        for k, s in enumerate(asked):
            was_h = law["hold"] == "all" or law["hold"][column[s]] == "1"
            was_d = law["defined"] == "all" or law["defined"][column[s]] == "1"
            if law["hold"] == "all" and defined[k] == "0":
                was_d = was_h = False   # a law of every space where it can be asked: undefined is no change
            now_h, now_d = hold[k] == "1", defined[k] == "1"
            if (was_h, was_d) != (now_h, now_d):
                what = "held, now fails" if was_h and not now_h else "failed, now holds" if now_h and not was_h else \
                    "is now defined" if now_d else "is no longer defined"
                changes.append(f"{s}: {law['lhs']} {law['rel']} {law['rhs']} {what}")
    return changes


def check(only, changed):
    data = json.loads(LAWS.read_text())
    spaces = data["spaces"]
    asked = changed_spaces(changed) if changed else (only or spaces)
    unknown = [s for s in asked if s not in spaces]
    if unknown:
        sys.exit(f"not a space of the matrix: {', '.join(unknown)}")
    if not asked:
        print("no space's headers changed: nothing to check")
        return 0
    asked = sorted(set(asked), key=spaces.index)   # the order rsc_laws prints its columns in
    lines = [f"{i}\t{l['lhs_form']}\t{l['rel']}\t{l['rhs_form']}" for i, l in enumerate(data["laws"])]
    with tempfile.TemporaryDirectory() as work:
        exe = build(pathlib.Path(work))
        r = subprocess.run([str(exe), "check", *asked], input="\n".join(lines) + "\n", capture_output=True, text=True)
        if r.returncode:
            sys.exit("rsc_laws check failed:\n" + r.stderr[-4000:])
    changes = changes_of(data, asked, r.stdout.splitlines())
    print(f"{len(data['laws'])} laws checked on {len(asked)} spaces ({', '.join(asked)}): {len(changes)} changes")
    for c in changes:
        print("  " + c)
    if changes:
        print("A law that now fails is a regression; one that now holds is a gain -- run `rsc_laws.py discover` to take it.")
    return 1 if changes else 0


def missing():
    with tempfile.TemporaryDirectory() as work:
        exe = build(pathlib.Path(work))
        r = subprocess.run([str(exe), "missing"], capture_output=True, text=True)
        if r.returncode:
            sys.exit("rsc_laws missing failed:\n" + r.stderr[-4000:])
    rows = [dict(zip(("space", "operation"), line.split("\t"))) for line in r.stdout.splitlines()]
    print(json.dumps(rows))
    return 0


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="mode", required=True)
    d = sub.add_parser("discover")
    d.add_argument("--size", type=int, default=7)
    c = sub.add_parser("check")
    c.add_argument("spaces", nargs="*")
    c.add_argument("--changed", metavar="REV")
    sub.add_parser("missing")
    a = ap.parse_args()
    if a.mode == "missing":
        return missing()
    if a.mode == "discover":
        discover(a.size)
        return 0
    return check(a.spaces, a.changed)


if __name__ == "__main__":
    sys.exit(main())
