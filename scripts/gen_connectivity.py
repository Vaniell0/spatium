#!/usr/bin/env python3
"""The connectivity matrix: every (scalar, space, probe) cell of
tests/connectivity/probes.hpp built alone, run, and graded.

  L0  compiles       -- a cell is its own translation unit, so a body that
                        does not compile is found, which a concept cannot do
  L1  runs, finite
  L2  its axiom or symmetry holds (the cell decides)
  L3  its numbers agree with the same cell over the reference scalar:
      double against Real50, every other scalar against double

Writes docs/connectivity.md (the table) and tests/connectivity/expected.hpp
(the cells that reached each level, which tests/test_connectivity.cpp holds
them to). With --check, writes nothing and fails if either would change.

Run inside `nix develop` from the repository root:
    python3 scripts/gen_connectivity.py [-j 8] [--check]
"""
import argparse
import concurrent.futures as cf
import math
import os
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parent.parent

# (name in the table, C++ type, the base type's epsilon, needs Boost)
SCALARS = [
    ("double", "double", 2.220446049250313e-16, False),
    ("float", "float", 1.1920929e-07, False),
    ("long double", "long double", 1.0842021724855044e-19, False),
    ("Real50", "spatium::Real50", 1e-50, True),
    ("Dual", "spatium::Dual<double>", 2.220446049250313e-16, False),
    ("Dual2", "spatium::Dual<spatium::Dual<double>>", 2.220446049250313e-16, False),
]
SPACES = ["E3", "S2", "H2", "S2xE1", "SPDLogE", "SPDAff", "TorusChart", "SphereLevelSet", "CylinderADL",
          # products of spaces green on their own: green together by symmetry
          "S2xH2", "H2xE3", "S2xS2", "SPDLogExS2", "SPDAffxS2", "CylxE3",
          # spaces given by their metric alone: exp and log derived
          "S2Derived", "H2Derived", "H2DerivedxS2",
          # the three families at other dimensions (S2, H2, E3 are above)
          "E1", "E4", "E8", "S1", "S3", "S4", "S7", "H1", "H4", "H8",
          # a surface whose function is a type: exact partials, geodesic exp and log
          "TorusTyped"]
# A dimension is a property of a space, not of a scalar: these run over a
# native, a reference and a derivative-carrying scalar, not all six.
DIMENSION_SPACES = {"E1", "E4", "E8", "S1", "S3", "S4", "S7", "H1", "H4", "H8"}
DIMENSION_SCALARS = {"double", "Real50", "Dual"}
# Which scalars a space runs over, where not all six. The typed torus's geodesic
# flow evaluates its function on nested Duals: over fifty digits a cell is
# minutes, and the compile of Dual2's depth is heavier still, so it rides on
# four, and double is held against long double (eps 1e-19) in place of Real50.
SPACE_SCALARS = {sp: DIMENSION_SCALARS for sp in DIMENSION_SPACES}
SPACE_SCALARS["TorusTyped"] = {"double", "float", "long double", "Dual"}
SPACE_REFERENCE = {"TorusTyped": "long double"}   # what double is held against, instead of Real50
PROBES = ["MetricAxioms", "DerivedDistance", "ExpLog", "Midpoint", "FrechetMean", "VerifyExpLog", "Derivative",
          "Infinity"]
REFERENCE = {"double": "Real50"}   # every other scalar is held against double
# Cells a unit test would take minutes to run: a derived space is a shooting
# of geodesic flows, and over fifty digits each of its cells costs one to ten
# minutes. They are held where the matrix is regenerated from scratch (the
# CI job, `--check`), not in tests/test_connectivity.cpp, which defines
# SPATIUM_CONNECTIVITY_HEAVY to run them too.
HEAVY_SPACES = {"S2Derived", "H2Derived", "H2DerivedxS2"}


def is_heavy(sname, space):
    return sname == "Real50" and space in HEAVY_SPACES


def tolerance(eps, scale):
    return math.sqrt(eps) * 64 * (scale + 1)


def cell_source(ctype, space, probe, sname):
    return f"""#include "probes.hpp"
int main() {{
    std::vector<double> sig;
    const int level = connectivity::run_cell<{ctype}, connectivity::{space}, connectivity::{probe}>(sig);
    connectivity::print_cell("{sname}", "{space}", "{probe}", level, sig);
}}
"""


def build_and_run(workdir, sname, ctype, boost, space, probe):
    stem = f"{sname.replace(' ', '_')}__{space}__{probe}"
    src = workdir / f"{stem}.cpp"
    exe = workdir / stem
    src.write_text(cell_source(ctype, space, probe, sname))
    cmd = ["g++", "-std=c++23", "-O1", "-w", f"-I{ROOT}/include", f"-I{ROOT}/tests/connectivity",
           f"-DSPATIUM_HAS_BOOST_MULTIPRECISION={1 if boost else 0}", str(src), "-o", str(exe)]
    try:
        c = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        return (sname, space, probe, -1, [], "compile timeout")
    if c.returncode != 0:
        first = next((l for l in c.stderr.splitlines() if "error" in l), "compile error")
        return (sname, space, probe, -1, [], first.split("error:")[-1].strip()[:140])
    try:
        r = subprocess.run([str(exe)], capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        return (sname, space, probe, 0, [], "run timeout")
    if r.returncode != 0 or "|" not in r.stdout:
        return (sname, space, probe, 0, [], f"run failed ({r.returncode})")
    parts = r.stdout.strip().split("|")
    sig = [float(x) for x in parts[4].split()] if len(parts) > 4 and parts[4].strip() else []
    return (sname, space, probe, int(parts[3]), sig, "")


def grade_l3(results):
    eps = {s[0]: s[2] for s in SCALARS}
    by = {(r[0], r[1], r[2]): r for r in results}
    graded = {}
    for key, (sname, space, probe, level, sig, note) in by.items():
        lvl = level
        ref_name = SPACE_REFERENCE.get(space, REFERENCE.get(sname, "double")) if sname == "double" \
            else REFERENCE.get(sname, "double")
        ref = by.get((ref_name, space, probe))
        if level == 2 and sname != ref_name and ref and ref[3] >= 1 and len(ref[4]) == len(sig) and sig:
            worse = max(eps[sname], eps[ref[0]])
            if all(abs(a - b) <= tolerance(worse, abs(b)) for a, b in zip(sig, ref[4])):
                lvl = 3
        graded[key] = (lvl, note)
    return graded


SYMBOL = {-1: "✗ L0", 0: "L0", 1: "L1", 2: "L2", 3: "L3"}


def render_md(graded, boost):
    names = [s[0] for s in SCALARS if boost or not s[3]]
    lines = [
        "# Connectivity matrix",
        "",
        "Generated by `scripts/gen_connectivity.py` -- do not edit. Every cell is one",
        "generic algorithm (probe) of the library instantiated on one space over one",
        "scalar, in its own translation unit, and graded:",
        "",
        "- **✗ L0** does not compile; **L0** compiles but fails or crashes when run",
        "- **L1** runs, results finite",
        "- **L2** its axiom or symmetry holds within the scalar's tolerance",
        "- **L3** its numbers agree with the same cell over the reference scalar",
        "  (double against Real50, every other scalar against double)",
        "- **·** not run: the spaces at other dimensions (E1, E4, E8, S1, S3, S4, S7,",
        "  H1, H4, H8) ride on double, Real50 and Dual -- a dimension belongs to the",
        "  space, not to the scalar; what depends on it by theorem is held by `tests/test_dimension_theorems.cpp`.",
        "  The typed torus rides on double, float, long double and Dual (a geodesic flow",
        "  of nested Duals over fifty digits is minutes a cell), double held against long double.",
        "",
        "Probes are in `tests/connectivity/probes.hpp`; `tests/test_connectivity.cpp`",
        "holds every cell to the level recorded here.",
        "",
    ]
    total = sum(1 for k in graded if k[0] in names)
    green = sum(1 for k, (l, _) in graded.items() if k[0] in names and l >= 2)
    lines.append(f"**{green} of {total} cells at L2 or above.**")
    lines.append("")
    lines.append("| space | probe | " + " | ".join(names) + " |")
    lines.append("|---|---|" + "---|" * len(names))
    notes = []
    for space in SPACES:
        for probe in PROBES:
            row = []
            for n in names:
                if space in SPACE_SCALARS and n not in SPACE_SCALARS[space]:
                    row.append("·")      # not run: this space rides on fewer scalars
                    continue
                lvl, note = graded.get((n, space, probe), (-1, "missing"))
                cell = SYMBOL[lvl]
                if note:
                    notes.append(f"- {space} / {probe} / {n}: {note}")
                    cell += f" [{len(notes)}]"
                row.append(cell)
            lines.append(f"| {space} | {probe} | " + " | ".join(row) + " |")
    if notes:
        lines += ["", "## Notes", ""]
        lines += [f"{i + 1}. {n[2:]}" for i, n in enumerate(notes)]
    return "\n".join(lines) + "\n"


def render_hpp(graded, boost):
    ctype = {s[0]: s[1] for s in SCALARS}
    lines = [
        "#pragma once",
        "// Generated by scripts/gen_connectivity.py -- do not edit; rerun it.",
        "// The cells that compile, and the level each reached, for",
        "// tests/test_connectivity.cpp. X(scalar type, scalar name, space, probe, level).",
        "#define CONNECTIVITY_CELLS(X) \\",
    ]
    body, heavy = [], []
    for (sname, space, probe), (lvl, _) in sorted(graded.items()):
        if lvl < 1:
            continue
        needs_boost = next(s[3] for s in SCALARS if s[0] == sname)
        line = f"    X({ctype[sname]}, \"{sname}\", {space}, {probe}, {lvl})"
        (heavy if is_heavy(sname, space) else body).append((needs_boost, line))
    lines += [b + " \\" for nb, b in body if not nb]
    lines.append("    CONNECTIVITY_CELLS_BOOST(X)")
    lines.append("")
    lines.append("#if SPATIUM_HAS_BOOST_MULTIPRECISION")
    lines.append("#define CONNECTIVITY_CELLS_BOOST(X) \\")
    lines += [b + " \\" for nb, b in body if nb]
    lines.append("")
    lines.append("#else")
    lines.append("#define CONNECTIVITY_CELLS_BOOST(X)")
    lines.append("#endif")
    lines.append("")
    lines.append("// Minutes each; run only with SPATIUM_CONNECTIVITY_HEAVY (the CI job regenerates them).")
    lines.append("#if SPATIUM_HAS_BOOST_MULTIPRECISION")
    lines.append("#define CONNECTIVITY_CELLS_HEAVY(X) \\")
    lines += [b + " \\" for nb, b in heavy]
    lines.append("")
    lines.append("#else")
    lines.append("#define CONNECTIVITY_CELLS_HEAVY(X)")
    lines.append("#endif")
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-j", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--no-boost", action="store_true")
    ap.add_argument("--only", help="build only the spaces whose name contains this (a trial: writes nothing)")
    a = ap.parse_args()
    boost = not a.no_boost
    scalars = [s for s in SCALARS if boost or not s[3]]
    with tempfile.TemporaryDirectory() as tmp:
        work = pathlib.Path(tmp)
        spaces = [sp for sp in SPACES if not a.only or a.only in sp]
        jobs = [(s[0], s[1], s[3], sp, pr) for s in scalars for sp in spaces for pr in PROBES
                if sp not in SPACE_SCALARS or s[0] in SPACE_SCALARS[sp]]
        with cf.ThreadPoolExecutor(a.j) as ex:
            results = list(ex.map(lambda j: build_and_run(work, *j), jobs))
    graded = grade_l3(results)
    if a.only:
        for (sname, space, probe), (lvl, note) in sorted(graded.items(), key=lambda kv: (kv[0][1], kv[0][2], kv[0][0])):
            print(f"{space:14} {probe:16} {sname:12} {SYMBOL[lvl]:5} {note}")
        return 0
    md, hpp = render_md(graded, boost), render_hpp(graded, boost)
    md_path, hpp_path = ROOT / "docs/connectivity.md", ROOT / "tests/connectivity/expected.hpp"
    if a.check:
        stale = [p for p, text in ((md_path, md), (hpp_path, hpp)) if not p.exists() or p.read_text() != text]
        for p in stale:
            print(f"stale: {p.relative_to(ROOT)} -- rerun scripts/gen_connectivity.py", file=sys.stderr)
        return 1 if stale else 0
    md_path.write_text(md)
    hpp_path.write_text(hpp)
    print(md.split("\n\n")[5] if len(md.split("\n\n")) > 5 else "")
    return 0


if __name__ == "__main__":
    sys.exit(main())
