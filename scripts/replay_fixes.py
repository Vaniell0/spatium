#!/usr/bin/env python3
"""Replay the history as breaks a program can run again.

    replay_fixes.py <work dir> [--rev main] [--only <sha>,<sha>] [-j N] [--keep]
    replay_fixes.py <work dir> --summary

A commit that changes the library (include/, rsc/, src/) together with a test file tests/test_*.cpp is
replayed twice: its test files are built alone against the library as the commit left it (after) and
against the library of its parent (before), each linked with Catch2's main and run. A test case that
fails before and passes after is a break with its check: the parent is the broken state, the case is
what tells it from the repaired one, and the lines the commit replaced are where it lived.

One row per commit goes to <work dir>/rows.jsonl; a commit already there is skipped. A test file that
does not build before reads an interface the parent lacks (a new capability, not a repair); one that
does not build or pass after is the harness's failure and is counted apart, never as a break.

Run inside `nix develop`: the compiler, Catch2, Boost and Eigen come from its environment.
"""
import collections
import concurrent.futures as cf
import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

REPO = pathlib.Path(__file__).resolve().parent.parent
LIBRARY = ("include/", "rsc/", "src/")
TEST = re.compile(r"^tests/test_[^/]*\.cpp$")
SKIP_TOP = {"gallery"}                                   # pictures, read by no test
BUILD_S, RUN_S = 1200, 300


def git(*args, **kw):
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True, check=True, **kw).stdout


def candidates(rev):
    """Commits of rev, oldest first, that change the library and a test file at once."""
    out = []
    for c in git("rev-list", "--no-merges", "--reverse", rev).split():
        parents = git("rev-list", "--parents", "-n", "1", c).split()[1:]
        if len(parents) != 1:
            continue
        changed = [l.split("\t") for l in git("diff-tree", "--no-commit-id", "--name-status", "-r", c).splitlines()]
        tests = [f[-1] for f in changed if f[0][0] in "AM" and TEST.match(f[-1])]
        if tests and any(f[-1].startswith(LIBRARY) for f in changed):
            out.append((c, parents[0], tests))
    return out


def places(parent, commit):
    """Library lines the commit replaced, in the parent's numbering: {file: [[first line, count], ...]}."""
    got = collections.defaultdict(list)
    name = None
    for l in git("diff", "-U0", parent, commit, "--", *[p.rstrip("/") for p in LIBRARY]).splitlines():
        if l.startswith("--- "):
            name = l[6:] if l.startswith("--- a/") else None
        elif l.startswith("@@") and name:
            m = re.match(r"@@ -(\d+)(?:,(\d+))?", l)
            got[name].append([int(m.group(1)), int(m.group(2) or 1)])
    return dict(got)


def unpack(commit, into, only=None, without=()):
    tops = [t for t in git("ls-tree", "--name-only", commit).split("\n") if t and t not in SKIP_TOP]
    paths = [t for t in tops if (only is None or t in only) and t not in without]
    into.mkdir(parents=True, exist_ok=True)
    tar = subprocess.Popen(["git", "-C", str(REPO), "archive", commit, "--", *paths], stdout=subprocess.PIPE)
    subprocess.run(["tar", "-x", "-C", str(into)], stdin=tar.stdout, check=True)
    if tar.wait() != 0:
        raise RuntimeError(f"git archive {commit} failed")


def flags(tree):
    eigen = [p + "/eigen3" for p in os.environ.get("NIX_CFLAGS_COMPILE", "").split() if "eigen" in p]
    return (["g++", "-std=c++23", "-O1", "-w", "-DSPATIUM_HAS_BOOST_MULTIPRECISION=1", "-DSPATIUM_HAS_EIGEN=1",
             "-DSPATIUM_HAS_IPC_TOOLKIT=0", f'-DSPATIUM_SOURCE_DIR="{tree}"',
             f'-DRSC_GRAPH_JSON="{tree}/rsc/spatium/graph.json"']
            + [f"-I{tree}/{d}" for d in ("include", "", "examples", "rsc/include")]
            + [x for e in eigen for x in ("-isystem", e)])


def run_one(tree, test):
    """Build one test file of a tree alone and run it: {state, cases: {name: passed}, note}."""
    exe = tree / (pathlib.Path(test).stem + ".run")
    units = [test] + [str(p.relative_to(tree)) for p in sorted((tree / "src/physics").glob("*.cpp"))]
    try:
        r = subprocess.run(flags(tree) + units + ["-o", str(exe), "-lCatch2Main", "-lCatch2"],
                           capture_output=True, text=True, cwd=tree, timeout=BUILD_S)
    except subprocess.TimeoutExpired:
        return {"state": "build timeout", "cases": {}, "note": ""}
    if r.returncode != 0:
        errs = [l for l in r.stderr.splitlines() if "error" in l]
        return {"state": "no build", "cases": {}, "note": (errs[0] if errs else r.stderr[-200:])[-240:]}
    report = exe.with_suffix(".xml")
    try:
        r = subprocess.run([str(exe), "-r", f"junit::out={report}"], capture_output=True, text=True, cwd=tree, timeout=RUN_S)
    except subprocess.TimeoutExpired:
        return {"state": "run timeout", "cases": {}, "note": ""}
    try:
        cases = {c.get("name"): not (c.findall("failure") or c.findall("error"))
                 for c in ET.parse(report).getroot().iter("testcase")}
    except (ET.ParseError, FileNotFoundError):
        return {"state": "crash", "cases": {}, "note": f"exit {r.returncode}: {r.stderr.strip()[-200:]}"}
    return {"state": "ran", "cases": cases, "note": "" if r.returncode in (0,) or not all(cases.values()) else f"exit {r.returncode}"}


def replay(commit, parent, tests, work, keep):
    t0 = time.monotonic()
    home = work / commit[:12]
    shutil.rmtree(home, ignore_errors=True)
    after, before = home / "after", home / "before"
    unpack(commit, after)
    unpack(parent, before, without=("tests",))
    unpack(commit, before, only=("tests",))
    row = {"commit": commit, "parent": parent, "date": git("log", "-1", "--format=%cs", commit).strip(),
           "subject": git("log", "-1", "--format=%s", commit).strip(), "places": places(parent, commit),
           "tests": {}, "breaks": []}
    for t in tests:
        a = run_one(after, t)
        b = run_one(before, t) if a["state"] == "ran" else {"state": "not run", "cases": {}, "note": ""}
        row["tests"][t] = {"after": a, "before": b}
        if a["state"] != "ran":
            continue
        good = {n for n, ok in a["cases"].items() if ok}
        if b["state"] == "ran":
            row["breaks"] += [f"{t}::{n}" for n in sorted(good) if b["cases"].get(n) is False]
        elif b["state"] in ("crash", "run timeout") and len(good) == len(a["cases"]):
            row["breaks"].append(f"{t}::({b['state']} before)")
    row["seconds"] = round(time.monotonic() - t0)
    if not keep:
        shutil.rmtree(home, ignore_errors=True)
    return row


def summary(rows):
    files = [(r, t, v) for r in rows for t, v in r["tests"].items()]
    cell = collections.Counter()
    for r, t, v in files:
        a, b = v["after"], v["before"]
        if a["state"] != "ran":
            cell["harness: after does not " + ("build" if "build" in a["state"] else "run")] += 1
        elif not all(a["cases"].values()):
            cell["harness: a case fails after"] += 1
        elif b["state"] == "no build":
            cell["before does not build (reads a new interface)"] += 1
        elif b["state"] == "ran" and all(b["cases"].get(n, False) for n in a["cases"]):
            cell["before passes every case"] += 1
        elif b["state"] == "ran":
            cell["before fails or lacks a case"] += 1
        else:
            cell["before: " + b["state"]] += 1
    broke = [r for r in rows if r["breaks"]]
    print(f"commits replayed {len(rows)}, test files {len(files)}")
    for k, n in cell.most_common():
        print(f"  {n:4d}  {k}")
    print(f"commits with a break that replays: {len(broke)}; cases that fail before and pass after: {sum(len(r['breaks']) for r in broke)}")
    for r in broke:
        print(f"  {r['commit'][:9]} {r['date']} {r['subject'][:70]}\n      {len(r['breaks'])} cases, e.g. {r['breaks'][0][:110]}"
              f"\n      replaced in {', '.join(sorted(r['places']))[:140]}")


def main():
    arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
    work = pathlib.Path(sys.argv[1]).resolve()
    work.mkdir(parents=True, exist_ok=True)
    out = work / "rows.jsonl"
    rows = [json.loads(l) for l in out.read_text().splitlines()] if out.exists() else []
    if "--summary" in sys.argv:
        return summary(rows)
    done = {r["commit"] for r in rows}
    only = arg("--only")
    todo = [c for c in candidates(arg("--rev", "main")) if c[0] not in done
            and (only is None or any(c[0].startswith(s) for s in only.split(",")))]
    print(f"to replay: {len(todo)} commits ({len(done)} already in {out})", flush=True)
    with cf.ThreadPoolExecutor(int(arg("-j", "4"))) as ex, out.open("a") as f:
        for row in ex.map(lambda c: replay(*c, work, "--keep" in sys.argv), todo):
            f.write(json.dumps(row) + "\n")
            f.flush()
            states = ", ".join(f"{pathlib.Path(t).stem}: {v['before']['state']}/{v['after']['state']}" for t, v in row["tests"].items())
            print(f"{row['seconds']:5d} s  {row['commit'][:9]} breaks {len(row['breaks'])}  {states[:150]}", flush=True)
    summary([json.loads(l) for l in out.read_text().splitlines()])


if __name__ == "__main__":
    main()
