#!/usr/bin/env python3
"""Who uses what: the consumer rule, counted.

The library's rule is that a module with no real caller counts as no code and no tests
(CLAUDE.md, Principles). This script finds, for every name that `gen_capabilities.py` lists as
declared by a header, whether anything outside that header mentions it:

    unused      nothing else in include/, src/, tests/, examples/, benchmarks/ or rsc/ mentions it
    tests-only  only files under tests/ do
    used        something else does

It is approximate and says so: the match is by identifier, not by overload or by call, so a name
shared with another thing reads as used, and a type used only through `auto` or a helper reads as
unused. That is good enough for a ratchet, which is what it is for: `--check` fails when a name
that has no use at all is added (or an old one loses its last use), and prints the names. It does
not fail on tests-only names, which are reported; making those a ratchet too is for after the
existing ones are cleared.

The names that are unused today are in docs/consumers-baseline.json, one per line so that the file
merges. Fixing one (give it a consumer, or delete it) removes it from the baseline with
`--update-baseline`; adding one without a consumer fails the check, and a person has to decide
whether the baseline is the place for it.

    python3 scripts/gen_consumers.py                  # the table per domain
    python3 scripts/gen_consumers.py --list           # and the names
    python3 scripts/gen_consumers.py --check          # CI
    python3 scripts/gen_consumers.py --update-baseline
"""
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gen_capabilities import INCLUDE_ROOT, ROOT, description, names  # noqa: E402

BASELINE = ROOT / "docs" / "consumers-baseline.json"
SEARCH = {"include": "include/spatium", "tests": "tests", "examples": "examples", "src": "src",
          "rsc": "rsc", "bench": "benchmarks"}
IDENT = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")


def identifiers(path):
    try:
        return set(IDENT.findall(path.read_text(encoding="utf-8", errors="ignore")))
    except OSError:
        return set()


def survey():
    """{(domain, header, name): 'unused' | 'tests-only' | 'used'}"""
    files = {}
    for cat, d in SEARCH.items():
        for dp, _, fs in os.walk(ROOT / d):
            for f in fs:
                if f.endswith((".hpp", ".cpp", ".cppm", ".h")):
                    p = Path(dp) / f
                    files[p] = (cat, identifiers(p))
    result = {}
    for path in sorted(INCLUDE_ROOT.rglob("*.hpp")):
        rel = path.relative_to(INCLUDE_ROOT).as_posix()
        if rel == "_export_macro.hpp":
            continue
        lines = path.read_text(encoding="utf-8").splitlines()
        if description(lines).lower().startswith("internal"):
            continue                      # a header that says it is not for use is not asked to have users
        domain = rel.split("/")[0] if "/" in rel else "(top level)"
        for name in names(lines):
            cats = {cat for p, (cat, ids) in files.items() if p != path and name in ids}
            result[(domain, rel, name)] = "used" if cats - {"tests"} else ("tests-only" if cats else "unused")
    return result


def summarize(result):
    by = defaultdict(lambda: defaultdict(int))
    for (domain, _, _), kind in result.items():
        by[domain][kind] += 1
    return by


def baseline_names(result):
    return sorted(f"{rel}: {name}" for (_, rel, name), kind in result.items() if kind == "unused")


def main():
    result = survey()
    unused_now = baseline_names(result)
    if "--update-baseline" in sys.argv:
        BASELINE.write_text(json.dumps({"unused": unused_now}, indent=0) + "\n", encoding="utf-8")
        print(f"wrote {BASELINE.relative_to(ROOT)}: {len(unused_now)} names with no use")
        return 0
    if "--check" in sys.argv:
        old = set(json.loads(BASELINE.read_text(encoding="utf-8"))["unused"]) if BASELINE.exists() else set()
        added = [n for n in unused_now if n not in old]
        gone = [n for n in old if n not in set(unused_now)]
        tests_only = sum(1 for k in result.values() if k == "tests-only")
        if added:
            print("names with no consumer anywhere (the library's rule: a module nothing calls counts as no code):", file=sys.stderr)
            for n in added:
                print(f"  {n}", file=sys.stderr)
            print("give each a caller (a function, a scene, a test through its consumer) or delete it; if the baseline "
                  "really is the place for one, run scripts/gen_consumers.py --update-baseline and say why in review.",
                  file=sys.stderr)
            return 1
        print(f"consumers: {len(unused_now)} names with no use (baseline {len(old)}, {len(gone)} fixed since: "
              f"run --update-baseline to keep the gain), {tests_only} used only by tests")
        return 0
    by = summarize(result)
    print(f"{'domain':<12}{'names':>7}{'unused':>8}{'tests-only':>12}")
    for d, c in sorted(by.items(), key=lambda kv: -kv[1]["unused"]):
        print(f"{d:<12}{sum(c.values()):>7}{c['unused']:>8}{c['tests-only']:>12}")
    total = len(result)
    print(f"{'all':<12}{total:>7}{sum(1 for k in result.values() if k == 'unused'):>8}"
          f"{sum(1 for k in result.values() if k == 'tests-only'):>12}")
    if "--list" in sys.argv:
        for kind in ("unused", "tests-only"):
            print(f"\n{kind}:")
            for (domain, rel, name), k in sorted(result.items()):
                if k == kind:
                    print(f"  {rel}: {name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
