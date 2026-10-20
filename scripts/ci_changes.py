#!/usr/bin/env python3
"""Which CI jobs a change needs.

A pull request that only edits documentation does not need nine builds and a matrix of
1168 cells, and a change to geometry/ does not need the matrix at all: its cells include
`probes.hpp`, whose include closure is core/, algebra/ and spaces/ (checked with
`g++ -MM`, 2026-10-08). This script reads the changed paths and says, for the two groups
of heavy jobs, whether they must run:

    code    every build and test job (anything that is not documentation)
    matrix  the connectivity matrix (only what its cells read)
    rsc     the RSC checks: rsc/core stays free of the library, rsc/spatium/graph.json is what its
            generator makes of the headers and the matrix, and the report and the witness pairs
            are what rsc_find.py makes of the graph

Used by the `changes` job of .github/workflows/ci.yml, which prints `code=`, `matrix=` and `rsc=`
lines for $GITHUB_OUTPUT. Without a usable base (a first push, a scheduled or manual run)
everything runs. `--selftest` holds the rules to examples; CI runs it.

    git diff --name-only BASE HEAD | python3 scripts/ci_changes.py
    python3 scripts/ci_changes.py --all
    python3 scripts/ci_changes.py --selftest
"""

import re
import sys

# Not built, not tested: the heavy jobs have nothing to say about these.
DOCS_ONLY = [
    r"^docs/(?!connectivity\.md$)",       # docs/connectivity.md is checked by the matrix itself
    r"\.md$",
    r"^LICENSE$",
    r"^gallery/",
    r"\.(png|svg|mp4|gif|jpg)$",
    r"^scripts/(gen_roadmap|gen_capabilities|gen_consumers|gen_dependency_graph|check_claude_md_layout|check_doc_file_refs)\.py$",
]

# What the matrix cells read: probes.hpp includes only algebra/, core/ and spaces/ (and the
# export macro), the generator compiles them itself, and the compiler comes from the flake.
MATRIX_INPUTS = [
    r"^include/spatium/(core|algebra|spaces)/",
    r"^include/spatium/_export_macro\.hpp$",
    r"^tests/connectivity/",
    r"^scripts/gen_connectivity\.py$",
    r"^docs/connectivity\.md$",
    r"^CMakeLists\.txt$",
    r"^cmake/",
    r"^flake\.(nix|lock)$",
    r"^\.github/workflows/",
    r"^scripts/ci_changes\.py$",
]

# What rsc/ reads or is: the graph of the library is generated from every concept under include/
# and from the matrix's expected levels (and the generator imports gen_connectivity.py), so any
# header change can change it, and so can the sampler and the pair machinery the finder runs routes on
# (rsc_find.py --check). A Markdown file under rsc/ changes nothing.
RSC_INPUTS = [
    r"^rsc/(?!.*\.md$)",
    r"^include/spatium/",
    r"^tests/connectivity/",
    r"^tests/symmetry/(pair|rsc_pair_support)\.hpp$",   # what the finder's runs and the witness pairs are built on
    r"^tests/symmetry/pairs/rsc_",                       # the witness pairs it writes
    r"^scripts/gen_connectivity\.py$",
    r"^CMakeLists\.txt$",
    r"^cmake/",
    r"^flake\.(nix|lock)$",
    r"^\.github/workflows/",
    r"^scripts/ci_changes\.py$",
]


def _any(patterns, path):
    return any(re.search(p, path) for p in patterns)


def classify(paths):
    paths = [p.strip() for p in paths if p.strip()]
    if not paths:
        return True, True, True               # an empty diff proves nothing: run everything
    code = any(not _any(DOCS_ONLY, p) for p in paths)
    matrix = any(_any(MATRIX_INPUTS, p) for p in paths)
    rsc = any(_any(RSC_INPUTS, p) for p in paths)
    return code, matrix, rsc


def selftest():
    # (changed paths, (code, matrix, rsc))
    cases = [
        (["docs/ROADMAP.md"], (False, False, False)),
        (["docs/roadmap/2026-10-08-01-x.md", "CLAUDE.md"], (False, False, False)),
        (["gallery/boom.png", "README.md"], (False, False, False)),
        (["include/spatium/geometry/polygon.hpp"], (True, False, True)),
        (["include/spatium/physics/mechanics/narrow_phase.hpp", "tests/symmetry/pairs.hpp"], (True, False, True)),
        (["include/spatium/algebra/polynomial.hpp"], (True, True, True)),
        (["include/spatium/spaces/spd.hpp", "docs/roadmap/x.md"], (True, True, True)),
        (["tests/connectivity/probes.hpp"], (True, True, True)),
        (["docs/connectivity.md"], (False, True, False)),
        ([".github/workflows/ci.yml"], (True, True, True)),
        (["CMakeLists.txt"], (True, True, True)),
        (["examples/donut_demo.cpp"], (True, False, False)),
        (["tests/test_rsc_graph.cpp"], (True, False, False)),
        (["rsc/README.md"], (False, False, False)),
        (["scripts/gen_roadmap.py"], (False, False, False)),
        (["scripts/gen_connectivity.py"], (True, True, True)),
        (["rsc/core/graph.hpp"], (True, False, True)),
        (["rsc/spatium/graph.json"], (True, False, True)),
        (["rsc/spatium/gen_graph.py", "rsc/README.md"], (True, False, True)),
        (["rsc/spatium/rsc_find.py"], (True, False, True)),
        (["tests/symmetry/rsc_pair_support.hpp"], (True, False, True)),
        (["tests/symmetry/pairs/rsc_Sphere_2_distance.cpp"], (True, False, True)),
        (["tests/symmetry/pairs/sphere_exp_log.cpp"], (True, False, False)),
        ([], (True, True, True)),
    ]
    bad = 0
    for paths, want in cases:
        got = classify(paths)
        if got != want:
            bad += 1
            print(f"FAIL {paths}: got {got}, want {want}")
    if bad:
        sys.exit(1)
    print(f"ci_changes: {len(cases)} cases hold")


def main():
    if "--selftest" in sys.argv:
        return selftest()
    if "--all" in sys.argv:
        code, matrix, rsc = True, True, True
    else:
        code, matrix, rsc = classify(sys.stdin.read().splitlines())
    print(f"code={'true' if code else 'false'}")
    print(f"matrix={'true' if matrix else 'false'}")
    print(f"rsc={'true' if rsc else 'false'}")


if __name__ == "__main__":
    main()
