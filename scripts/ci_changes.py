#!/usr/bin/env python3
"""Which CI jobs a change needs.

A pull request that only edits documentation does not need nine builds and a matrix of
1168 cells, and a change to geometry/ does not need the matrix at all: its cells include
`probes.hpp`, whose include closure is core/, algebra/ and spaces/ (checked with
`g++ -MM`, 2026-10-08). This script reads the changed paths and says, for the two groups
of heavy jobs, whether they must run:

    code    every build and test job (anything that is not documentation)
    matrix  the connectivity matrix (only what its cells read)

Used by the `changes` job of .github/workflows/ci.yml, which prints `code=` and `matrix=`
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


def _any(patterns, path):
    return any(re.search(p, path) for p in patterns)


def classify(paths):
    paths = [p.strip() for p in paths if p.strip()]
    if not paths:
        return True, True                     # an empty diff proves nothing: run everything
    code = any(not _any(DOCS_ONLY, p) for p in paths)
    matrix = any(_any(MATRIX_INPUTS, p) for p in paths)
    return code, matrix


def selftest():
    cases = [
        (["docs/ROADMAP.md"], (False, False)),
        (["docs/roadmap/2026-10-08-01-x.md", "CLAUDE.md"], (False, False)),
        (["gallery/boom.png", "README.md"], (False, False)),
        (["include/spatium/geometry/polygon.hpp"], (True, False)),
        (["include/spatium/physics/mechanics/narrow_phase.hpp", "tests/symmetry/pairs.hpp"], (True, False)),
        (["include/spatium/algebra/polynomial.hpp"], (True, True)),
        (["include/spatium/spaces/spd.hpp", "docs/roadmap/x.md"], (True, True)),
        (["tests/connectivity/probes.hpp"], (True, True)),
        (["docs/connectivity.md"], (False, True)),
        ([".github/workflows/ci.yml"], (True, True)),
        (["CMakeLists.txt"], (True, True)),
        (["rsc/README.md"], (False, False)),
        (["scripts/gen_roadmap.py"], (False, False)),
        (["scripts/gen_connectivity.py"], (True, True)),
        (["rsc/include/dispatcher.hpp"], (True, False)),
        ([], (True, True)),
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
        code, matrix = True, True
    else:
        code, matrix = classify(sys.stdin.read().splitlines())
    print(f"code={'true' if code else 'false'}")
    print(f"matrix={'true' if matrix else 'false'}")


if __name__ == "__main__":
    main()
