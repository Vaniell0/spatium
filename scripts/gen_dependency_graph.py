#!/usr/bin/env python3
"""Generate docs/dependency-graph.dot from real #include <spatium/...> edges.

Aggregates to domain level (the first path component under
include/spatium/, e.g. `algebra`, `geometry`, `physics`) rather than
file level -- a file-level graph over 100+ headers is too noisy to read
or to eyeball for a stale layering claim. Files directly under
include/spatium/ (Point, Morphism, the umbrella headers) are their own
pseudo-domain, "spatium".

This exists so docs/architecture.md's layering description can be
checked against the real graph instead of drifting silently -- see
docs/conventions.md and the CI step that runs this script and diffs the
output against what's committed.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
INCLUDE_ROOT = ROOT / "include" / "spatium"
OUTPUT = ROOT / "docs" / "dependency-graph.dot"

INCLUDE_RE = re.compile(r'#\s*include\s*<spatium/([^>]+)>')


def domain_of(relative_path: str) -> str:
    # Files directly under include/spatium/ (point.hpp, morphism.hpp, and
    # the escalating convenience umbrellas core.hpp/mesh.hpp/physics.hpp,
    # each of which deliberately spans several domains at once -- see
    # their own header comments) are the "spatium" pseudo-domain. A file
    # named e.g. mesh.hpp is NOT the mesh/ domain's own umbrella -- it's
    # a broader bundle that happens to be named after the biggest new
    # thing it adds, so no stem-matching special case here.
    parts = relative_path.split("/")
    return parts[0] if len(parts) > 1 else "spatium"


def main() -> int:
    edges = set()
    for path in sorted(INCLUDE_ROOT.rglob("*.hpp")):
        rel = path.relative_to(INCLUDE_ROOT).as_posix()
        src_domain = domain_of(rel)
        text = path.read_text(encoding="utf-8", errors="ignore")
        for match in INCLUDE_RE.finditer(text):
            dst_domain = domain_of(match.group(1))
            if dst_domain != src_domain:
                edges.add((src_domain, dst_domain))

    lines = ["digraph spatium_dependencies {", "    rankdir=LR;"]
    for src, dst in sorted(edges):
        lines.append(f'    "{src}" -> "{dst}";')
    lines.append("}")
    output = "\n".join(lines) + "\n"

    if "--check" in sys.argv:
        current = OUTPUT.read_text(encoding="utf-8") if OUTPUT.exists() else ""
        if current != output:
            sys.stderr.write(
                f"{OUTPUT} is stale -- run `python3 scripts/gen_dependency_graph.py` "
                "and commit the result.\n"
            )
            return 1
        return 0

    OUTPUT.write_text(output, encoding="utf-8")
    print(f"Wrote {OUTPUT} ({len(edges)} domain-level edges).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
