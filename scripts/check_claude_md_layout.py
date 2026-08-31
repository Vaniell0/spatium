#!/usr/bin/env python3
"""Check that CLAUDE.md's Directory layout section mentions every real
include/spatium/ subdirectory.

This exists because it's exactly the class of staleness found by hand
2026-08-28: CLAUDE.md's `physics/` entry still described the pre-split
layout after physics/atomic/ was carved out, and the `algebra/` entry
was missing several files added over the prior months. A directory that
holds headers but isn't named anywhere in CLAUDE.md is a silent version
of the same problem -- this catches that mechanically instead of relying
on the next person (or agent) to notice by hand.

Not a full content check (doesn't verify each directory's *description*
is accurate, only that the directory is mentioned at all) -- narrower
and cheaper than that, deliberately.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
INCLUDE_ROOT = ROOT / "include" / "spatium"
CLAUDE_MD = ROOT / "CLAUDE.md"


def real_content_dirs() -> set[str]:
    dirs = set()
    for hpp in INCLUDE_ROOT.rglob("*.hpp"):
        rel_dir = hpp.parent.relative_to(INCLUDE_ROOT).as_posix()
        if rel_dir != ".":
            dirs.add(rel_dir)
    return dirs


def directory_layout_section(text: str) -> str:
    match = re.search(r"## Directory layout\n(.*?)\n## ", text, re.S)
    return match.group(1) if match else ""


def main() -> int:
    dirs = real_content_dirs()
    section = directory_layout_section(CLAUDE_MD.read_text(encoding="utf-8"))
    if not section:
        sys.stderr.write("CLAUDE.md: couldn't find a '## Directory layout' section.\n")
        return 1

    missing = sorted(d for d in dirs if f"include/spatium/{d}/" not in section)
    # vendor/ is third-party (stb_image_write), not project structure to document.
    missing = [d for d in missing if not d.startswith("vendor")]

    if missing:
        sys.stderr.write(
            "CLAUDE.md's Directory layout doesn't mention:\n"
            + "\n".join(f"  include/spatium/{d}/" for d in missing)
            + "\nAdd a line for each, or this directory's own header comment "
            "if it explains itself.\n"
        )
        return 1

    print(f"CLAUDE.md's Directory layout covers all {len(dirs)} real content directories.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
