#!/usr/bin/env python3
"""The history entries of the roadmap: docs/roadmap/*.md, one file per entry.

An entry is `docs/roadmap/<date>-<nn>-<slug>.md`:

    ---
    date: 2026-10-07
    title: Riemannian descent on a space with intrinsic coordinates
    stage: 2.1                      # optional: the plan stage it belongs to
    closes: [16, 25]                # optional: findings rows this entry closed (deleted)
    files: [include/spatium/algebra/calculus.hpp, tests/test_calculus.cpp]   # optional
    ---
    The paragraph(s), as they should read.

The files are the roadmap's history, listed by GitHub in name order (the date first).
They are NOT rendered into docs/ROADMAP.md, and that is deliberate: a rendering
committed with every change is a shared generated file that two changes in flight both
edit at its end, and conflict on -- which is the problem the files exist to end. Instead

    scripts/gen_roadmap.py            prints the entries as one page (`> page.md` to keep it)
    scripts/gen_roadmap.py --check    validates them (run in CI)

`--check` fails when
  * an entry lacks a title or a date, or its date is not the one in its file name;
  * a listed file does not exist in the tree;
  * a findings row an entry says it closed is still in docs/findings-2026-09.md.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ENTRIES = ROOT / "docs" / "roadmap"
FINDINGS = ROOT / "docs" / "findings-2026-09.md"


def parse(path):
    text = path.read_text(encoding="utf-8")
    m = re.match(r"---\n(.*?)\n---\n(.*)\Z", text, re.S)
    if not m:
        raise ValueError(f"{path.name}: no front matter between --- lines")
    meta = {}
    for line in m.group(1).split("\n"):
        line = line.split("  #")[0].rstrip()
        if not line.strip():
            continue
        k, _, v = line.partition(":")
        v = v.strip()
        if v.startswith("["):
            v = [x.strip() for x in v.strip("[]").split(",") if x.strip()]
        meta[k.strip()] = v
    return meta, m.group(2).strip()


def problems(path, meta):
    out = []
    name = path.name
    if not meta.get("title"):
        out.append(f"{name}: no title")
    d = meta.get("date", "")
    if not re.fullmatch(r"\d{4}-\d\d-\d\d", d):
        out.append(f"{name}: date '{d}' is not YYYY-MM-DD")
    elif not name.startswith(d + "-"):
        out.append(f"{name}: the file name must start with its date {d}")
    for f in meta.get("files", []):
        if not (ROOT / f).exists():
            out.append(f"{name}: lists {f}, which is not in the tree")
    findings = FINDINGS.read_text(encoding="utf-8") if FINDINGS.exists() else ""
    for n in meta.get("closes", []):
        if re.search(rf"^\|\s*{re.escape(str(n))}\s*\|", findings, re.M):
            out.append(f"{name}: says it closed findings row {n}, which is still in {FINDINGS.name}")
    return out


def render():
    paragraphs, errors = [], []
    for path in sorted(ENTRIES.glob("*.md")):
        try:
            meta, body = parse(path)
        except ValueError as e:
            errors.append(str(e))
            continue
        errors += problems(path, meta)
        paragraphs.append(f"**{meta.get('title', '?')} ({meta.get('date', '?')}).** {body}")
    return "\n\n".join(paragraphs), errors


def main():
    block, errors = render()
    if "--check" in sys.argv:
        for b in errors:
            print(b, file=sys.stderr)
        if errors:
            return 1
        print(f"docs/roadmap: {len(list(ENTRIES.glob('*.md')))} entries are valid")
        return 0
    for e in errors:
        print("warning:", e, file=sys.stderr)
    print(block)
    return 0


if __name__ == "__main__":
    sys.exit(main())
