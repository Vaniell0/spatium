#!/usr/bin/env python3
"""The history entries of the roadmap: docs/roadmap/*.md, one file per entry.

An entry is `docs/roadmap/<date>-<nn>-<slug>.md`:

    ---
    date: 2026-10-07
    title: Riemannian descent on a space with intrinsic coordinates
    stage: 2.1                      # optional: the plan stage it belongs to
    closes: [16, 25]                # optional: findings rows this entry closed (deleted)
    files: [include/spatium/algebra/calculus.hpp, tests/test_calculus.cpp]   # optional
    open:                           # optional: what the entry did NOT do, one item per line
      - restrict, quotient and compactify
      - a Sobol sequence
    resolves: [2026-10-07-12#2]     # optional: open items of earlier entries this one did
    ---
    The paragraph(s), as they should read.

The files are the roadmap's history, listed by GitHub in name order (the date first).
They are NOT rendered into docs/ROADMAP.md, and that is deliberate: a rendering
committed with every change is a shared generated file that two changes in flight both
edit at its end, and conflict on -- which is the problem the files exist to end. Instead

    scripts/gen_roadmap.py            prints the entries as one page (`> page.md` to keep it)
    scripts/gen_roadmap.py --check    validates them (run in CI)
    scripts/gen_roadmap.py --open     lists what is still not done: every `open` item no later
                                      entry `resolves`, grouped by plan stage

An open item's id is `<date>-<nn>#<k>` (the entry's file-name prefix and the item's number
from 1): `resolves` names those, and `--check` fails on one that names nothing. The point is
that "what was not done" stops being a sentence at the end of a paragraph, where it was lost
each time, and becomes a list that is computed.

`--check` fails when
  * an entry lacks a title or a date, or its date is not the one in its file name;
  * a listed file does not exist in the tree;
  * a findings row an entry says it closed is still in docs/findings-2026-09.md;
  * a `resolves` id names no entry, or an item number the entry does not have.
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
    meta, block = {}, None
    for line in m.group(1).split("\n"):
        if block is not None and re.match(r"^\s+-\s", line):
            meta[block].append(line.strip()[1:].strip())
            continue
        block = None
        line = line.split("  #")[0].rstrip()
        if not line.strip():
            continue
        k, _, v = line.partition(":")
        v = v.strip()
        if not v:
            block = k.strip()
            meta[block] = []
            continue
        if v.startswith("["):
            v = [x.strip() for x in v.strip("[]").split(",") if x.strip()]
        meta[k.strip()] = v
    return meta, m.group(2).strip()


def item_ids(path, meta):
    """The ids of an entry's open items: <date>-<nn>#<k>."""
    return [f"{path.stem[:13]}#{k}" for k in range(1, len(meta.get("open", [])) + 1)]


def problems(path, meta, all_ids=frozenset()):
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
    for rid in meta.get("resolves", []):
        if rid not in all_ids:
            out.append(f"{name}: resolves {rid}, which is not an open item of any entry")
    findings = FINDINGS.read_text(encoding="utf-8") if FINDINGS.exists() else ""
    for n in meta.get("closes", []):
        if re.search(rf"^\|\s*{re.escape(str(n))}\s*\|", findings, re.M):
            out.append(f"{name}: says it closed findings row {n}, which is still in {FINDINGS.name}")
    return out


def load():
    entries, errors = [], []
    for path in sorted(ENTRIES.glob("*.md")):
        try:
            meta, body = parse(path)
        except ValueError as e:
            errors.append(str(e))
            continue
        entries.append((path, meta, body))
    ids = {i for path, meta, _ in entries for i in item_ids(path, meta)}
    for path, meta, _ in entries:
        errors += problems(path, meta, ids)
    return entries, errors


def render():
    entries, errors = load()
    paragraphs = [f"**{m.get('title', '?')} ({m.get('date', '?')}).** {body}" for _, m, body in entries]
    return "\n\n".join(paragraphs), errors


def open_items():
    """[(id, stage, item, entry title)] not resolved by any entry."""
    entries, errors = load()
    resolved = {r for _, m, _ in entries for r in m.get("resolves", [])}
    out = []
    for path, meta, _ in entries:
        for i, item in zip(item_ids(path, meta), meta.get("open", [])):
            if i not in resolved:
                out.append((i, meta.get("stage", "-"), item, meta.get("title", "?")))
    return out, errors


def main():
    if "--open" in sys.argv:
        items, errors = open_items()
        for e in errors:
            print("warning:", e, file=sys.stderr)
        by_stage = {}
        for i, stage, item, title in items:
            by_stage.setdefault(stage, []).append((i, item, title))
        print(f"{len(items)} open items, not resolved by a later entry\n")
        for stage in sorted(by_stage, key=lambda x: [int(p) if p.isdigit() else 0 for p in re.split(r"[.]", x)]):
            print(f"stage {stage}")
            for i, item, title in by_stage[stage]:
                print(f"  {i:20s} {item}   [{title}]")
            print()
        return 0
    block, errors = render()
    if "--check" in sys.argv:
        for b in errors:
            print(b, file=sys.stderr)
        if errors:
            return 1
        items, _ = open_items()
        print(f"docs/roadmap: {len(list(ENTRIES.glob('*.md')))} entries are valid, {len(items)} open items")
        return 0
    for e in errors:
        print("warning:", e, file=sys.stderr)
    print(block)
    return 0


if __name__ == "__main__":
    sys.exit(main())
