#!/usr/bin/env python3
"""The roadmap's history, by area: docs/roadmap/<area>/<slug>.md, one file per entry.

An area is a folder with a README.md: where the area stands, written by hand, and its
reasoning and earlier history. The areas follow the layers and domains of the library (ir,
document, evaluation, lowering, views, rsc, algebra, spaces, geometry, physics, relativity,
discrete, sound, tooling), so the roadmap grows where the work is, not in a list of dates.
An entry is a file named for what it did, not when:

    ---
    date: 2026-10-07
    title: Riemannian descent on a space with intrinsic coordinates
    stage: 2.1                      # optional: the plan stage it belongs to
    closes: [16, 25]                # optional: findings rows this entry closed (deleted)
    files: [include/spatium/algebra/calculus.hpp, tests/test_calculus.cpp]   # optional
    open:                           # optional: what the entry did NOT do, one item per line
      - restrict, quotient and compactify
      - a Sobol sequence
    resolves: [algebra/sobol-where-halton-fails#2]   # optional: open items of earlier entries this one did
    ---
    The paragraph(s), as they should read.

Nothing is rendered into a committed file, and that is deliberate: a page regenerated with
every change is a shared file that two changes in flight both edit, and conflict on -- which
is what one file per entry exists to end. A new entry is a new file in its area. Instead

    scripts/gen_roadmap.py                    prints every area's entries, by date
    scripts/gen_roadmap.py --area evaluation  one area's entries and what it has not done
    scripts/gen_roadmap.py --check            validates them (run in CI)
    scripts/gen_roadmap.py --open             what is still not done, by area: every `open`
                                              item no later entry `resolves`

An open item's id is `<area>/<slug>#<k>` (the entry's path in docs/roadmap without .md, and
the item's number from 1): `resolves` names those, and `--check` fails on one that names
nothing. "What was not done" is then a list that is computed, not a sentence at the end of a
paragraph, where it was lost each time.

`--check` fails when
  * an entry is not in an area folder, or an area folder has no README.md;
  * an entry lacks a title, or its date is not YYYY-MM-DD;
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


def entry_id(path):
    """<area>/<slug>: the entry's path in docs/roadmap without .md."""
    return f"{path.parent.name}/{path.stem}"


def item_ids(path, meta):
    """The ids of an entry's open items: <area>/<slug>#<k>."""
    return [f"{entry_id(path)}#{k}" for k in range(1, len(meta.get("open", [])) + 1)]


def problems(path, meta, all_ids=frozenset()):
    out = []
    name = entry_id(path)
    if not meta.get("title"):
        out.append(f"{name}: no title")
    d = meta.get("date", "")
    if not re.fullmatch(r"\d{4}-\d\d-\d\d", d):
        out.append(f"{name}: date '{d}' is not YYYY-MM-DD")
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


def areas():
    return sorted(p.name for p in ENTRIES.iterdir() if p.is_dir())


def load(area=None):
    entries, errors = [], []
    for stray in sorted(ENTRIES.glob("*.md")):
        errors.append(f"{stray.name}: an entry belongs in an area folder, docs/roadmap/<area>/")
    for a in areas():
        if not (ENTRIES / a / "README.md").exists():
            errors.append(f"{a}/: an area has a README.md, where it stands")
    paths = [p for p in ENTRIES.glob("*/*.md") if p.name != "README.md"]
    for path in sorted(paths, key=lambda p: (p.parent.name, p.stem)):
        try:
            meta, body = parse(path)
        except ValueError as e:
            errors.append(str(e))
            continue
        entries.append((path, meta, body))
    ids = {i for path, meta, _ in entries for i in item_ids(path, meta)}
    for path, meta, _ in entries:
        errors += problems(path, meta, ids)
    if area:
        entries = [e for e in entries if e[0].parent.name == area]
    entries.sort(key=lambda e: (e[0].parent.name, e[1].get("date", ""), e[0].stem))
    return entries, errors


def render(area=None):
    entries, errors = load(area)
    out, current = [], None
    for path, m, body in entries:
        if path.parent.name != current:
            current = path.parent.name
            out.append(f"## {current}")
        out.append(f"**{m.get('title', '?')} ({m.get('date', '?')}).** {body}")
    return "\n\n".join(out), errors


def open_items(area=None):
    """[(id, area, item, entry title)] not resolved by any entry."""
    entries, errors = load()
    resolved = {r for _, m, _ in entries for r in m.get("resolves", [])}
    out = []
    for path, meta, _ in entries:
        if area and path.parent.name != area:
            continue
        for i, item in zip(item_ids(path, meta), meta.get("open", [])):
            if i not in resolved:
                out.append((i, path.parent.name, item, meta.get("title", "?")))
    return out, errors


def main():
    area = sys.argv[sys.argv.index("--area") + 1] if "--area" in sys.argv else None
    if area and area not in areas():
        print(f"no area {area}; the areas are {', '.join(areas())}", file=sys.stderr)
        return 1
    if "--open" in sys.argv or area:
        items, errors = open_items(area)
        for e in errors:
            print("warning:", e, file=sys.stderr)
        if area:
            block, _ = render(area)
            print(block + "\n")
        by_area = {}
        for i, a, item, title in items:
            by_area.setdefault(a, []).append((i, item, title))
        print(f"{len(items)} open items, not resolved by a later entry\n")
        for a in sorted(by_area):
            print(a)
            for i, item, title in by_area[a]:
                print(f"  {i}  {item}   [{title}]")
            print()
        return 0
    block, errors = render()
    if "--check" in sys.argv:
        for b in errors:
            print(b, file=sys.stderr)
        if errors:
            return 1
        items, _ = open_items()
        entries, _ = load()
        print(f"docs/roadmap: {len(entries)} entries in {len(areas())} areas are valid, {len(items)} open items")
        return 0
    for e in errors:
        print("warning:", e, file=sys.stderr)
    print(block)
    return 0


if __name__ == "__main__":
    sys.exit(main())
