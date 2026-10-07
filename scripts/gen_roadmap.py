#!/usr/bin/env python3
"""Generate the history entries of docs/ROADMAP.md from docs/roadmap/*.md.

An entry of the roadmap's history -- what was done, what was measured, what was
not done -- is one file, `docs/roadmap/<date>-<nn>-<slug>.md`:

    ---
    date: 2026-10-07
    title: Riemannian descent on a space with intrinsic coordinates
    stage: 2.1                      # optional: the plan stage it belongs to
    closes: [16, 25]                # optional: findings rows this entry closed (deleted)
    files: [include/spatium/algebra/calculus.hpp, tests/test_calculus.cpp]   # optional
    ---
    The paragraph(s), as they should read in ROADMAP.md.

and this script renders them, in file-name order, between the markers

    <!-- roadmap-entries:begin -->  ...  <!-- roadmap-entries:end -->

of docs/ROADMAP.md as `**Title (date).** body`. Why files: every change used to
append a paragraph at the same place of one 1700-line file, so any two changes in
flight conflicted at every merge; now each change adds its own file, and the file
is checked rather than trusted.

`--check` (run in CI) fails when
  * the generated block is stale;
  * an entry lacks a title or a date, or its date is not the one in its file name;
  * a listed file does not exist in the tree;
  * a findings row an entry says it closed is still in docs/findings-2026-09.md.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ENTRIES = ROOT / "docs" / "roadmap"
ROADMAP = ROOT / "docs" / "ROADMAP.md"
FINDINGS = ROOT / "docs" / "findings-2026-09.md"
BEGIN = "<!-- roadmap-entries:begin -->"
END = "<!-- roadmap-entries:end -->"


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
    text = ROADMAP.read_text(encoding="utf-8")
    if BEGIN not in text or END not in text:
        print(f"docs/ROADMAP.md has no {BEGIN} ... {END} markers", file=sys.stderr)
        return 1
    head, rest = text.split(BEGIN, 1)
    _, tail = rest.split(END, 1)
    new = head + BEGIN + "\n\n" + block + "\n\n" + END + tail
    if "--check" in sys.argv:
        bad = list(errors)
        if new != text:
            bad.append("docs/ROADMAP.md is stale; run scripts/gen_roadmap.py")
        for b in bad:
            print(b, file=sys.stderr)
        if bad:
            return 1
        print(f"docs/ROADMAP.md is up to date ({len(list(ENTRIES.glob('*.md')))} entries)")
        return 0
    for e in errors:
        print("warning:", e, file=sys.stderr)
    ROADMAP.write_text(new, encoding="utf-8")
    print(f"wrote {ROADMAP.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
