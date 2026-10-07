#!/usr/bin/env python3
"""Generate docs/capabilities.md: what the library can do, found from the code.

For every header under include/spatium/ it records the first sentence of its
own descriptive comment and the names it declares at namespace level (types,
concepts, free functions), grouped by domain. Nothing here is written by
hand, so it cannot fall behind the way a hand-kept API reference does: a new
header or a new public name shows up in the diff of this file, and CI runs
`--check` the way it does for the connectivity matrix and the dependency
graph.

It is a finder, not a manual. The sentence says what the header is for; the
names say where to grep next. A header with no descriptive comment is listed
at the top as such, and that count is the number to keep going down.

Heuristics, stated because they are heuristics: a description is the first
run of two or more `//` lines in the first 150 lines; names are declarations
at column 0 (the code style puts namespace-level items there), outside
`*_detail` namespaces; member functions are not listed -- their class is.
"""
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
INCLUDE_ROOT = ROOT / "include" / "spatium"
OUTPUT = ROOT / "docs" / "capabilities.md"

TYPE_RE = re.compile(r"^(?:struct|class|enum class|enum)\s+(\w+)")
CONCEPT_RE = re.compile(r"^concept\s+(\w+)\s*=")
USING_RE = re.compile(r"^using\s+(\w+)\s*=")
FUNC_RE = re.compile(r"^(?:\[\[\w+\]\]\s*)?(?:(?:inline|constexpr|static|friend|explicit|consteval)\s+)*"
                     r"[\w:<>,\*&\s~]+?[\s\*&>](\w+)\s*\(")
SKIP_FIRST_WORDS = {"if", "for", "while", "return", "switch", "static_assert", "using", "namespace",
                    "template", "typedef", "else", "case", "do", "throw", "requires", "operator", "friend"}
MAX_NAMES = 14


SENTENCE_END = re.compile(r"(?<!\be\.g)(?<!\bi\.e)(?<!\bvs)(?<!\betc)(?<!\bcf)[.?!](?=\s|$)")


def clean(text):
    text = re.sub(r"[─━═]+", " ", text)
    text = re.sub(r"^\s*-*\s*spatium/\S+\s*-*\s*", "", text)
    return re.sub(r"\s+", " ", text).strip()


def description(lines):
    """First sentence of the header-level comment: the first run of >= 2
    comment lines that comes before any declaration. A comment that belongs
    to a declaration is not the header's description."""
    run = []
    for line in lines[:200]:
        s = line.strip()
        if s.startswith("//"):
            run.append(s[2:].strip())
            continue
        if len(run) >= 2:
            break
        run = []
        if not s or s.startswith("#") or s.startswith("}") or re.match(
                r"^(SPATIUM_EXPORT\s+)?(inline\s+)?namespace\b", s) or s.startswith("SPATIUM_EXPORT"):
            continue
        return ""           # a declaration came first: the header has no description of its own
    text = clean(" ".join(t for t in run if t))
    if not text:
        return ""
    m = SENTENCE_END.search(text)
    sentence = text[: m.start() + 1] if m else text
    return sentence[:200] + ("..." if len(sentence) > 200 else "")


def names(lines):
    found = []
    in_detail = False
    for line in lines:
        if not line or line[0] in " \t#/}":
            if in_detail and line.startswith("}") and "detail" in line:
                in_detail = False
            continue
        if re.match(r"^namespace\s+\w*detail\w*\s*\{", line):
            in_detail = True
            continue
        if in_detail:
            continue
        stripped = line
        stripped = re.sub(r"^template\s*<.*?>\s*", "", stripped)
        stripped = re.sub(r"^requires\s.*$", "", stripped)
        if not stripped:
            continue
        for rx in (TYPE_RE, CONCEPT_RE, USING_RE):
            m = rx.match(stripped)
            if m and not re.match(r"^\w+\s+std::", stripped):
                found.append(m.group(1))
                break
        else:
            first = stripped.split(None, 1)[0]
            if first in SKIP_FIRST_WORDS or first.startswith("SPATIUM_"):
                continue
            m = FUNC_RE.match(stripped)
            if m and m.group(1) not in SKIP_FIRST_WORDS:
                found.append(m.group(1))
    seen, out = set(), []
    for n in found:
        if n not in seen and not n.startswith("_"):
            seen.add(n)
            out.append(n)
    return out


def render():
    by_domain = defaultdict(list)
    undescribed = []
    total = 0
    for path in sorted(INCLUDE_ROOT.rglob("*.hpp")):
        rel = path.relative_to(INCLUDE_ROOT).as_posix()
        if rel == "_export_macro.hpp":
            continue
        total += 1
        lines = path.read_text(encoding="utf-8").splitlines()
        desc = description(lines)
        decl = names(lines)
        domain = rel.split("/")[0] if "/" in rel else "(top level)"
        by_domain[domain].append((rel, desc, decl))
        if not desc:
            undescribed.append(rel)

    out = []
    out.append("# Capabilities\n")
    out.append("Generated by `scripts/gen_capabilities.py` from the headers themselves; do not edit by hand. "
               "CI runs it with `--check`. A finder, not a manual: the sentence says what a header is for, "
               "the names say where to look next.\n")
    out.append(f"{total} headers; {len(undescribed)} without a descriptive comment"
               + (": " + ", ".join(f"`{u}`" for u in undescribed) if undescribed else "") + ".\n")
    for domain in sorted(by_domain):
        out.append(f"\n## {domain}\n")
        out.append("| header | what it is for | declares |")
        out.append("|---|---|---|")
        for rel, desc, decl in by_domain[domain]:
            shown = ", ".join(f"`{n}`" for n in decl[:MAX_NAMES])
            if len(decl) > MAX_NAMES:
                shown += f", +{len(decl) - MAX_NAMES} more"
            cell = (desc or "*no description*").replace("|", "\\|")
            out.append(f"| `{rel}` | {cell} | {shown} |")
    return "\n".join(out) + "\n"


def main():
    text = render()
    if "--check" in sys.argv:
        current = OUTPUT.read_text(encoding="utf-8") if OUTPUT.exists() else ""
        if current != text:
            print("docs/capabilities.md is stale; run scripts/gen_capabilities.py", file=sys.stderr)
            return 1
        print("docs/capabilities.md is up to date")
        return 0
    OUTPUT.write_text(text, encoding="utf-8")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
