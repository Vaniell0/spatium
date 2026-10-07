#!/usr/bin/env python3
"""Generate docs/capabilities/: what the library can do, found from the code.

For every header under include/spatium/ it records the first sentence of its
own descriptive comment and the names it declares at namespace level (types,
concepts, free functions), one page per domain. Nothing here is written by
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
OUT_DIR = ROOT / "docs" / "capabilities"

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


def page_name(domain):
    return ("top-level" if domain == "(top level)" else domain.replace("/", "-")) + ".md"


def render():
    """{file name: text}: one page per domain, and an index of links. Per-domain pages, not one
    file, so that two changes in different domains never edit the same lines (one shared
    generated table conflicted at every merge)."""
    by_domain = defaultdict(list)
    for path in sorted(INCLUDE_ROOT.rglob("*.hpp")):
        rel = path.relative_to(INCLUDE_ROOT).as_posix()
        if rel == "_export_macro.hpp":
            continue
        lines = path.read_text(encoding="utf-8").splitlines()
        domain = rel.split("/")[0] if "/" in rel else "(top level)"
        by_domain[domain].append((rel, description(lines), names(lines)))

    pages = {}
    for domain, rows in by_domain.items():
        undescribed = [r for r, d, _ in rows if not d]
        out = [f"# Capabilities: `{domain}`\n",
               "Generated by `scripts/gen_capabilities.py` from the headers themselves; do not edit by hand. "
               "A finder, not a manual: the sentence says what a header is for, the names say where to look next.\n",
               f"{len(rows)} headers; {len(undescribed)} without a descriptive comment"
               + (": " + ", ".join(f"`{u}`" for u in undescribed) if undescribed else "") + ".\n",
               "| header | what it is for | declares |", "|---|---|---|"]
        for rel, desc, decl in rows:
            shown = ", ".join(f"`{n}`" for n in decl[:MAX_NAMES])
            if len(decl) > MAX_NAMES:
                shown += f", +{len(decl) - MAX_NAMES} more"
            cell = (desc or "*no description*").replace("|", "\\|")
            out.append(f"| `{rel}` | {cell} | {shown} |")
        pages[page_name(domain)] = "\n".join(out) + "\n"
    index = ["# Capabilities\n",
             "What the library can do, found from the headers and regenerated by `scripts/gen_capabilities.py` "
             "(CI runs it with `--check`). One page per domain:\n"]
    for domain in sorted(by_domain):
        index.append(f"- [`{domain}`]({page_name(domain)})")
    pages["README.md"] = "\n".join(index) + "\n"
    return pages


def main():
    pages = render()
    if "--check" in sys.argv:
        bad = []
        for name, text in pages.items():
            f = OUT_DIR / name
            if not f.exists() or f.read_text(encoding="utf-8") != text:
                bad.append(f"docs/capabilities/{name} is stale")
        if OUT_DIR.exists():
            for f in OUT_DIR.glob("*.md"):
                if f.name not in pages:
                    bad.append(f"docs/capabilities/{f.name} is not generated any more")
        for b in bad:
            print(b, file=sys.stderr)
        if bad:
            print("run scripts/gen_capabilities.py", file=sys.stderr)
            return 1
        print(f"docs/capabilities: {len(pages) - 1} pages are up to date")
        return 0
    OUT_DIR.mkdir(exist_ok=True)
    for f in OUT_DIR.glob("*.md"):
        if f.name not in pages:
            f.unlink()
    for name, text in pages.items():
        (OUT_DIR / name).write_text(text, encoding="utf-8")
    print(f"wrote docs/capabilities/ ({len(pages)} files)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
