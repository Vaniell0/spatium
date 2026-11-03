---
date: 2026-10-07
title: The capability index per domain, and CLAUDE.md as a map
stage: cross-cutting
files: [scripts/gen_capabilities.py, CLAUDE.md, docs/capabilities/README.md]
open:
  - the CI step's display name still says docs/capabilities.md
---
Friction of the same kind as the roadmap's, found by the same day of merging: after each merge every other open pull request conflicted again, in two places -- the one generated table `docs/capabilities.md`, where two changes adding a header in neighbouring rows edit adjacent lines, and the directory bullets of `CLAUDE.md`, one line of up to 3.4 KB per directory that every change that added a header to `spaces/` or `algebra/` edited. The first is split by domain: `scripts/gen_capabilities.py` now writes `docs/capabilities/<domain>.md` and an index of links, so changes in different domains never touch the same file (`--check` is unchanged in CI and now also fails on a page the generator no longer writes). The second is a fact about what `CLAUDE.md` is for: it is read in every session and had grown to 33 KB of per-header detail that the finder now carries, generated from the headers and so never stale; its directory section is a map again, one or two lines per directory with a pointer to the finder, 15 KB in all, the design rules kept where they were. Not done: the CI step's display name still says `docs/capabilities.md` (changing a workflow file needs the `workflow` scope of the push token, which the HTTPS route lacks); the step runs the right command.
