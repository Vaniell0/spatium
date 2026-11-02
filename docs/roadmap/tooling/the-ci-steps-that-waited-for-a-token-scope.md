---
date: 2026-10-07
title: The CI steps that waited for a token scope
stage: cross-cutting
files: [scripts/gen_roadmap.py, scripts/gen_capabilities.py]
---
Two pieces of CI were left unfinished while SSH to GitHub was unreachable from the mobile connection the work was done on: a push over HTTPS with the `gh` token is refused a change to a workflow file for want of the `workflow` scope. With SSH back they are made. A `roadmap` job runs `scripts/gen_roadmap.py --check` (front matter, dates, listed files, closed findings rows), and the display name of the capabilities step, which still named the single file `docs/capabilities.md` that the per-domain pages replaced, now names the directory; the step ran the right command throughout.
