---
date: 2026-11-03
title: the roadmap is kept by area, and each area says what is done, what is taken up now and what would be good later
files: [scripts/gen_roadmap.py, docs/ROADMAP.md, docs/roadmap/tooling/README.md]
---

The roadmap grew as files named by date, `docs/roadmap/<date>-<nn>-<slug>.md`, with 21 of them on one day. Beside them was a `ROADMAP.md` of 1794 lines, its history in date order and its reasoning in topics. What an area had done, and what it had left, could be read only by reading everything.

It is now kept by area: ir, document, evaluation, lowering, views, rsc, algebra, spaces, geometry, physics, relativity, discrete, sound, tooling, and history for the steering index of September 2026. Each area is a folder whose README says where it stands, what is taken up now and what would be good later, above its reasoning and earlier history, which moved there from `ROADMAP.md` as it was written. An entry is a file named for what it did, and its date is in its front matter. A change still adds one file, so two changes in flight still never edit the same lines.

`ROADMAP.md` is now the layer map and a table of the areas, 58 lines. `gen_roadmap.py --area <area>` prints an area's entries by date and its open items, and `--open` groups everything not done by area. An open item is now named `<area>/<slug>#k`, and every `resolves` was rewritten to it. All 91 entries and their 236 open items moved, and the 1444 file references in the docs still resolve.
