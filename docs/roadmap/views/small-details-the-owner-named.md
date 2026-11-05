---
date: 2026-11-05
title: the details the owner named: function keys, a link copied in place, clearing, a slow answer
files: [web/src/keys.ts, web/test/keys.test.ts, web/src/core.ts, web/src/view/cells.ts, web/src/view/history.ts, web/src/page.ts, site/index.html, site/gallery/index.html, site/style.css]
open:
  - the keyboard's log with a base fills two placeholders; moving between them is MathLive's Tab
  - the canvas still has a gear per card; the canvas is rebuilt later
---

Five details the owner named, each checked in headless Chromium.

- **Function keys.** In MathLive's own layout, sin and ln were the shifted face of π and e, and lg and log with a base were missing; the key icons shrank on Shift. The keyboard now has one layout of its own (`keys.ts`) with sin, cos, ln, lg, log with a base, √, a square and a power as keys of their own, and no shift layer. A node test reads every function key's call, filled in, through the editor's parser. Pressed on the page, sin then 1 gives ≈ 0.8414709848079, ln then 2 gives ≈ 0.69314718055995, and lg then 100 gives ≈ 2, proved to 15 digits. The keyboard stays open.
- **Copying a link.** The icon turns into a check for a moment and a screen reader hears "Link copied". No tooltip appears and no toast rises under the keyboard; the bar has no `title` tooltips left.
- **Clearing.** Esc twice clears the line, its question going into the history. Settings › Saved work › Forget history and canvas asks once, in place, then forgets both.
- **A slow answer.** A question asked again keeps its last answer, dimmed, until the new one comes, and the core already drops a superseded question. One still unanswered after ten seconds says "Still working", with Stop, which ends the worker and starts a new one.
