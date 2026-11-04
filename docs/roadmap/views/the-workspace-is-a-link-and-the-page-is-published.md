---
date: 2026-11-03
title: the workspace is a link, and the page is published on GitHub Pages
files: [site/index.html]
resolves: [views/the-page-is-windows-of-questions#1, views/the-page-is-windows-of-questions#2]
open:
  - a link opened in a browser replaces the workspace that browser kept; the one it had is not offered back
  - the page is published by hand (site/publish.sh), not by CI, since a workflow cannot be pushed over HTTPS without its scope
---

The address names the workspace: `#w=` and the windows as JSON, deflate-raw (`CompressionStream`), base64url. A window is written as what it holds, an array led by a letter -- `["f", rows...]` a formula, `["t", text]` a text, `["g", id]` a gallery card, `["r"]` the reference, `["G"]` the gallery, `+` after the letter when its graph is expanded -- never an answer, which the core computes again. A window the page cannot open (an unknown letter, a card no longer in the gallery) is left out. A moment after each change the address is rewritten in place (`history.replaceState`, no new history entry) and the same JSON is kept in the browser, so a reload or a new visit opens the workspace as it was left. A page opens, in this order: the workspace a `#w=` names, a question `#solve?q=` names (as text, the old link), the workspace this browser kept, two questions to start from. A link the page cannot read is said to be unreadable and the page goes on to the next. The bar's link button copies the address.

The two questions the page starts from are a link of 82 characters, `#w=` included (their JSON is 69); with the gallery and the reference opened and a graph expanded, 95. Written with keys (`{"kind": "ask", "mode": "formula", "rows": [...]}`) the same two were 141, which is why a window is an array.

Checked in a headless browser on the page as `site/publish.sh --dir` lays it out: a reload opens the four windows as they were, the expanded graph included; `#solve?q=2%20%2B%202` opens one window answered 4 and the address becomes its `#w=`; a `#w=` pasted into an open page replaces its windows, leaving out an unknown card and an unknown letter; `#w=garbage!!` is said to be unreadable and the two questions open; at 400 px the bar wraps and the page does not scroll sideways.

The page is served from the `gh-pages` branch at `vaniell0.github.io/spatium/`. `site/publish.sh` builds the core (`nix develop .#wasm`), lays out what the page loads -- the page, the gallery, the reference, the vendored editor and KaTeX, the core's `.js` and `.wasm` -- and commits that tree on top of the remote's `gh-pages` without a checkout, naming the commit it was built from; `--dir` lays it out only, to look at first. Nothing the page does not load is published: not the core's sources, not the checks.
