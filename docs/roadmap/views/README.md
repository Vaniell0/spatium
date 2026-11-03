# Views: the page

## Where it stands

The page (`site/`), published at `vaniell0.github.io/spatium/` from the `gh-pages` branch by `site/publish.sh`, shows what the core answers and composes nothing: windows of questions asked as rows of a LaTeX editor (MathLive) and answered live by the core in a worker, with their steps, certificates and graphs, the workspace kept in its link (`#w=`) and in the browser, the gallery written by `examples/gallery.cpp`, the reference written from the reader's tables. The core it calls is the library built to WebAssembly (`site/wasm/`), held against the host on every question in `questions.txt`.

## Now

- the page, step two: sliders for the parameters a graph holds at their middles, graphs of answers in pieces, the exact view
- testing the page by many questions, from a generator shared with the synthesis data

## Later

- free windows; on a phone, sheets with a switch
- scenes in a window, from a TypeScript builder, drawn in WebGL2

The entries of this area, by date, and what they did not do: `python3 scripts/gen_roadmap.py --area views`.
