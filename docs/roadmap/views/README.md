# Views: the page

## Where it stands

The page (`site/`), published at `vaniell0.github.io/spatium/` from each commit of main CI passed, shows what the core answers and composes nothing: cards, each one question written as the rows of a LaTeX editor (MathLive) with its own settings, answered live by the core in a worker, with their steps, certificates and graphs, the workspace kept in its link (`#w=`) and in the browser, the gallery written by `examples/gallery.cpp`, the reference written from the reader's tables. The core it calls is the library built to WebAssembly (`site/wasm/`), held against the host on every question in `questions.txt`.

## Now

- the page as a beta: free cells on the canvas with wires, a scene cell from the TypeScript builder, the gallery of renders
- testing the page by many questions, from a generator shared with the synthesis data

## Later

- free windows; on a phone, sheets with a switch
- scenes in a window, from a TypeScript builder, drawn in WebGL2

The entries of this area, by date, and what they did not do: `python3 scripts/gen_roadmap.py --area views`.
