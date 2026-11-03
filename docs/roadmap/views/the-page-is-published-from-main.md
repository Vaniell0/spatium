---
date: 2026-11-03
title: the page is published on GitHub Pages from each commit of main CI passed
files: [.github/workflows/pages.yml, site/bundle.sh]
resolves: [views/the-workspace-is-a-link-and-the-page-is-published#2]
---

GitHub Pages builds from GitHub Actions: `.github/workflows/pages.yml` runs when CI completes on a push to main and passed (or by hand), checks out the commit CI ran on, lays the page out with `site/bundle.sh _site` -- the core built in `nix develop .#wasm`, the page, the gallery, the reference, the vendored editor and KaTeX -- and deploys it (`upload-pages-artifact`, `deploy-pages`). A commit CI failed is not published, and nothing is committed to a branch for it.

`site/bundle.sh <dir>` is the layout step of the script it replaces, `site/publish.sh`, without the commit on a `gh-pages` branch: a local server on its directory shows the page as it will be served. Its output is the layout the link was checked on in the browser, file for file, less `.nojekyll`, which a page built by Actions does not need.
