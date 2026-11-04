#!/usr/bin/env bash
# Lay the page out as it is served: build the core and the page's package (web/), and copy what
# the page loads into <dir>.
#   site/bundle.sh <dir>
# Run from anywhere in the repository; .github/workflows/pages.yml publishes <dir> on GitHub Pages,
# and a local server on <dir> shows the page as it will be.
set -euo pipefail
out=$(realpath -m "${1:?usage: site/bundle.sh <dir>}")
cd "$(dirname "$0")/.."

nix develop .#wasm --command bash -c "site/wasm/build.sh && tsc -p web"

mkdir -p "$out/wasm"
cp site/index.html site/style.css site/core-worker.js "$out/"
mkdir -p "$out/pkg" && cp -r site/pkg/src "$out/pkg/"   # the package, compiled; not its tests
cp -r site/data site/gallery site/vendor "$out/"
cp site/wasm/spatium_core.js site/wasm/spatium_core.wasm "$out/wasm/"
echo "the page is laid out in $out"
