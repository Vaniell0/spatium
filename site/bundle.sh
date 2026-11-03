#!/usr/bin/env bash
# Lay the page out as it is served: build the core and copy what the page loads into <dir>.
#   site/bundle.sh <dir>
# Run from anywhere in the repository; .github/workflows/pages.yml publishes <dir> on GitHub Pages,
# and a local server on <dir> shows the page as it will be.
set -euo pipefail
out=$(realpath -m "${1:?usage: site/bundle.sh <dir>}")
cd "$(dirname "$0")/.."

nix develop .#wasm --command site/wasm/build.sh

mkdir -p "$out/wasm"
cp site/index.html site/style.css site/app.js site/core-worker.js "$out/"
cp -r site/data site/gallery site/vendor "$out/"
cp site/wasm/spatium_core.js site/wasm/spatium_core.wasm "$out/wasm/"
echo "the page is laid out in $out"
