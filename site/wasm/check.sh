#!/usr/bin/env bash
# Build the core both ways and hold the WebAssembly one against the host one.
# Run inside `nix develop .#wasm` from the repository root.
set -euo pipefail
cd "$(dirname "$0")/../.."
site/wasm/build.sh
out="${TMPDIR:-/tmp}/spatium_core_host"
c++ -std=c++23 -O2 -I include -I . site/wasm/core.cpp -o "$out"
node site/wasm/check.mjs "$out"
