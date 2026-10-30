#!/usr/bin/env bash
# Build the core the page calls to WebAssembly: site/wasm/spatium_core.js and .wasm.
# Run inside `nix develop .#wasm` from the repository root. libc++ keeps std::move_only_function
# (which the field IR holds its closures in) behind -fexperimental-library.
set -euo pipefail
cd "$(dirname "$0")/../.."
em++ -std=c++23 -O2 -fexperimental-library -Wno-inline-namespace-reopened-noninline -I include -I . site/wasm/core.cpp -o site/wasm/spatium_core.js \
    -sMODULARIZE=1 -sEXPORT_NAME=SpatiumCore -sENVIRONMENT=web,node \
    -sALLOW_MEMORY_GROWTH=1 \
    -sEXPORTED_FUNCTIONS=_spatium_respond,_spatium_preview \
    -sEXPORTED_RUNTIME_METHODS=cwrap
ls -l site/wasm/spatium_core.js site/wasm/spatium_core.wasm
