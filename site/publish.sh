#!/usr/bin/env bash
# Publish the page to the gh-pages branch, which GitHub Pages serves at vaniell0.github.io/spatium/.
# From a clean checkout at the repository root:
#   site/publish.sh [remote]     build the core, commit what the page loads on gh-pages, push
#   site/publish.sh --dir <dir>  build the core and lay the page out in <dir> only, to look at first
# The commit is written without a checkout (a tree from a temporary index) on top of the remote's
# gh-pages, and says which commit of the source it was built from.
set -euo pipefail
cd "$(dirname "$0")/.."

dir="" remote="origin"
if [[ ${1:-} == --dir ]]; then dir=${2:?"--dir needs a directory"}; else remote=${1:-origin}; fi
if [[ -z $dir && -n $(git status --porcelain) ]]; then
    echo "the tree has changes; publish from a clean checkout" >&2
    exit 1
fi

nix develop .#wasm --command site/wasm/build.sh

out=${dir:-$(mktemp -d)}
[[ -n $dir ]] || trap 'rm -rf "$out"' EXIT
mkdir -p "$out/wasm"
cp site/index.html site/style.css site/app.js site/core-worker.js "$out/"
cp -r site/data site/gallery site/vendor "$out/"
cp site/wasm/spatium_core.js site/wasm/spatium_core.wasm "$out/wasm/"
touch "$out/.nojekyll"  # served as the files are, with no Jekyll build
[[ -z $dir ]] || { echo "the page is laid out in $dir"; exit 0; }

source=$(git rev-parse --short HEAD)
parent=$(git ls-remote "$remote" refs/heads/gh-pages | cut -f1)
[[ -z $parent ]] || git fetch -q "$remote" "$parent"
scratch=$(mktemp -d)
trap 'rm -rf "$out" "$scratch"' EXIT
export GIT_INDEX_FILE=$scratch/index
git --git-dir="$(git rev-parse --absolute-git-dir)" --work-tree="$out" add -A -f "$out"
tree=$(git write-tree)
unset GIT_INDEX_FILE
if [[ -n $parent && $(git rev-parse "$parent^{tree}") == "$tree" ]]; then
    echo "gh-pages already holds this page (site from $source)"
    exit 0
fi
commit=$(git commit-tree "$tree" ${parent:+-p "$parent"} -m "site from $source")
git push "$remote" "$commit:refs/heads/gh-pages"
echo "published $commit (site from $source)"
