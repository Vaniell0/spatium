---
date: 2026-10-08
title: A change pays only for the jobs it can touch
stage: 7
files: [.github/workflows/ci.yml, scripts/ci_changes.py]
open:
  - the matrix is one job of 15 minutes on four cores; splitting it into shards would cut the wall time but not the CPU, and needs the generator to take a shard and the check to merge them
  - the nine builds repeat the same compile in different compilers and option sets on every code change; a cache of compiled units (ccache on the non-module jobs) is not done, and a smaller set on a pull request with the full set on main and at night is a policy to decide
  - the filter of the matrix is a list of directories read off the include closure of probes.hpp by hand; if the cells ever include another directory the nightly run is what notices
  - with a protected branch and required checks, a skipped job must be allowed to count as passed; the branch has no protection yet
---

A pull request used to pay for everything: nine builds that compile the whole test suite in different compilers and option sets (4 to 9 minutes each) and the connectivity matrix (15 minutes, the longest job and so the length of every pull request), about 47 runner-minutes, whatever the change was. `scripts/ci_changes.py` now says what a change needs, and a first job, `changes`, passes it on. Documentation alone (docs/, markdown, images, the scripts that have their own cheap job) skips every build and the matrix; a change to anything that is built runs the builds; and the matrix runs only for what its cells read. Those are core/, algebra/ and spaces/, because the cells include `probes.hpp` and its include closure (read with `g++ -MM`) is nothing else, plus the probes, the generator, the compiler pinned in the flake, and the workflow itself. Of the pull requests merged on 2026-10-07 two were documentation only and would have run nothing heavy; one changed geometry and physics and would have skipped the matrix. A full run happens every night and on request, because filters are a claim about what depends on what, and the night is what notices when the claim is wrong. A diff that cannot be computed (a first push, a force-push) runs everything. The rules hold to examples in `ci_changes.py --selftest`, which the first job runs.
