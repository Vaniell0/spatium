---
date: 2026-10-08
title: A job that waits for a mirror has a time limit
stage: 7
files: [.github/workflows/ci.yml]
open:
  - nothing yet runs the heavy checks only when the files they read changed (the connectivity matrix rebuilds 1168 cells whatever the diff), so every pull request pays the whole CI
  - the apt toolchain installs are repeated in six jobs and not cached
---

On 2026-10-07 the job "cmake build without Nix (gcc15)" stopped answering inside `apt-get update` (the last line of its log is a `Get:5 ... noble-security InRelease`) and sat there until GitHub cancelled it after six hours, on the pull request that held the registry of pairs and, an hour earlier, on the push to main of the NaN fix. Neither was the code: the same job had passed in 15 minutes on the pull request of that fix. Nothing in the workflow had a time limit, so the longest wait that GitHub allows was the only one, and a pull request waited on it with the laptop idle. Every job now has a `timeout-minutes` well above the longest run seen (10 for the documentation checks, 30 for the consumer builds, 60 for the builds and the matrix), and `apt-get update` retries three times with a 30 second timeout on its connections. A rerun of the failed job is the answer to a cancellation; this makes the cancellation come in minutes.
