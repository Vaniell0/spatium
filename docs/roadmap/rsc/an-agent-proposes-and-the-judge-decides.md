---
date: 2026-11-13
title: An agent proposes a branch for a task, and the judge decides on it
files:
  - rsc/spatium/rsc_propose.py
open:
  - an attempt is two builds from cold, the proposer's and the judge's, about 25 of its 50 minutes; handing the proposer a built tree would make both incremental
  - one attempt per task and no retry with the judge's reason fed back
  - the proposer is a general model told the task in words; nothing learned from the verdicts yet
---

`rsc_propose.py --task N` takes a task of `rsc/spatium/tasks.json` and gives it to Claude Code run headless in a work tree of its own on a new branch, with the task, how to reproduce it, what the judge will check and the commands to build and test on the E-cores. The proposer may read, edit, build, run the law check and commit on its branch; pushing, the GitHub CLI, switching branches and the user's hooks are refused it, and a budget and a timeout bound it. When it stops, a branch with commits goes to `rsc_judge.py` inside the flake's shell, and the report -- the proposer's last message, its cost and turns, the commits and the verdict -- is written beside the work tree. Nothing it does leaves the machine.

The first attempt took the law check's task on H2, d(p,p) = d(exp(p,w), exp(p,w)): the proposer worked 1913 s for $0.71 and committed a distance computed as cosh d - 1 by Lagrange's identity, two tests, a roadmap entry and the law snapshot rediscovered (834 laws to 637); the judge accepted it on every step in 1010 s. A first run of the judge from this script failed at the build, outside the flake's shell where the tests find no Catch2, which is why it is called inside it.
