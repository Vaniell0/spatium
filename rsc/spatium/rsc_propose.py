#!/usr/bin/env python3
"""rsc_propose: an agent hands in a branch for a task, and the judge decides on it.

    python3 rsc/spatium/rsc_propose.py --task N [--base main] [--model sonnet] [--budget 5]

The task is an entry of the base's rsc/spatium/tasks.json. The proposer is Claude Code run headless
in a work tree of its own, ~/tmp/propose/<task>, on a new branch task/<task> from the base; it is
told the task, how to reproduce it, what the judge will check, and that it may commit on that
branch and nowhere else. It cannot push: pushing and the GitHub CLI are refused to it, and nothing
here pushes. When it is done, a branch with commits goes to rsc_judge.py, and the verdict is written
beside the work tree with what the proposer said and what it cost. Run from the repository root.
"""
import argparse
import json
import pathlib
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[2]
WORK = pathlib.Path.home() / "tmp" / "propose"

PROMPT = """You are working in a git work tree of the spatium library ({tree}), on the branch {branch},
made from {base} for one task. Read CLAUDE.md at the root of this tree before you change anything.

The task, found by rsc/spatium/rsc_laws.py (laws of the library's spaces found by enumeration):

{task}

What it means: {meaning}

Reproduce it first:
    nix develop -c python3 rsc/spatium/rsc_laws.py {reproduce}

Then find the cause in the library, fix it there, and add a test under tests/ that fails before your
fix and passes after it. Build and run the tests on the E-cores, one build at a time:
    nix develop -c cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug -DSPATIUM_BOOST=ON -DSPATIUM_BUILD_VIEWER=OFF -DSPATIUM_BUILD_EXAMPLES=OFF
    nix develop -c taskset -c 4-11 ninja -C build -j4 spatium_tests
    taskset -c 4-11 build/tests/spatium_tests "[the tag of your test]"
When the fix works, take its gain into the snapshot of the laws:
    nix develop -c taskset -c 4-11 python3 rsc/spatium/rsc_laws.py discover
and write a roadmap entry under docs/roadmap/<area>/ in the form of the others there.

A judge will decide on your branch, in this order: the tests are not weaker (every test case of the
base still there by name, no fewer assertions); it builds; every test passes; the task's law now
holds on its space (or its operation is no longer missing), checked against the base's snapshot; no
law stops holding on a space your change touched; the roadmap, file references, capabilities and
consumers checks pass (scripts/gen_roadmap.py --check, scripts/check_doc_file_refs.py,
scripts/gen_capabilities.py --check, scripts/gen_consumers.py --check).

Commit your work on {branch} with git add and git commit, in English, the message saying what the
code does and why, with the numbers you measured, ending with a blank line and
    Claude-Session: proposer for task {slug}
You may commit on this branch only. Do not push, do not use gh, do not touch other branches or other
work trees, do not write outside {tree} except under ~/.cache. If you cannot fix it, say why in your
last message and commit nothing. Your last message is what the judge's report shows of you: say what
the cause was, what you changed, and what you measured."""

MEANING = {
    "breaks a law of its neighbours": "the law holds on at least 80% of the spaces where it can be asked, and not on this one",
    "missing operation": "the tangent type of this space lacks this arithmetic",
}


def git(*args, cwd=ROOT):
    r = subprocess.run(["git", *args], cwd=cwd, capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f"git {' '.join(args)}: {r.stderr.strip()}")
    return r.stdout


def slug_of(task, index):
    return f"{task['space'].lower().replace('+', 'x')}-{index}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", type=int, required=True)
    ap.add_argument("--base", default="main")
    ap.add_argument("--model", default="sonnet")
    ap.add_argument("--budget", type=float, default=5.0, help="dollars of API the proposer may spend")
    ap.add_argument("--timeout", type=int, default=3600, help="seconds the proposer may take")
    a = ap.parse_args()

    task = json.loads(git("show", f"{a.base}:rsc/spatium/tasks.json"))[a.task]
    if "triaged" in task:
        sys.exit(f"task {a.task} is triaged, no defect: {task['triaged']}")
    slug = slug_of(task, a.task)
    branch = f"task/{slug}"
    tree = WORK / slug
    WORK.mkdir(parents=True, exist_ok=True)
    if tree.exists():
        sys.exit(f"{tree} exists: a proposal for this task was made; remove it with `git worktree remove` first")
    git("worktree", "add", "-b", branch, str(tree), a.base)

    reproduce = f"check {task['space']}" if task["kind"] == "breaks a law of its neighbours" else "missing"
    prompt = PROMPT.format(tree=tree, branch=branch, base=a.base, slug=slug, task=json.dumps(task, indent=1),
                           meaning=MEANING.get(task["kind"], task["kind"]), reproduce=reproduce)
    cmd = ["claude", "-p", prompt, "--model", a.model, "--output-format", "json",
           "--permission-mode", "acceptEdits",
           "--allowedTools", "Read", "Edit", "Write", "Grep", "Glob",
           "Bash(nix develop:*)", "Bash(taskset:*)", "Bash(git add:*)", "Bash(git commit:*)", "Bash(git diff:*)",
           "Bash(git status:*)", "Bash(git log:*)", "Bash(python3:*)", "Bash(build/tests/spatium_tests:*)",
           "--disallowedTools", "Bash(git push:*)", "Bash(gh:*)", "Bash(git checkout:*)", "Bash(git switch:*)",
           "Bash(git worktree:*)", "Bash(git reset:*)", "Bash(git rebase:*)",
           "--settings", json.dumps({"disableAllHooks": True}),
           "--strict-mcp-config", "--no-session-persistence", "--max-budget-usd", str(a.budget)]
    t0 = time.time()
    try:
        r = subprocess.run(cmd, cwd=tree, capture_output=True, text=True, timeout=a.timeout)
        said, code = r.stdout, r.returncode
    except subprocess.TimeoutExpired as e:
        said, code = (e.stdout or b"").decode() if isinstance(e.stdout, bytes) else (e.stdout or ""), "timeout"
    proposer = {"seconds": round(time.time() - t0), "exit": code}
    try:
        out = json.loads(said)
        proposer.update({k: out.get(k) for k in ("result", "total_cost_usd", "num_turns", "is_error")})
    except (json.JSONDecodeError, TypeError):
        proposer["result"] = said[-3000:]

    commits = git("log", "--oneline", f"{a.base}..{branch}").strip().splitlines()
    report = {"task": task, "index": a.task, "branch": branch, "proposer": proposer, "commits": commits}
    if commits:
        # inside the flake's shell, as every build here is: outside it the tests do not find Catch2
        j = subprocess.run(["nix", "develop", "-c", "python3", str(ROOT / "rsc" / "spatium" / "rsc_judge.py"), branch,
                            "--task", str(a.task), "--base", a.base], cwd=ROOT, capture_output=True, text=True)
        try:
            report["verdict"] = json.loads(j.stdout)
        except json.JSONDecodeError:
            report["verdict"] = {"accept": False, "reason": "the judge failed: " + (j.stdout + j.stderr)[-2000:]}
    else:
        report["verdict"] = {"accept": False, "reason": "no commit: nothing to judge"}
    path = WORK / f"{slug}.json"
    path.write_text(json.dumps(report, indent=1) + "\n")
    print(f"{branch}: {'accepted' if report['verdict']['accept'] else 'refused'} -- {report['verdict']['reason']}")
    print(f"report: {path}")
    return 0 if report["verdict"]["accept"] else 1


if __name__ == "__main__":
    sys.exit(main())
