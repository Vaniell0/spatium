#!/usr/bin/env python3
"""rsc_judge: whether a branch does a task, and nothing else changes for the worse.

    python3 rsc/spatium/rsc_judge.py BRANCH --task N [--base main] [--out verdict.json]

A task is an entry of rsc/spatium/tasks.json on the base: what rsc_laws.py found wrong -- a space
breaking laws its neighbours keep, or arithmetic a tangent lacks. Whoever proposes a fix, a person
or an agent, hands in a branch; the judge decides on it, in this order, and stops at the first
step that fails:

  size        the change against the base: files touched, lines added and removed,
              S = touched / (net growth + 1), reported, never a reason to refuse
  no weaker   the tests may grow and may not shrink: every TEST_CASE of the base is still there by
              name, and TEST_CASEs and assertions under tests/ are counted on both sides
  builds      the test binary builds
  tests       every test case passes but the connectivity matrix's (which has its own job in CI)
  the task    against the base's snapshot of the laws, whatever the branch's own says:
              for a law: rsc_laws.py check on the task's space reports the task's law as now
              holding, and no law there stops holding; for an operation: rsc_laws.py missing no
              longer lists it
  no regress  rsc_laws.py check on every space whose headers the change touched reports no law
              that stops holding
  generated   the checks CI runs on what is generated from the tree: roadmap, file references,
              capabilities, consumers

The branch is judged in a work tree of its own, ~/tmp/judge/tree, with a build directory that stays
between verdicts, so a second branch rebuilds only what it changes; nothing in the checkout this
runs from is touched. The verdict is JSON: every step with its outcome, its time and what it saw,
and `accept` true only if all passed. Run inside `nix develop` from the repository root.
"""
import argparse
import json
import pathlib
import re
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[2]
JUDGE = pathlib.Path.home() / "tmp" / "judge"
TREE = JUDGE / "tree"
BUILD = JUDGE / "build"
CPUS = "4-11"   # the E-cores: a verdict is long analysis, not a measured benchmark
JOBS = "4"


def sh(cmd, cwd, input=None, timeout=None):
    r = subprocess.run(cmd, cwd=cwd, input=input, capture_output=True, text=True, timeout=timeout)
    return r.returncode, r.stdout, r.stderr


def git(*args, cwd=ROOT):
    code, out, err = sh(["git", *args], cwd)
    if code:
        raise RuntimeError(f"git {' '.join(args)}: {err.strip()}")
    return out


# ── what can be decided from text alone (test_rsc_judge.py holds these) ─────────

def size_of(numstat):
    """Files touched, lines added and removed, and S, from `git diff --numstat`."""
    touched = added = removed = 0
    for line in numstat.splitlines():
        a, r, _ = line.split("\t", 2)
        touched += 1
        if a != "-":
            added += int(a)
            removed += int(r)
    growth = max(added - removed, 0)
    return {"touched": touched, "added": added, "removed": removed, "S": round(touched / (growth + 1), 4)}


ASSERTION = re.compile(r"\b(CHECK|REQUIRE|CHECK_THAT|REQUIRE_THAT|CHECK_FALSE|REQUIRE_FALSE|CHECK_NOTHROW|REQUIRE_NOTHROW|"
                       r"CHECK_THROWS\w*|REQUIRE_THROWS\w*|STATIC_REQUIRE|STATIC_CHECK)\s*\(")


def strength(sources):
    """TEST_CASEs and assertions in a set of test sources."""
    cases = sum(len(re.findall(r"\bTEST_CASE\s*\(", s)) for s in sources)
    asserts = sum(len(ASSERTION.findall(s)) for s in sources)
    return {"test_cases": cases, "assertions": asserts}


def case_names(sources):
    return {m for s in sources for m in re.findall(r'\bTEST_CASE\s*\(\s*"((?:[^"\\]|\\.)*)"', s)}


def weaker(before_sources, after_sources):
    """What makes the tests weaker: a test case of the base gone (by name), or fewer of either count.
    A count alone can be kept by removing a test and adding an empty one; a name cannot."""
    gone = sorted(case_names(before_sources) - case_names(after_sources))
    b, a = strength(before_sources), strength(after_sources)
    return gone, [k for k in b if a[k] < b[k]], b, a


def law_changes(output):
    """The lines of rsc_laws.py check: (space, law, what) for each change it lists."""
    out = []
    for line in output.splitlines():
        m = re.match(r"\s+(\S+): (.+) (held, now fails|failed, now holds|is now defined|is no longer defined)$", line)
        if m:
            out.append((m[1], m[2], m[3]))
    return out


def task_done(task, changes):
    """For a law task: its law now holds on its space, and no law there stops holding."""
    lost = [f"{law} {what}" for s, law, what in changes if s == task["space"] and what == "held, now fails"]
    if lost:
        return False, "a law stops holding on the task's space: " + lost[0]
    gained = {law for s, law, what in changes if s == task["space"] and what == "failed, now holds"}
    if task["claim"] in gained:
        return True, f"{task['claim']} now holds on {task['space']}, with {len(gained) - 1} more"
    return False, f"{task['claim']} still does not hold on {task['space']}"


def regressions(changes):
    return [f"{s}: {law} {what}" for s, law, what in changes if what in ("held, now fails", "is no longer defined")]


# ── the steps that need the tree ───────────────────────────────────────────────

def prepare(branch):
    JUDGE.mkdir(parents=True, exist_ok=True)
    if not TREE.exists():
        git("worktree", "add", "--detach", str(TREE), branch)
    else:
        git("checkout", "--detach", "--force", branch, cwd=TREE)
        git("clean", "-fdq", cwd=TREE)
    if not (BUILD / "build.ninja").exists():
        code, _, err = sh(["cmake", "-S", str(TREE), "-B", str(BUILD), "-G", "Ninja", "-DCMAKE_BUILD_TYPE=Debug",
                           "-DSPATIUM_BOOST=ON", "-DSPATIUM_BUILD_VIEWER=OFF", "-DSPATIUM_BUILD_EXAMPLES=OFF"], TREE)
        if code:
            raise RuntimeError("cmake: " + err[-2000:])


def test_sources(rev):
    names = git("ls-tree", "-r", "--name-only", rev, "tests/").split()
    return [git("show", f"{rev}:{n}") for n in names if n.endswith((".cpp", ".hpp"))]


def judge(branch, base, task):
    steps = []
    def step(name, ok, seen, t0):
        steps.append({"step": name, "ok": ok, "seen": seen, "seconds": round(time.time() - t0, 1)})
        return ok

    t = time.time()
    merge_base = git("merge-base", base, branch).strip()
    size = size_of(git("diff", "--numstat", merge_base, branch))
    step("size", True, size, t)

    t = time.time()
    gone, fewer, before, after = weaker(test_sources(merge_base), test_sources(branch))
    seen = {"before": before, "after": after}
    if gone:
        seen["test cases gone"] = gone[:10]
    if not step("no weaker", not gone and not fewer, seen, t):
        return steps

    t = time.time()
    prepare(branch)
    code, out, err = sh(["taskset", "-c", CPUS, "ninja", "-C", str(BUILD), "-j", JOBS, "spatium_tests"], TREE)
    if not step("builds", code == 0, (out + err)[-1500:] if code else "spatium_tests", t):
        return steps

    t = time.time()
    code, out, err = sh(["taskset", "-c", CPUS, str(BUILD / "tests" / "spatium_tests"), "~[connectivity]"], TREE)
    tail = [l for l in out.splitlines() if "passed" in l or "failed" in l][-2:]
    if not step("tests", code == 0, tail if code == 0 else (out[-2500:]), t):
        return steps

    # Both law checks compare with the base's snapshot: a branch that does a task also takes its gain
    # with `rsc_laws.py discover`, and against its own snapshot the gain would be no change at all.
    snapshot = TREE / "rsc" / "spatium" / "laws.json"
    snapshot.write_text(git("show", f"{merge_base}:rsc/spatium/laws.json"))
    laws = ["python3", "rsc/spatium/rsc_laws.py", "check"]
    if task["kind"] == "breaks a law of its neighbours":
        t = time.time()
        code, out, err = sh(["taskset", "-c", CPUS, *laws, task["space"]], TREE)
        ok, why = task_done(task, law_changes(out))
        if not step("the task", ok, why, t):
            return steps

    if task["kind"] == "missing operation":
        t = time.time()
        code, out, err = sh(["python3", "rsc/spatium/rsc_laws.py", "missing"], TREE)
        still = code == 0 and {"space": task["space"], "operation": task["claim"]} in json.loads(out or "[]")
        if not step("the task", code == 0 and not still,
                    f"{task['space']} still lacks {task['claim']}" if still else err[-400:] if code else
                    f"{task['space']} has {task['claim']}", t):
            return steps

    t = time.time()
    code, out, err = sh(["taskset", "-c", CPUS, *laws, "--changed", merge_base], TREE)
    lost = regressions(law_changes(out))
    if task["kind"] == "breaks a law of its neighbours":
        # the task's own space was judged above, where gains are the point
        lost = [x for x in lost if not x.startswith(task["space"] + ":")]
    if not step("no regress", not lost, lost or out.strip().splitlines()[:1], t):
        return steps

    git("checkout", "--", "rsc/spatium/laws.json", cwd=TREE)   # the branch's own snapshot back

    t = time.time()
    failed = []
    for cmd in (["python3", "scripts/gen_roadmap.py", "--check"], ["python3", "scripts/check_doc_file_refs.py"],
                ["python3", "scripts/gen_capabilities.py", "--check"], ["python3", "scripts/gen_consumers.py", "--check"]):
        code, out, err = sh(cmd, TREE)
        if code:
            failed.append(f"{cmd[1]}: {(out + err).strip()[-400:]}")
    step("generated", not failed, failed or "roadmap, file references, capabilities, consumers", t)
    return steps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("branch")
    ap.add_argument("--task", type=int, required=True, help="index into the base's rsc/spatium/tasks.json")
    ap.add_argument("--base", default="main")
    ap.add_argument("--out")
    a = ap.parse_args()
    tasks = json.loads(git("show", f"{a.base}:rsc/spatium/tasks.json"))
    task = tasks[a.task]
    steps = judge(a.branch, a.base, task)
    accept = all(s["ok"] for s in steps) and steps[-1]["step"] == "generated"
    first = next((s for s in steps if not s["ok"]), None)
    verdict = {"branch": a.branch, "base": a.base, "task": task, "accept": accept,
               "reason": "every step passed" if accept else f"{first['step']}: {first['seen']}" if first else "stopped",
               "steps": steps}
    text = json.dumps(verdict, indent=1)
    if a.out:
        pathlib.Path(a.out).write_text(text + "\n")
    print(text)
    return 0 if accept else 1


if __name__ == "__main__":
    sys.exit(main())
