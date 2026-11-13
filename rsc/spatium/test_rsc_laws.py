#!/usr/bin/env python3
"""The parts of rsc_laws that need no compiler: what a snapshot comparison reports, what becomes a task.

    python3 rsc/spatium/test_rsc_laws.py

The evaluation itself is held by `rsc_laws.py check` in CI; this holds the rules applied to its
output, on a snapshot small enough to read.
"""
import pathlib
import sys
import unittest

HERE = pathlib.Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
import rsc_laws as rl  # noqa: E402

SPACES = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]


def law(lhs, rhs, hold, defined, rel="=", size=4):
    return {"lhs": lhs, "rel": rel, "rhs": rhs, "size": size, "hold": hold, "defined": defined}


DATA = {
    "spaces": SPACES, "tuples": 1, "max_size": 4, "terms": 0, "merged": 0, "missing": [],
    "laws": [
        law("d(q,p)", "d(p,q)", "all", "all"),
        law("q", "exp(p,log(p,q))", "1111111110", "1111111111"),   # nine of ten: J is a suspect
        law("u", "log(p,exp(p,u))", "1111100000", "1111111111"),   # half: a class, not a defect
    ],
}


class Changes(unittest.TestCase):
    def test_same_bits_are_no_change(self):
        self.assertEqual(rl.changes_of(DATA, ["B", "J"], ["0\t11\t11", "1\t10\t11", "2\t10\t11"]), [])

    def test_a_law_that_stops_holding_is_reported(self):
        out = rl.changes_of(DATA, ["B"], ["0\t1\t1", "1\t0\t1", "2\t1\t1"])
        self.assertEqual(out, ["B: q = exp(p,log(p,q)) held, now fails"])

    def test_a_law_that_starts_holding_is_reported_too(self):
        out = rl.changes_of(DATA, ["J"], ["0\t1\t1", "1\t1\t1", "2\t0\t1"])
        self.assertEqual(out, ["J: q = exp(p,log(p,q)) failed, now holds"])

    def test_a_law_of_every_space_may_be_undefined_on_one(self):
        self.assertEqual(rl.changes_of(DATA, ["C"], ["0\t0\t0", "1\t1\t1", "2\t1\t1"]), [])

    def test_a_law_of_every_space_failing_where_defined_is_reported(self):
        out = rl.changes_of(DATA, ["C"], ["0\t0\t1", "1\t1\t1", "2\t1\t1"])
        self.assertEqual(out, ["C: d(q,p) = d(p,q) held, now fails"])


class Tasks(unittest.TestCase):
    def test_a_space_breaking_what_its_neighbours_keep_is_a_task(self):
        tasks = rl.tasks_of(DATA)
        self.assertEqual([(t["space"], t["claim"]) for t in tasks], [("J", "q = exp(p,log(p,q))")])

    def test_a_law_of_half_the_spaces_draws_a_class_and_no_task(self):
        tasks = rl.tasks_of({**DATA, "laws": DATA["laws"][2:]})
        self.assertEqual(tasks, [])

    def test_a_missing_operation_is_a_task(self):
        tasks = rl.tasks_of({**DATA, "laws": [], "missing": [{"space": "B", "operation": "unary minus on a tangent"}]})
        self.assertEqual([(t["kind"], t["space"]) for t in tasks], [("missing operation", "B")])


if __name__ == "__main__":
    unittest.main()
