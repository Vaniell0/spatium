#!/usr/bin/env python3
"""The parts of the judge that need no tree and no compiler: size, strength, what a check reports, a task done.

    python3 rsc/spatium/test_rsc_judge.py
"""
import pathlib
import sys
import unittest

HERE = pathlib.Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
import rsc_judge as rj  # noqa: E402

TASK = {"kind": "breaks a law of its neighbours", "space": "H2", "claim": "d(p,p) = d(exp(p,w),exp(p,w))"}

CHECK = """834 laws checked on 2 spaces (H2, S2xH2): 3 changes
  H2: d(p,p) = d(exp(p,w),exp(p,w)) failed, now holds
  H2: 0 = d(q,q) failed, now holds
  S2xH2: d(p,q) = d(q,p) held, now fails
A law that now fails is a regression; one that now holds is a gain -- run `rsc_laws.py discover` to take it.
"""


class Size(unittest.TestCase):
    def test_s_is_touched_over_net_growth(self):
        self.assertEqual(rj.size_of("10\t2\ta.hpp\n0\t5\tb.cpp\n-\t-\tc.png"),
                         {"touched": 3, "added": 10, "removed": 7, "S": round(3 / 4, 4)})

    def test_a_change_that_shrinks_divides_by_one(self):
        self.assertEqual(rj.size_of("1\t9\ta.hpp")["S"], 1.0)


class Strength(unittest.TestCase):
    def test_cases_and_assertions_are_counted(self):
        src = 'TEST_CASE("a") { CHECK(x); REQUIRE_THAT(y, z); CHECK_FALSE(w); }\nTEST_CASE("b") { STATIC_REQUIRE(t); }'
        self.assertEqual(rj.strength([src]), {"test_cases": 2, "assertions": 4})

    def test_a_word_that_only_contains_check_is_not_one(self):
        self.assertEqual(rj.strength(["bool checked(int); recheck(1);"])["assertions"], 0)


class Laws(unittest.TestCase):
    def test_the_check_output_is_read_line_by_line(self):
        self.assertEqual(rj.law_changes(CHECK), [
            ("H2", "d(p,p) = d(exp(p,w),exp(p,w))", "failed, now holds"),
            ("H2", "0 = d(q,q)", "failed, now holds"),
            ("S2xH2", "d(p,q) = d(q,p)", "held, now fails"),
        ])

    def test_the_task_is_done_when_its_law_now_holds(self):
        ok, why = rj.task_done(TASK, rj.law_changes(CHECK))
        self.assertTrue(ok, why)

    def test_the_task_is_not_done_when_its_law_still_fails(self):
        ok, _ = rj.task_done(TASK, [("H2", "0 = d(q,q)", "failed, now holds")])
        self.assertFalse(ok)

    def test_a_law_lost_on_the_task_space_undoes_the_task(self):
        changes = rj.law_changes(CHECK) + [("H2", "q = exp(p,log(p,q))", "held, now fails")]
        ok, why = rj.task_done(TASK, changes)
        self.assertFalse(ok)
        self.assertIn("stops holding", why)

    def test_a_regression_anywhere_is_listed(self):
        self.assertEqual(rj.regressions(rj.law_changes(CHECK)), ["S2xH2: d(p,q) = d(q,p) held, now fails"])


if __name__ == "__main__":
    unittest.main()
