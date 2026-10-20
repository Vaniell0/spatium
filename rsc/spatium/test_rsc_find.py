#!/usr/bin/env python3
"""The parts of the finder that need no compiler: what a route becomes, which is the reference, what is independent.

    python3 rsc/spatium/test_rsc_find.py

The run itself (compile, execute, compare) is held by `rsc_find.py --all --check` in CI; this holds the
rules it applies, on a graph small enough to read.
"""
import json
import pathlib
import sys
import tempfile
import unittest

HERE = pathlib.Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
import emit            # noqa: E402
import rsc_find as rf  # noqa: E402

# a type with a closed-form distance, a derivation, and the same over Real50; a product of two such
GRAPH = {
    "version": 1,
    "nodes": [
        {"id": "S", "kind": "type", "cell": "connectivity::S2"},
        {"id": "S@Real50", "kind": "type", "cell": "connectivity::S2", "lift_of": "S", "scalar": "Real50"},
        {"id": "E", "kind": "type", "cell": "connectivity::E3"},
        {"id": "P", "kind": "type", "cell": "connectivity::S2xE1", "factors": ["S", "E"]},
    ],
    "edges": [
        {"kind": "operation", "op": "distance", "variant": "member", "output": "distance", "inputs": ["S"], "needs": [], "level": 3},
        {"kind": "operation", "op": "log_map", "variant": "member", "output": "log_map", "inputs": ["S"], "needs": [], "level": 3},
        {"kind": "operation", "op": "metric", "variant": "member", "output": "metric", "inputs": ["S"], "needs": [], "level": 3},
        {"kind": "operation", "op": "distance", "variant": "derived", "output": "distance", "inputs": ["S"], "needs": ["log_map", "metric"], "level": 3},
        {"kind": "operation", "op": "distance", "variant": "member", "output": "distance", "inputs": ["S@Real50"], "needs": [], "level": 3},
        {"kind": "operation", "op": "distance", "variant": "lift:Real50", "output": "distance", "inputs": ["S"], "needs": ["S@Real50#distance"], "level": 3, "entry": True},
        {"kind": "operation", "op": "distance", "variant": "member", "output": "distance", "inputs": ["E"], "needs": [], "level": 3},
        {"kind": "operation", "op": "distance", "variant": "member", "output": "distance", "inputs": ["P"], "needs": [], "level": 3},
        {"kind": "operation", "op": "distance", "variant": "components", "output": "distance", "inputs": ["P"], "needs": ["S#distance", "E#distance"], "level": 0},
    ],
}


def graph():
    d = tempfile.mkdtemp()
    p = pathlib.Path(d) / "graph.json"
    p.write_text(json.dumps(GRAPH))
    return emit.Graph(p)


class Emission(unittest.TestCase):
    def test_a_derived_route_calls_its_needs_in_order(self):
        g = graph()
        code = emit.emit_route(g, [1, 2, 3])
        self.assertIn("inline auto e0(const auto& s, const auto&... a) { ++rsc_pair::calls(); return rsc_pair::val(s.log_map(a...)); }", code)
        self.assertIn("e2(const auto& s, const auto& p, const auto& q)", code)
        self.assertIn("const auto v = e0(s, p, q); return rsc_pair::val(sqrt(e1(s, p, v, v)));", code)
        self.assertEqual(emit.entry_name([1, 2, 3]), "e2")

    def test_a_lift_converts_up_and_back_and_counts_no_call(self):
        g = graph()
        code = emit.emit_route(g, [4, 5])
        lift = [l for l in code.splitlines() if l.startswith("inline auto e1")][0]
        self.assertNotIn("calls()", lift)
        self.assertIn("using R = spatium::Real50;", lift)
        self.assertIn("connectivity::S2::make<R>()", lift)
        self.assertIn("rsc_pair::down(e0(sl, rsc_pair::up_pt<R>(sl, p), rsc_pair::up_pt<R>(sl, q)))", lift)

    def test_a_product_calls_its_factors_on_their_own_spaces(self):
        g = graph()
        code = emit.emit_route(g, [0, 6, 8])
        comp = [l for l in code.splitlines() if l.startswith("inline auto e2")][0]
        self.assertIn("e0(s.space1,", comp)
        self.assertIn("e1(s.space2,", comp)
        self.assertIn("sqrt(d1 * d1 + d2 * d2)", comp)


class Rules(unittest.TestCase):
    def test_a_lift_is_not_a_node_but_what_it_computes_is(self):
        g = graph()
        self.assertEqual(g.nodes_of([4, 5]), frozenset({("S", "distance", "member")}))   # the same chain as the closed form
        self.assertEqual(g.nodes_of([4, 5]), g.nodes_of([0]))

    def test_independent_means_no_shared_node(self):
        g = graph()
        derived, comps = g.nodes_of([1, 2, 3]), g.nodes_of([0, 6, 8])
        member = g.nodes_of([7])
        self.assertFalse(derived & member)
        self.assertFalse(comps & member)
        self.assertTrue(comps & g.nodes_of([0]))      # the product from its factors uses S's closed form

    def test_the_reference_is_a_closed_form_at_l3_else_the_real50_lift(self):
        g = graph()
        k, why = rf.choose_reference(g, [[3], [0], [4, 5]])
        self.assertEqual(k, 1)
        self.assertIn("closed form", why)
        k, why = rf.choose_reference(g, [[1, 2, 3], [4, 5]])
        self.assertEqual(k, 1)
        self.assertEqual(why, "Real50 lift")
        k, why = rf.choose_reference(g, [[1, 2, 3], [1, 2, 3, 0]])
        self.assertEqual(k, 0)
        self.assertIn("shortest", why)

    def test_an_edge_with_no_emission_is_refused(self):
        g = graph()
        g.edges.append({"kind": "operation", "op": "distance", "variant": "teleport", "output": "distance", "inputs": ["S"], "needs": [], "level": 0})
        with self.assertRaises(ValueError):
            emit.emit_route(g, [len(g.edges) - 1])


if __name__ == "__main__":
    unittest.main()
