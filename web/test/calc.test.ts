// The calculator's history: a question asked again moves up instead of being kept twice, spacing
// does not make two questions different, a blank line is never kept, and it holds at most MOST.
// @ts-ignore: node's modules, present where the test runs
import { test } from "node:test";
// @ts-ignore
import assert from "node:assert/strict";
import { type Calc, type Entry, MOST, commit, empty } from "../src/calc.js";

const entry = (rows: string[], t = 0): Entry => ({ written: { rows }, result: "x", strength: "exact", note: "", ok: true, t });
const start = (): Calc => ({ line: empty(), history: [] });

test("the newest is on top, and one asked again moves up", () => {
  let c = commit(start(), entry(["1+1"], 1));
  c = commit(c, entry(["2x^2-3x-5=0"], 2));
  c = commit(c, entry(["1 + 1"], 3));
  assert.deepEqual(c.history.map((e) => e.t), [3, 2]);
});

test("a system is one question, its rows in order", () => {
  let c = commit(start(), entry(["x+2y=5", "3x-y=1"]));
  c = commit(c, entry(["3x-y=1", "x+2y=5"]));
  assert.equal(c.history.length, 2);
});

test("a blank line is not kept, and the history holds at most MOST", () => {
  assert.equal(commit(start(), entry(["", "  "])).history.length, 0);
  let c = start();
  for (let k = 0; k < MOST + 5; ++k) c = commit(c, entry([`${k}`], k));
  assert.equal(c.history.length, MOST);
  assert.equal(c.history[0].t, MOST + 4);
});
