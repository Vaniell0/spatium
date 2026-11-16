// The canvas's automatic placement: columns by depth, no overlaps, stable output; names read
// and defined by a line.
// @ts-ignore: node's modules, present where the test runs
import { test } from "node:test";
// @ts-ignore
import assert from "node:assert/strict";
import { definedBy, layout, mentioned, nodesOf } from "../src/canvas.js";
import { graph, math, type Doc } from "../src/doc.js";

const torus: Doc = { preset: "canvas", cells: [
  math({ text: "R in [0.5, 4]" }), math({ text: "r in [0.1, 1]" }), math({ text: "h: (r + h)^2 - 2 * r^2 = 0 on [0, 1]" }),
  math({ text: "S0 = torus(R, r)" }), math({ text: "S1 = offset(S0, h)" }) ] };
const quadratic: Doc = { preset: "canvas", cells: [math({ text: "p in [-10, 10]" }), math({ text: "x: 2x^2 - 3x - p = 0" }), graph("x")] };

test("a line's names: defined and read", () => {
  assert.deepEqual(definedBy("x, y: x + 2y = 5, 3x - y = 1"), ["x", "y"]);
  assert.deepEqual(definedBy("I = integral exp(-p * t^2) dt on [0, 1]"), ["I"]);
  assert.deepEqual(definedBy("p in [0.5, 4]"), ["p"]);
  assert.deepEqual(definedBy("2 + 2"), []);
  assert.ok(mentioned("S1 = offset(S0, h)").includes("S0") && mentioned("S1 = offset(S0, h)").includes("h") && !mentioned("S1 = offset(S0, h)").includes("S1"));
});

test("nodes are wired by name to the nearest earlier definition", () => {
  const nodes = nodesOf(torus);
  assert.deepEqual(nodes.map((n) => n.id), ["formula-1", "formula-2", "formula-3", "formula-4", "formula-5"]);
  assert.deepEqual(nodes[2].reads, ["formula-2"]);
  assert.deepEqual(nodes[3].reads.sort(), ["formula-1", "formula-2"]);
  assert.deepEqual(nodes[4].reads.sort(), ["formula-3", "formula-4"]);
  assert.deepEqual(nodesOf(quadratic)[2].reads, ["formula-2"]);
});

const overlaps = (p: { x: number; y: number; w: number; h: number }[]) =>
  p.some((a, i) => p.some((b, j) => i < j && a.x < b.x + b.w && b.x < a.x + a.w && a.y < b.y + b.h && b.y < a.y + a.h));

test("layout: columns by depth, no overlaps, every coordinate on the 8 px grid, stable", () => {
  for (const doc of [torus, quadratic]) {
    const nodes = nodesOf(doc);
    const placed = layout(nodes, () => 120);
    const by = new Map(placed.map((p) => [p.id, p]));
    for (const n of nodes) for (const r of n.reads) assert.ok(by.get(n.id)!.column > by.get(r)!.column && by.get(n.id)!.x > by.get(r)!.x + by.get(r)!.w);
    assert.ok(!overlaps(placed));
    assert.ok(placed.every((p) => p.x % 8 === 0 && p.y % 8 === 0));
    assert.deepEqual(layout(nodes, () => 120), placed);
  }
  const t = layout(nodesOf(torus), () => 120);
  assert.deepEqual(t.map((p) => p.column), [0, 0, 1, 1, 2]);
});

test("layout: a pinned node keeps its place and the rest clear it", () => {
  const nodes = nodesOf(torus);
  nodes[1].at = [48, 48];
  const placed = layout(nodes, () => 120);
  const pinned = placed.find((p) => p.id === "formula-2")!;
  assert.deepEqual([pinned.x, pinned.y], [48, 48]);
  assert.ok(!overlaps(placed));
});
