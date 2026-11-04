// The document and its link, held two ways: a link is read back by this package and, apart from
// it, inflated by node's zlib, which must give the same JSON.
// @ts-ignore: node's modules, present where the test runs
import { test } from "node:test";
// @ts-ignore
import assert from "node:assert/strict";
// @ts-ignore
import { inflateRawSync } from "node:zlib";
// node's Buffer, with no type package for node in the tree
declare const Buffer: { from(s: string, encoding: string): Uint8Array };
import { type Doc, fromJSON, fromLink, graph, math, toJSON, toLink, setting, splitLines } from "../src/doc.js";

const doc: Doc = {
  preset: "canvas",
  cells: [
    math({ latex: "p\\in[0.5,4]" }),
    math({ latex: "x^2-p=0,\\ x\\in[0,10]" }, { numbers: "exact", digits: 17 }),
    math({ text: "2x^2 - 3x - 5 = 0" }),
    graph("x", { window: [-2, -1, 5, 3] }),
  ],
};

test("a document is its JSON read back", () => {
  assert.deepEqual(fromJSON(toJSON(doc)), doc);
});

test("a link reads back to the document, and zlib reads it as the same JSON", async () => {
  const code = await toLink(doc);
  assert.match(code, /^[A-Za-z0-9_-]+$/);
  assert.deepEqual(await fromLink(code), doc);
  const inflated = inflateRawSync(Buffer.from(code, "base64url")).toString("utf8");
  assert.deepEqual(JSON.parse(inflated), toJSON(doc));
});

test("a link of the windows page is read as its cells", () => {
  const old = [1, ["f", "2x^2-3x-5=0"], ["f+", "p\\in[0.5,4]", "x^2-p=0"], ["t", "2 + 2"], ["r"], ["G"], ["g", "scene"]];
  assert.deepEqual(fromJSON(old), {
    preset: "canvas",
    cells: [math({ latex: "2x^2-3x-5=0" }), math({ latex: "p\\in[0.5,4]" }), math({ latex: "x^2-p=0" }), math({ text: "2 + 2" })],
  });
});

test("what a link cannot say is left out, and a link that is not one is refused", async () => {
  const read = fromJSON([2, "calculator", ["l", "1+1", { digits: 7, numbers: "exact", window: [3, 0, 1, 1], steps: "yes" }], ["q", "?"], 5]);
  assert.deepEqual(read, { preset: "calculator", cells: [math({ latex: "1+1" }, { numbers: "exact" })] });
  assert.equal(fromJSON([3, "canvas"]), null);
  assert.equal(fromJSON({ cells: [] }), null);
  assert.equal(await fromLink("garbage!!"), null);
});

test("a setting left out is the preset's", () => {
  const cell = math({ latex: "1/3" }, { digits: 6 });
  assert.equal(setting(doc, cell, "digits"), 6);
  assert.equal(setting(doc, cell, "numbers"), "both");
  assert.equal(setting({ preset: "calculator", cells: [] }, math({ latex: "1" }), "steps"), true);
});

test("a pasted text of several lines is a cell per line", () => {
  assert.deepEqual(splitLines("$$\\begin{aligned} p &\\in [0.5, 4] \\\\ x^2 - p &= 0 \\end{aligned}$$"), ["p \\in [0.5, 4]", "x^2 - p = 0"]);
  assert.deepEqual(splitLines("a = 2\n\n3a"), ["a = 2", "3a"]);
  assert.deepEqual(splitLines("\\[ \\frac{1}{2} \\]"), ["\\frac{1}{2}"]);
  // one line stays one line, a fraction's braces and all
  assert.deepEqual(splitLines("\\frac{3+\\frac{7}{12}}{\\frac{19}{84}}"), ["\\frac{3+\\frac{7}{12}}{\\frac{19}{84}}"]);
});
