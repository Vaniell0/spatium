// A page of exercises' notation, read by the formula editor's own parser once it is written in
// the parser's: a division by a colon keeps the order of operations, a decimal comma is a point.
// @ts-ignore: node's modules, present where the test runs
import { test } from "node:test";
// @ts-ignore
import assert from "node:assert/strict";
import { notation } from "../src/answers.js";

// the editor's parser, vendored beside the page
// @ts-ignore: a vendored module with no types
const ls = await import(new URL("../../vendor/latex-syntax/latex-syntax.min.js", import.meta.url).href);

test("a : b is a quotient, tighter than a sum", () => {
  const q = String.raw`\left(1.5:0.003\right):62\frac{1}{20}+1.364:0.124`;
  const tree = ls.parse(notation(q));
  assert.equal(tree[0], "Add");
  assert.equal(tree[2][0], "Divide");
  // without it the colon is a relation looser than the sum, and the sum is divided
  assert.equal(ls.parse(q)[0], "Colon");
});

test("a decimal comma is a decimal point", () => {
  assert.deepEqual(ls.parse(notation("3{,}2+1")), ["Add", { num: "3.2" }, 1]);
  assert.deepEqual(ls.parse("3{,}2+1"), ["Add", 32, 1]);
});

test("a command's own colon is left alone", () => {
  assert.equal(notation("a\\:b"), "a\\:b");
});
