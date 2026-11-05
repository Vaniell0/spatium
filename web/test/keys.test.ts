// Every function key, its call filled in, is a formula the editor's parser reads: no key inserts
// what the page would then refuse (the owner's sin, ln, lg, log that did not work).
// @ts-ignore: node's modules, present where the test runs
import { test } from "node:test";
// @ts-ignore
import assert from "node:assert/strict";
import { notation } from "../src/answers.js";
import { functionKeys } from "../src/keys.js";

// @ts-ignore: a vendored module with no types
const ls = await import(new URL("../../vendor/latex-syntax/latex-syntax.min.js", import.meta.url).href);

test("every function key's call is read, with no error in it", () => {
  for (const k of functionKeys) {
    const tree = JSON.stringify(ls.parse(notation(k.example)));
    assert.ok(!tree.includes("Error"), `${k.label}: ${tree}`);
  }
});
