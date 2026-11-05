// The page has a place for every word of the library's language, or names it a gap: the words
// are read from site/data/reference.json, which the library writes beside gen/language.ts, so a
// word the TypeScript lists lost, or one the registry has and the library no longer does, fails
// here; the gaps are printed, and their list only shrinks.
// @ts-ignore: node's modules, present where the test runs
import { test } from "node:test";
// @ts-ignore
import assert from "node:assert/strict";
// @ts-ignore
import { readFileSync } from "node:fs";
import { gaps, homes } from "../src/registry.js";

const reference = JSON.parse(readFileSync(new URL("../../data/reference.json", import.meta.url), "utf8"));
const key: Record<string, string> = { operators: "written", relations: "written" };

test("every word of the library has a home on the page, and no other word does", () => {
  for (const [list, words] of Object.entries(homes)) {
    const library = (reference[list] as Record<string, string>[]).map((w) => w[key[list] ?? "name"]);
    assert.deepEqual(Object.keys(words).sort(), [...library].sort(), list);
  }
});

// today's gaps; a page change that shows one removes it here, and a new word with no place fails
const accepted = [
  "statements: shape",
  "arithmetics: float", "arithmetics: decimal", "arithmetics: real50",
  "shapes: torus", "shapes: cylinder", "shapes: sphere", "shapes: flake", "shapes: offset", "shapes: scatter", "shapes: union",
];

test("the gaps are the accepted ones", () => {
  const now = gaps();
  console.log(now.join("\n"));
  assert.deepEqual(now.map((g) => g.split(": ").slice(0, 2).join(": ")), accepted);
});
