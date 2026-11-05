// Where the page shows each word of the library's language. The words are the library's own
// (gen/language.ts, written by examples/site_reference.cpp from the C++ tables); each list here is
// a record over its words, so a word the library gains is a type error in this file until it is
// given a place on the page, or named a gap with what would show it. The coverage test prints the
// gaps, and the list of them only shrinks.

import * as L from "./gen/language.js";

export type Home = { at: string } | { gap: string };

// a list the page draws whole from site/data/reference.json, so every word of it is there by
// construction
function every<K extends string>(words: readonly K[], at: string): Record<K, Home> {
  return Object.fromEntries(words.map((w) => [w, { at }])) as Record<K, Home>;
}

const sheet = { at: "calculator: the answer's sheet" };
const axis = { at: "calculator: the inequality's axis" };
const noShape = { gap: "no page draws a shape; the canvas scene node will" };
const noControl = { gap: "no control computes in it; the canvas Inspector's Numbers will" };

export const homes = {
  functions: every(L.functions, "reference: Functions, and the hint under the caret"),
  operators: every(L.operators, "reference: Operators"),
  constants: every(L.constants, "reference: Constants, and the hint under the caret"),
  statements: {
    value: sheet,
    parameter: { at: "calculator: a parameter row, and the plot's slider" },
    root: sheet,
    integral: sheet,
    shape: noShape,
  } satisfies Record<L.StatementKind, Home>,
  relations: { "=": sheet, "<": axis, "<=": axis, ">": axis, ">=": axis } satisfies Record<L.Relation, Home>,
  routes: {
    "evaluate": sheet,
    "expression": sheet,
    "closed form": sheet,
    "linear system": sheet,
    "interval method": axis,
    "every root, proved": sheet,
    "root on a bracket": sheet,
    "integral": sheet,
    "synthesis": sheet,
  } satisfies Record<L.Route, Home>,
  arithmetics: {
    double: { at: "canvas: a card's Numbers, decimal; the calculator computes in it" },
    float: noControl,
    rational: { at: "canvas: a card's Numbers, exact; the calculator where every operation is one on fractions" },
    decimal: noControl,
    real50: noControl,
  } satisfies Record<L.Arithmetic, Home>,
  shapes: { torus: noShape, cylinder: noShape, sphere: noShape, flake: noShape, offset: noShape, scatter: noShape, union: noShape } satisfies
    Record<L.ShapeOperation, Home>,
};

// every word with no place yet, as "list: word: what would show it"
export function gaps(): string[] {
  return Object.entries(homes).flatMap(([list, words]) =>
    Object.entries(words as Record<string, Home>).flatMap(([w, h]) => ("gap" in h ? [`${list}: ${w}: ${h.gap}`] : [])));
}
