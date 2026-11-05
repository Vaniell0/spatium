// The calculator's own state: the line being written, and the history of what was asked, kept in
// this browser. An entry is what its answer said in one line, so the history draws without asking
// the core; opening one asks for its whole sheet again.

import type { Written } from "./doc.js";

export interface Entry {
  written: Written;   // the question, as it was written in the line
  result: string;     // the answer's result, LaTeX
  strength: string;   // how strongly it holds: "exact", "proved to 16 digits"
  note: string;       // its first note on how, or ""
  ok: boolean;
  t: number;          // when it was asked, ms
}

export interface Calc {
  line: Written;
  history: Entry[];   // the newest first
}

const KEY = "spatium.calc";
export const MOST = 50;

export const empty = (): Written => ({ rows: [""] });

// what a question says, with its spacing set aside: two entries that say the same are one
const said = (w: Written) => ("rows" in w ? w.rows : [w.text]).map((r) => r.replace(/\s+/g, "")).filter(Boolean).join("\n");

export const blank = (w: Written) => said(w) === "";

// an entry on top of the history; one that says the same moves up instead of being kept twice
export function commit(calc: Calc, entry: Entry): Calc {
  if (blank(entry.written)) return calc;
  const rest = calc.history.filter((e) => said(e.written) !== said(entry.written));
  return { ...calc, history: [entry, ...rest].slice(0, MOST) };
}

export function load(): Calc {
  try {
    const kept = JSON.parse(localStorage.getItem(KEY) ?? "null") as Partial<Calc> | null;
    if (kept && Array.isArray(kept.history) && kept.line && typeof kept.line === "object")
      return { line: kept.line, history: kept.history.filter((e) => e && typeof e.result === "string" && e.written).slice(0, MOST) };
  } catch { /* nothing kept */ }
  return { line: empty(), history: [] };
}

export function save(calc: Calc) {
  try { localStorage.setItem(KEY, JSON.stringify(calc)); } catch { /* not kept */ }
}
