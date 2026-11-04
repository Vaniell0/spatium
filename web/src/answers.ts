// What the core answers for each cell. A cell's question is its own rows, or, where the policy
// reads above, every row up to it: a row of a formula as the MathJSON the editor's parser
// writes, a text in the language as it is. A graph cell draws what the cell it names answers
// with (`graph`); its window is the view's, a uniform of the shader, never a new question.
// Answers are kept by question, so a cell whose question did not change is not asked again.

import type { Answer, Core } from "./core.js";
import { type Cell, type Doc, graph, presets, setting } from "./doc.js";

export type Parse = (latex: string) => unknown;

// A formula as a page of exercises writes it, in the notation the editor's parser reads: a decimal
// comma (3{,}2, as the editor writes it under a comma locale) as a point, and a division written
// a : b as a \div b, which the parser would otherwise read as a relation looser than a sum.
export function notation(latex: string): string {
  return latex.replace(/\{,\}/g, ".").replace(/(?<!\\):/g, "\\div ");
}

// the rows the core reads a math cell as, its empty rows left out
export function rowsOf(cell: Cell, parse: Parse): unknown[] {
  if (cell.kind !== "math") return [];
  if ("text" in cell.written) return cell.written.text.trim() ? [{ text: cell.written.text }] : [];
  return cell.written.rows.filter((r) => r.trim()).map(parse);
}

// The word a caret stands after, for a hint from the reference: a function the editor set
// (\sin, \operatorname{erf}), or letters typed one by one; a single letter is a name, not a word.
export function wordBefore(text: string): string {
  const t = text.replace(/(\\left)?\(\s*$|\s+$/, "");
  return (/\\(?:operatorname\{|mathrm\{)?([A-Za-z]+)\}?$/.exec(t) ?? /(?:^|[^A-Za-z\\])([A-Za-z]{2,})$/.exec(t))?.[1] ?? "";
}

// The name a math cell's answer defines: the last statement of the document it was read as.
export function nameOf(answer: Answer | undefined): string | null {
  const doc = answer?.document as { statements?: { name?: string }[] } | undefined;
  const last = doc?.statements?.at(-1);
  return typeof last?.name === "string" ? last.name : null;
}

export class Answers {
  private kept = new Map<string, Answer>();   // question -> answer

  constructor(private readonly core: Core, private readonly parse: Parse) {}

  // A cell the page adds on its own (the policy's "graph") follows the cell it came from; it is
  // the same object from run to run, so a view and the core's queue know it again.
  private derived = new WeakMap<Cell, Cell>();

  // The answer of every cell, in order, as each arrives: `shown` is called once per cell, with
  // null for a cell that has nothing to ask, `above` the rows asked before the cell's own (a
  // refused `row` counts them), and `from` the cell a derived one follows.
  async run(doc: Doc, shown: (cell: Cell, answer: Answer | null, above: number, from?: Cell) => void): Promise<void> {
    const rows: unknown[] = [];
    const policy = presets[doc.preset];
    const answers = new Map<string, Answer>();   // name -> the answer that defines it
    const autoGraph = policy.auto.includes("graph");
    const draw = async (cell: Cell & { kind: "graph" }, from?: Cell) => {
      const of = answers.get(cell.of);
      const request = of?.graph as { document: string } | undefined;
      if (!request) { shown(cell, of ? { refused: `${cell.of} has no graph` } : null, 0, from); return true; }
      const drawn = await this.asked(cell, "plot", JSON.stringify(request));
      if (drawn === null) return false;   // a newer run took this one's place
      shown(cell, drawn, 0, from);
      return true;
    };
    for (const cell of doc.cells) {
      if (cell.kind === "graph") {
        if (!(await draw(cell))) return;
        continue;
      }
      const own = rowsOf(cell, this.parse);
      if (policy.reads === "own") rows.length = 0;
      const above = rows.length;
      if (!own.length) { shown(cell, null, above); continue; }
      rows.push(...own);
      const exact = setting(doc, cell, "numbers") === "exact";
      const answer = await this.asked(cell, "rows", JSON.stringify(exact ? { rows, engine: { arithmetic: "rational" } } : { rows }));
      if (answer === null) return;   // a newer run took this one's place
      const name = nameOf(answer);
      if (name) answers.set(name, answer);
      shown(cell, answer, above);
      if (autoGraph && name && answer.graph) {
        let g = this.derived.get(cell);
        if (!g || (g.kind === "graph" && g.of !== name)) this.derived.set(cell, (g = graph(name)));
        if (!(await draw(g as Cell & { kind: "graph" }, cell))) return;
      }
    }
  }

  private async asked(cell: Cell, kind: "rows" | "plot", payload: string): Promise<Answer | null> {
    const key = `${kind} ${payload}`;
    const hit = this.kept.get(key);
    if (hit) return hit;
    const answer = await this.core.ask(cell, kind, payload);
    if (answer) this.kept.set(key, answer);
    return answer;
  }
}
