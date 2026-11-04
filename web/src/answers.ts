// What the core answers for each cell. A cell is read in the document the cells above it make,
// so its question is every row up to it: a LaTeX cell as the MathJSON its parser writes, a cell
// in the language as its text. A graph cell draws what the cell it names answers with
// (`graph`); its window is the view's, a uniform of the shader, never a new question. Answers are kept by question, so a cell
// nothing above it changed is not asked again.

import type { Answer, Core } from "./core.js";
import { type Cell, type Doc, graph, presets } from "./doc.js";

export type Parse = (latex: string) => unknown;

// A formula as a page of exercises writes it, in the notation the editor's parser reads: a decimal
// comma (3{,}2, as the editor writes it under a comma locale) as a point, and a division written
// a : b as a \div b, which the parser would otherwise read as a relation looser than a sum.
export function notation(latex: string): string {
  return latex.replace(/\{,\}/g, ".").replace(/(?<!\\):/g, "\\div ");
}

// the row the core reads a math cell as, or null for an empty one
export function rowOf(cell: Cell, parse: Parse): unknown | null {
  if (cell.kind !== "math") return null;
  if ("text" in cell.written) return cell.written.text.trim() ? { text: cell.written.text } : null;
  return cell.written.latex.trim() ? parse(cell.written.latex) : null;
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
  // null for a cell that has nothing to ask, and `from` the cell a derived one follows.
  async run(doc: Doc, shown: (cell: Cell, answer: Answer | null, from?: Cell) => void): Promise<void> {
    const rows: unknown[] = [];
    const answers = new Map<string, Answer>();   // name -> the answer that defines it
    const autoGraph = presets[doc.preset].auto.includes("graph");
    const draw = async (cell: Cell & { kind: "graph" }, from?: Cell) => {
      const of = answers.get(cell.of);
      const request = of?.graph as { document: string } | undefined;
      if (!request) { shown(cell, of ? { refused: `${cell.of} has no graph` } : null, from); return true; }
      const drawn = await this.asked(cell, "plot", JSON.stringify(request));
      if (drawn === null) return false;   // a newer run took this one's place
      shown(cell, drawn, from);
      return true;
    };
    for (const cell of doc.cells) {
      if (cell.kind === "graph") {
        if (!(await draw(cell))) return;
        continue;
      }
      const row = rowOf(cell, this.parse);
      if (row === null) { shown(cell, null); continue; }
      rows.push(row);
      const answer = await this.asked(cell, "rows", JSON.stringify({ rows }));
      if (answer === null) return;   // a newer run took this one's place
      const name = nameOf(answer);
      if (name) answers.set(name, answer);
      shown(cell, answer);
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
