// The canvas as a graph of nodes: each cell of the document is a node, wired to the nodes whose
// names its text reads, and placed in layers left to right by dependency. Nothing here touches
// the page; the view (view/canvas.ts) draws what these functions say, and a test holds them.

import type { Cell, Doc } from "./doc.js";

export type Kind = "formula" | "graph";

export interface Node {
  id: string;          // "formula-1", "graph-2": a label, not a name of the language
  kind: Kind;
  cell: Cell;
  text: string;        // the source as the person edits it: the language, or a formula row's LaTeX
  defines: string[];   // the names its line defines
  reads: string[];     // the names it reads that other nodes define
  w: number;
  at?: [number, number];   // pinned by hand; placed automatically otherwise
}

export const WIDTH: Record<Kind, number> = { formula: 232, graph: 360 };

// the names a line of the language defines: `r = 0.5`, `p in [a, b]`, `x: f = g`, `x, y: ...`,
// `I = integral ...`; a shorthand line with no name (`2 + 2`, `x^2 = 2`) defines none here
export function definedBy(text: string): string[] {
  const line = text.split("\n").find((l) => l.trim() && !l.trim().startsWith("#")) ?? "";
  const m = /^\s*([A-Za-z_]\w*(?:\s*,\s*[A-Za-z_]\w*)*)\s*(?::|=(?!=)|\bin\b)/.exec(line);
  return m ? m[1].split(",").map((s) => s.trim()) : [];
}

// the identifiers a text mentions, apart from the names it defines itself
export function mentioned(text: string): string[] {
  const own = new Set(definedBy(text));
  const seen = new Set<string>();
  for (const m of text.replace(/#.*$/gm, "").matchAll(/[A-Za-z_]\w*/g)) if (!own.has(m[0])) seen.add(m[0]);
  return [...seen];
}

// the document's cells as nodes, wired by name in the order written: a node reads the nearest
// earlier definition of a name
export function nodesOf(doc: Doc): Node[] {
  const count: Record<Kind, number> = { formula: 0, graph: 0 };
  const defined = new Map<string, string>();   // name -> the id that defines it
  return doc.cells.map((cell) => {
    const kind: Kind = cell.kind === "graph" ? "graph" : "formula";
    const id = `${kind}-${++count[kind]}`;
    const text = cell.kind === "graph" ? "" : "text" in cell.written ? cell.written.text : cell.written.rows.join("\n");
    const defines = cell.kind === "graph" ? [] : definedBy(text);
    const wants = cell.kind === "graph" ? [cell.of] : mentioned(text);
    const reads = [...new Set(wants.filter((n) => defined.has(n)).map((n) => defined.get(n)!))];
    for (const n of defines) defined.set(n, id);
    const at = (cell.settings as { at?: [number, number] }).at;
    return { id, kind, cell, text, defines, reads, w: WIDTH[kind], at };
  });
}

export interface Placed { id: string; x: number; y: number; w: number; h: number; column: number }

const snap = (v: number) => Math.round(v / 8) * 8;
const GAP_X = 64, GAP_Y = 32, LEFT = 48, TOP = 48, CLEAR = 16;

// Layered, left to right by dependency: a node's column is the length of the longest chain of
// wires into it; within a column nodes follow the mean y of their parents (the barycentre), 32
// apart (16 sideways); a column's x is the left edge plus the widths of the columns before it, 64 apart; every
// coordinate snaps to 8. A pinned node keeps its place, and the rest are laid out around it with
// 16 px clearance. The same input gives the same output.
export function layout(nodes: Node[], heights: (id: string) => number): Placed[] {
  const by = new Map(nodes.map((n) => [n.id, n]));
  const column = new Map<string, number>();
  const depth = (id: string, path: Set<string>): number => {
    if (column.has(id)) return column.get(id)!;
    if (path.has(id)) return 0;   // a cycle: told elsewhere, not laid out twice
    path.add(id);
    const d = Math.max(-1, ...(by.get(id)?.reads ?? []).map((r) => depth(r, path))) + 1;
    column.set(id, d);
    return d;
  };
  for (const n of nodes) depth(n.id, new Set());
  const columns = Math.max(0, ...column.values()) + 1;
  const widths = Array.from({ length: columns }, (_, c) => Math.max(0, ...nodes.filter((n) => column.get(n.id) === c).map((n) => n.w)));
  const xs = widths.map((_, c) => LEFT + widths.slice(0, c).reduce((s, w) => s + w + GAP_X, 0));
  const placed: Placed[] = [];
  const free = (x: number, y: number, w: number, h: number) =>
    !placed.some((p) => x < p.x + p.w + CLEAR && x + w + CLEAR > p.x && y < p.y + p.h + GAP_Y && y + h + GAP_Y > p.y);
  for (const n of nodes) if (n.at) placed.push({ id: n.id, x: n.at[0], y: n.at[1], w: n.w, h: heights(n.id), column: column.get(n.id)! });
  for (let c = 0; c < columns; ++c) {
    const here = nodes.filter((n) => column.get(n.id) === c && !n.at);
    const wish = (n: Node) => {
      const parents = n.reads.map((r) => placed.find((p) => p.id === r)).filter((p): p is Placed => Boolean(p));
      return parents.length ? parents.reduce((s, p) => s + p.y + p.h / 2, 0) / parents.length : TOP;
    };
    here.sort((a, b) => wish(a) - wish(b) || nodes.indexOf(a) - nodes.indexOf(b));
    for (const n of here) {
      const h = heights(n.id), x = xs[c];
      const want = snap(Math.max(TOP, wish(n) - h / 2));
      let y = want;
      // the first free slot nearest the wish: down, then up, in 8 px steps
      for (let k = 0; !free(x, y, n.w, h) && k < 4000; ++k) y = k % 2 ? want - 8 * Math.ceil(k / 2) : want + 8 * Math.ceil(k / 2);
      if (y < TOP) { y = want; while (!free(x, y, n.w, h)) y += 8; }
      placed.push({ id: n.id, x, y, w: n.w, h, column: c });
    }
  }
  return placed;
}

// wires by name, as the view draws them: from the defining node to the reader
export function wires(nodes: Node[]): [from: string, to: string][] {
  return nodes.flatMap((n) => n.reads.map((r) => [r, n.id] as [string, string]));
}

// a cubic Bezier from (x1, y1) to (x2, y2), its control points d to the right and left,
// d = max(64, |x2 - x1| / 2)
export function wirePath(x1: number, y1: number, x2: number, y2: number): string {
  const d = Math.max(64, Math.abs(x2 - x1) / 2);
  return `M${x1} ${y1}C${x1 + d} ${y1} ${x2 - d} ${y2} ${x2} ${y2}`;
}
