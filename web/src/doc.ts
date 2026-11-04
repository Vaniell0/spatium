// The document a page holds: cells, each a line written in LaTeX or in the language, or a view
// of another cell read by its name. The calculator and the canvas are two presets of one
// document under two policies (below): what a person may do, and what the page does on its own.

export type Written = { latex: string } | { text: string };

export type Numbers = "exact" | "decimal" | "both";
export type Digits = 6 | 12 | 17;

// How a cell is computed and shown; a setting left out is the preset's.
export interface Settings {
  numbers?: Numbers;
  digits?: Digits;
  steps?: boolean;            // the chain of steps under an answer
  method?: string;            // a route the core names; none is the core's own choice
  window?: [number, number, number, number];  // a graph's view, x0 y0 x1 y1, as panned and zoomed
}

export type Cell =
  | { kind: "math"; written: Written; settings: Settings }
  | { kind: "graph"; of: string; settings: Settings };

export type Preset = "calculator" | "canvas";

export interface Doc {
  preset: Preset;
  cells: Cell[];
}

// A preset is a policy: what a person may do with the cells, and what the page does to what is
// typed on its own. The two pages differ only in their policies.
export interface Policy {
  settings: (keyof Settings)[];   // the settings a cell offers
  defaults: Required<Pick<Settings, "numbers" | "digits" | "steps">>;
  may: {
    add: "input" | "any";         // a calculator adds lines to type in; the canvas, any kind of cell
    place: boolean;               // move and size cells freely
    reorder: boolean;             // move a cell above the cells it reads is refused either way
  };
  auto: Auto[];                   // done to the cells as they are typed, in this order
}

// What a page does on its own. Cells it adds are derived from the answers, not kept in the
// document: they come back from the same input, and a link carries only what was typed.
//   split    a pasted text of several lines (LaTeX \\ or newlines) becomes a cell per line
//   graph    an answer that names a graph is followed by that graph
//   sliders  a name an answer reads and nothing defines becomes a parameter to slide
export type Auto = "split" | "graph" | "sliders";

export const presets: Record<Preset, Policy> = {
  calculator: {
    settings: [],
    defaults: { numbers: "both", digits: 12, steps: true },
    may: { add: "input", place: false, reorder: false },
    auto: ["split", "graph", "sliders"],
  },
  canvas: {
    settings: ["digits", "steps"],
    defaults: { numbers: "both", digits: 12, steps: true },
    may: { add: "any", place: true, reorder: true },
    auto: ["split"],
  },
};

export function setting<K extends keyof Policy["defaults"]>(doc: Doc, cell: Cell, key: K): Policy["defaults"][K] {
  return (cell.settings[key] as Policy["defaults"][K] | undefined) ?? presets[doc.preset].defaults[key];
}

export const math = (written: Written, settings: Settings = {}): Cell => ({ kind: "math", written, settings });
export const graph = (of: string, settings: Settings = {}): Cell => ({ kind: "graph", of, settings });

// ── the document as a link: #w= ─────────────────────────────────────────
// [2, preset, cells...], a cell an array led by a letter: ["l", latex], ["t", text], ["g", of],
// with its settings as a third element when it has any. JSON, deflate-raw, base64url.
// A link of the windows page before it ([1, windows...]) is read as the cells it held.

type Packed = [string, string, Settings?];

function pack(cell: Cell): Packed {
  const s = Object.keys(cell.settings).length ? cell.settings : undefined;
  const head: [string, string] =
    cell.kind === "graph" ? ["g", cell.of] : "latex" in cell.written ? ["l", cell.written.latex] : ["t", cell.written.text];
  return s ? [...head, s] : head;
}

const DIGITS: readonly number[] = [6, 12, 17];
const NUMBERS: readonly string[] = ["exact", "decimal", "both"];

// only the settings a cell can have, with the values they can take: a link is read, not trusted
function settingsOf(x: unknown): Settings {
  if (!x || typeof x !== "object") return {};
  const o = x as Record<string, unknown>, s: Settings = {};
  if (typeof o.numbers === "string" && NUMBERS.includes(o.numbers)) s.numbers = o.numbers as Numbers;
  if (typeof o.digits === "number" && DIGITS.includes(o.digits)) s.digits = o.digits as Digits;
  if (typeof o.steps === "boolean") s.steps = o.steps;
  if (typeof o.method === "string") s.method = o.method;
  const w = o.window;
  if (Array.isArray(w) && w.length === 4 && w.every((v) => typeof v === "number" && isFinite(v)) && w[0] < w[2] && w[1] < w[3])
    s.window = [w[0], w[1], w[2], w[3]];
  return s;
}

function unpack(x: unknown): Cell | null {
  if (!Array.isArray(x) || typeof x[0] !== "string" || typeof x[1] !== "string") return null;
  const settings = settingsOf(x[2]);
  if (x[0] === "l") return math({ latex: x[1] }, settings);
  if (x[0] === "t") return math({ text: x[1] }, settings);
  if (x[0] === "g") return graph(x[1], settings);
  return null;
}

export function toJSON(doc: Doc): unknown[] {
  return [2, doc.preset, ...doc.cells.map(pack)];
}

// The document a JSON value names, or null when it names none. Cells it cannot read are left out.
export function fromJSON(x: unknown): Doc | null {
  if (!Array.isArray(x)) return null;
  if (x[0] === 2 && (x[1] === "calculator" || x[1] === "canvas"))
    return { preset: x[1], cells: x.slice(2).map(unpack).filter((c): c is Cell => c !== null) };
  if (x[0] === 1) {
    // the windows page: each row of a formula window a cell, a text window one; the reference
    // and the gallery are panels of the page now, not cells
    const cells: Cell[] = [];
    for (const w of x.slice(1)) {
      if (!Array.isArray(w) || typeof w[0] !== "string") continue;
      const letter = w[0].replace(/\+$/, "");
      if (letter === "f") for (const r of w.slice(1)) if (typeof r === "string" && r.trim()) cells.push(math({ latex: r }));
      if (letter === "t" && typeof w[1] === "string") cells.push(math({ text: w[1] }));
    }
    return { preset: "canvas", cells };
  }
  return null;
}

const base64url = (bytes: Uint8Array) =>
  btoa(String.fromCharCode(...bytes)).replace(/\+/g, "-").replace(/\//g, "_").replace(/=+$/, "");
const unbase64url = (s: string) => Uint8Array.from(atob(s.replace(/-/g, "+").replace(/_/g, "/")), (c) => c.charCodeAt(0));
const through = (bytes: Uint8Array<ArrayBuffer> | string, stream: CompressionStream | DecompressionStream) =>
  new Response(new Blob([bytes]).stream().pipeThrough(stream)).arrayBuffer();

export async function toLink(doc: Doc): Promise<string> {
  return base64url(new Uint8Array(await through(JSON.stringify(toJSON(doc)), new CompressionStream("deflate-raw"))));
}

export async function fromLink(code: string): Promise<Doc | null> {
  try {
    return fromJSON(JSON.parse(new TextDecoder().decode(await through(unbase64url(code), new DecompressionStream("deflate-raw")))));
  } catch {
    return null;
  }
}

// ── what the page does on its own ────────────────────────────────────────

const ENVIRONMENT = /\\(begin|end)\{(aligned|align\*?|gathered|gather\*?|equation\*?|array|split)\}(\{[^}]*\})?/g;

// A pasted text as the lines it holds: display delimiters and the environments that only lay
// lines out are dropped, `\\` and newlines end a line, `&` (an alignment point) is dropped.
export function splitLines(pasted: string): string[] {
  const body = pasted.trim().replace(/^\$\$?|\$\$?$/g, "").replace(/^\\\[|\\\]$/g, "").replace(ENVIRONMENT, "");
  return body.split(/\\\\|\n/).map((line) => line.replace(/&/g, "").trim()).filter((line) => line.length > 0);
}
