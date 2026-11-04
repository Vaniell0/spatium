// The cells of a document, one under another: on the left the name a cell defines, then what is
// written in it -- a row of the formula editor, a line of the language, or the graph of a name --
// with the settings its policy offers, and under it the answer. What a person may do and what
// the page does on its own come from the document's policy (doc.ts), never from here.

import type { Answer } from "../core.js";
import { type Cell, type Doc, type Settings, graph, math, presets, setting, splitLines } from "../doc.js";
import { nameOf } from "../answers.js";
import { answerView } from "./answer.js";
import { el } from "./dom.js";
import { graphView, type View } from "./graph.js";

export interface Editor {
  ready: boolean;
  toLatex(text: string): string;    // a line of the language as a formula (MathLive's ascii-math reader)
  toText(latex: string): string;    // a formula as a line, when the core has not read it yet
}

export interface SheetHost {
  editor: Editor;
  source: () => boolean;             // formulas as LaTeX source
  changed(): void;                   // the document changed: ask again and keep it
  focused(field: HTMLElement): void; // a formula row took the focus (the keyboard goes there)
  completion?: (input: HTMLTextAreaElement) => void;
}

interface Shown { cell: Cell; box: HTMLElement; answer: Answer | null; name: string | null }

export class Sheet {
  readonly root = el("div", { class: "sheet" });
  private shown = new Map<Cell, Shown>();

  constructor(private doc: Doc, private readonly host: SheetHost) {
    this.build();
  }

  setDoc(doc: Doc) {
    this.doc = doc;
    this.build();
  }

  get document() { return this.doc; }

  // every cell drawn again from the document; answers come as they are asked
  build(focus?: Cell) {
    this.shown.clear();
    this.root.replaceChildren();
    if (!this.doc.cells.length) this.doc.cells.push(math({ latex: "" }));
    for (const cell of this.doc.cells) this.root.append(this.cellBox(cell));
    this.root.append(this.adder());
    if (focus) this.focus(focus);
  }

  // an answer arrived for a cell, or for a derived cell after the cell `from`
  show(cell: Cell, answer: Answer | null, from?: Cell) {
    let s = this.shown.get(cell);
    if (!s && from) {
      const parent = this.shown.get(from);
      if (!parent) return;
      const box = this.cellBox(cell, true);
      parent.box.after(box);
      s = this.shown.get(cell)!;
    }
    if (!s) return;
    s.answer = answer;
    s.name = cell.kind === "math" ? nameOf(answer ?? undefined) : null;
    const name = s.box.querySelector(".name")!;
    name.textContent = cell.kind === "graph" ? "graph" : s.name ? `${s.name} =` : "";
    const out = s.box.querySelector<HTMLElement>(".answer")!;
    if (!answer) { out.replaceChildren(); return; }
    if (cell.kind === "graph") {
      out.replaceChildren(graphView(answer, {
        view: cell.settings.window as View | undefined,
        height: 260,
        onView: (v) => { if (!this.isDerived(cell)) { cell.settings.window = v; this.host.changed(); } },
      }));
      return;
    }
    const typed = "text" in cell.written ? cell.written.text : "";
    out.replaceChildren(answerView(answer, typed, {
      digits: setting(this.doc, cell, "digits"),
      steps: setting(this.doc, cell, "steps"),
      source: this.host.source(),
    }));
    s.box.classList.toggle("refused", Boolean(answer.refused));
  }

  // derived cells the run no longer gave are dropped before a new run
  clearDerived() {
    for (const [cell, s] of this.shown) if (this.isDerived(cell)) { s.box.remove(); this.shown.delete(cell); }
  }

  private isDerived(cell: Cell) { return !this.doc.cells.includes(cell); }

  private cellBox(cell: Cell, derived = false): HTMLElement {
    const policy = presets[this.doc.preset];
    const box = el("section", { class: derived ? "cell derived" : "cell", "data-kind": cell.kind });
    const bar = el("div", { class: "cellbar" });
    const input = el("div", { class: "input" });
    bar.append(input);
    if (!derived) {
      if (cell.kind === "math") {
        bar.append(this.asSwitch(cell));
        for (const key of policy.settings) bar.append(this.settingControl(cell, key));
      }
      if (this.doc.cells.length > 1)
        bar.append(el("button", { type: "button", class: "icon small", "aria-label": "Remove the cell", text: "×",
          onclick: () => this.remove(cell) }));
    }
    this.fillInput(cell, input, derived);
    box.append(el("div", { class: "gutter" }, el("span", { class: "name mono" })),
      el("div", { class: "body" }, bar, el("div", { class: "answer", "aria-live": "polite" })));
    this.shown.set(cell, { cell, box, answer: null, name: null });
    return box;
  }

  private fillInput(cell: Cell, input: HTMLElement, derived: boolean) {
    if (cell.kind === "graph") {
      input.append(el("span", { class: "mono of", text: derived ? `graph of ${cell.of}` : `graph(${cell.of})` }));
      return;
    }
    if ("latex" in cell.written && this.host.editor.ready) {
      const mf = document.createElement("math-field") as HTMLElement & { value: string; mathVirtualKeyboardPolicy: string };
      mf.mathVirtualKeyboardPolicy = "manual";
      mf.setAttribute("aria-label", "Formula");
      mf.value = cell.written.latex;
      mf.addEventListener("input", () => { (cell.written as { latex: string }).latex = mf.value; this.host.changed(); });
      mf.addEventListener("focusin", () => this.host.focused(mf));
      mf.addEventListener("keydown", (e) => this.rowKeys(e, cell, mf.value), { capture: true });
      mf.addEventListener("paste", (e) => this.pasted(e, cell), { capture: true });
      input.append(mf);
      return;
    }
    const written = cell.written as { text?: string; latex?: string };
    const t = el("textarea", { class: "text", rows: "1", spellcheck: "false", autocapitalize: "off", "aria-label": "Line in the language",
      placeholder: "2x^2 - 3x - 5 = 0" });
    t.value = written.text ?? written.latex ?? "";
    t.addEventListener("input", () => {
      cell.written = { text: t.value };
      this.host.completion?.(t);
      this.host.changed();
    });
    t.addEventListener("keydown", (e) => this.rowKeys(e, cell, t.value));
    t.addEventListener("paste", (e) => this.pasted(e, cell));
    input.append(t);
  }

  // Enter opens a cell below; Backspace in an empty cell removes it
  private rowKeys(e: KeyboardEvent, cell: Cell, value: string) {
    if (e.defaultPrevented) return;
    if (e.key === "Enter" && !e.shiftKey) {
      e.preventDefault();
      this.insertAfter(cell, [math(cell.kind === "math" && "text" in cell.written ? { text: "" } : { latex: "" })]);
    } else if (e.key === "Backspace" && value === "" && this.doc.cells.length > 1) {
      e.preventDefault();
      this.remove(cell);
    }
  }

  // a pasted text of several lines: a cell per line (the policy's "split")
  private pasted(e: ClipboardEvent, cell: Cell) {
    if (!presets[this.doc.preset].auto.includes("split")) return;
    const lines = splitLines(e.clipboardData?.getData("text/plain") ?? "");
    if (lines.length < 2 || cell.kind !== "math") return;
    e.preventDefault();
    e.stopPropagation();
    const asLatex = "latex" in cell.written;
    const cells = lines.map((l) => math(asLatex ? { latex: l } : { text: l }));
    const at = this.doc.cells.indexOf(cell);
    const empty = asLatex ? !(cell.written as { latex: string }).latex.trim() : !(cell.written as { text: string }).text.trim();
    this.doc.cells.splice(empty ? at : at + 1, empty ? 1 : 0, ...cells);
    this.build(cells[cells.length - 1]);
    this.host.changed();
  }

  // LaTeX and the language are two ways to write one cell: switching carries what is written over
  private asSwitch(cell: Cell & { kind: "math" }) {
    const isLatex = "latex" in cell.written;
    const b = el("button", { type: "button", class: "quiet small", text: isLatex ? "TeX" : "Text",
      title: isLatex ? "Written as a formula; switch to a line of the language" : "Written in the language; switch to a formula" });
    b.addEventListener("click", () => {
      if ("latex" in cell.written) {
        const read = this.shown.get(cell)?.answer?.question;
        const line = typeof read === "string" && read.trim() ? read.split("\n").at(-1)! : this.host.editor.toText(cell.written.latex);
        cell.written = { text: line };
      } else {
        cell.written = { latex: this.host.editor.toLatex(cell.written.text) };
      }
      this.build(cell);
      this.host.changed();
    });
    if (!this.host.editor.ready && !isLatex) b.disabled = true;
    return b;
  }

  private settingControl(cell: Cell, key: keyof Settings): HTMLElement {
    if (key === "digits") {
      const now = setting(this.doc, cell, "digits");
      const next = now === 6 ? 12 : now === 12 ? 17 : 6;
      return el("button", { type: "button", class: "quiet small", text: `${now} digits`, title: "Significant digits shown",
        onclick: () => { cell.settings.digits = next; this.rebuildOne(cell); this.host.changed(); } });
    }
    if (key === "steps") {
      const on = setting(this.doc, cell, "steps");
      return el("button", { type: "button", class: "quiet small", "aria-pressed": String(on), text: "steps",
        onclick: () => { cell.settings.steps = !on; this.rebuildOne(cell); this.host.changed(); } });
    }
    return el("span");
  }

  private rebuildOne(cell: Cell) {
    const s = this.shown.get(cell);
    if (!s) return;
    const box = this.cellBox(cell);
    s.box.replaceWith(box);
    this.show(cell, s.answer);
  }

  private adder(): HTMLElement {
    const policy = presets[this.doc.preset];
    const row = el("div", { class: "adder" });
    const add = (label: string, make: () => Cell) =>
      row.append(el("button", { type: "button", class: "quiet", text: label, onclick: () => {
        const cell = make();
        this.doc.cells.push(cell);
        this.build(cell);
        this.host.changed();
      } }));
    add("+ Formula", () => math({ latex: "" }));
    if (policy.may.add === "any") {
      add("+ Text", () => math({ text: "" }));
      const names = [...this.shown.values()].map((s) => s.name).filter((n): n is string => Boolean(n));
      const last = names.at(-1);
      if (last) add(`+ Graph of ${last}`, () => graph(last));
    }
    return row;
  }

  insertAfter(cell: Cell, cells: Cell[]) {
    const at = this.doc.cells.indexOf(cell);
    this.doc.cells.splice(at + 1, 0, ...cells);
    this.build(cells[0]);
    this.host.changed();
  }

  private remove(cell: Cell) {
    const at = this.doc.cells.indexOf(cell);
    if (at < 0) return;
    this.doc.cells.splice(at, 1);
    this.build(this.doc.cells[Math.max(0, at - 1)]);
    this.host.changed();
  }

  focus(cell: Cell) {
    const s = this.shown.get(cell);
    const f = s?.box.querySelector<HTMLElement>("math-field, textarea");
    if (f) requestAnimationFrame(() => f.focus());
  }
}
