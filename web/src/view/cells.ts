// The cells of a document as cards, one under another. A card is one question: the rows of a
// formula -- a system, or names and what is asked of them -- or a text in the language, with
// its answer under it and the graph the answer names. Its gear opens the settings its policy
// offers. What a person may do and what the page does on its own come from the document's
// policy (doc.ts), never from here.

import type { Answer } from "../core.js";
import { type Cell, type Doc, graph, math, presets, setting, splitLines } from "../doc.js";
import { nameOf, wordBefore } from "../answers.js";
import { answerView, foldable } from "./answer.js";
import { el } from "./dom.js";
import { graphView, type GraphOptions, type View } from "./graph.js";

export interface CardsHost {
  editor: boolean;                    // the formula editor loaded; without it a card is written as text
  source: () => boolean;              // formulas as LaTeX source
  changed(): void;                    // the document changed: ask again and keep it
  focused(field: HTMLElement): void;  // a row took the focus (the keyboard goes there)
  hint(word: string): [name: string, caption: string, example: string] | null;  // what the reference says of a word
  curve: (cell: Cell) => GraphOptions["curve"];                                 // a graph's curve in doubles
  ask: (text: string) => void;                                                   // a question the sheet offers, as a new card
  commit?: () => void;                                                           // the calculator's line is done: into the history
}

type MathField = HTMLElement & { value: string; position: number; mathVirtualKeyboardPolicy: string;
  getValue(start: number, end: number, format: string): string };

interface Shown { box: HTMLElement; answer: Answer | null; above: number; name: string | null }

// the choices of a setting in a card's gear, in its words
const CHOICES = {
  numbers: { title: "Numbers", values: [["decimal", "decimal"], ["exact", "exact"]] },
  steps: { title: "Steps", values: [[true, "shown"], [false, "hidden"]] },
} as const;

const GEAR = '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M12 15.5a3.5 3.5 0 1 0 0-7 3.5 3.5 0 0 0 0 7Z"/><path d="M19.4 13.5a7.6 7.6 0 0 0 0-3l2-1.6-2-3.4-2.4 1a7.4 7.4 0 0 0-2.6-1.5L14 2.5h-4l-.4 2.5A7.4 7.4 0 0 0 7 6.5l-2.4-1-2 3.4 2 1.6a7.6 7.6 0 0 0 0 3l-2 1.6 2 3.4 2.4-1a7.4 7.4 0 0 0 2.6-1.5l.4 2.5h4l.4-2.5a7.4 7.4 0 0 0 2.6-1.5l2.4 1 2-3.4Z"/></svg>';

const filled = (cell: Cell) => cell.kind === "math" && ("rows" in cell.written ? cell.written.rows.some((r) => r.trim()) : cell.written.text.trim() !== "");

export class Cards {
  readonly root = el("div", { class: "cards" });
  private shown = new Map<Cell, Shown>();

  constructor(private doc: Doc, private readonly host: CardsHost) {
    this.build();
  }

  setDoc(doc: Doc) {
    this.doc = doc;
    this.build();
  }

  get document() { return this.doc; }

  // the calculator's commit, where this document is the calculator's line
  private get line() { return this.doc.preset === "calculator" ? this.host.commit : undefined; }

  // the answer a card shows now, if any
  answerOf(cell: Cell): Answer | null { return this.shown.get(cell)?.answer ?? null; }

  // every card drawn again from the document; answers come as they are asked
  build(focus?: Cell, row = 0) {
    this.shown.clear();
    this.root.replaceChildren();
    if (!this.doc.cells.length) this.doc.cells.push(math({ rows: [""] }));
    for (const cell of this.doc.cells) this.root.append(this.card(cell));
    if (presets[this.doc.preset].may.add !== "none") this.root.append(this.adder());
    if (focus) this.focus(focus, row);
  }

  // an answer arrived for a card, or for a graph the page adds after the card `from`
  show(cell: Cell, answer: Answer | null, above = 0, from?: Cell) {
    if (from) {
      const slot = this.shown.get(from)?.box.querySelector(".also");
      // closed unless the person keeps plots open; drawn the first time it opens
      const body = el("div", { class: "sec-body" });
      if (slot && answer) slot.replaceChildren(foldable("Plot", null, false, body,
        () => body.append(graphView(answer, { height: 260, curve: this.host.curve(cell) }))));
      return;
    }
    const s = this.shown.get(cell);
    if (!s) return;
    Object.assign(s, { answer, above, name: nameOf(answer ?? undefined) });
    const out = s.box.querySelector<HTMLElement>(".answer")!;
    s.box.classList.toggle("refused", Boolean(answer?.refused));
    // a row the core could not read, underlined where it is
    const written = cell.kind === "math" && "rows" in cell.written ? cell.written.rows : [];
    const at = written.flatMap((r, k) => (r.trim() ? [k] : []));
    const bad = typeof answer?.row === "number" ? at[answer.row - 1 - above] : -1;
    s.box.querySelectorAll(".rows > *").forEach((f, k) => f.classList.toggle("bad", k === bad));
    if (!answer) { out.replaceChildren(); return; }
    if (cell.kind === "graph") {
      out.replaceChildren(graphView(answer, { view: cell.settings.window as View | undefined, height: 260, curve: this.host.curve(cell),
        onView: (v) => { cell.settings.window = v; this.host.changed(); } }));
      return;
    }
    const typed = "text" in cell.written ? cell.written.text : "";
    out.replaceChildren(answerView(answer, typed, { steps: setting(this.doc, cell, "steps"), source: this.host.source(), ask: this.host.ask }));
  }

  // graphs the page added are dropped before a new run
  clearDerived() {
    for (const slot of this.root.querySelectorAll(".also")) slot.replaceChildren();
  }

  private card(cell: Cell): HTMLElement {
    const policy = presets[this.doc.preset];
    const box = el("article", { class: "card", "data-kind": cell.kind });
    const rows = el("div", { class: "rows" });
    const head = el("div", { class: "head" }, rows);
    const hint = el("p", { class: "hint", "aria-live": "polite" });
    head.append(hint);
    this.fill(cell, rows, hint);
    box.append(head);
    if (policy.settings.length || this.doc.cells.length > 1) {
      const panel = el("div", { class: "options", hidden: true });
      const gear = el("button", { type: "button", class: "gear", "aria-label": "This card's settings", "aria-expanded": "false" });
      gear.innerHTML = GEAR;
      gear.addEventListener("click", () => {
        const open = panel.hidden;
        panel.hidden = !open;
        gear.setAttribute("aria-expanded", String(open));
        if (open) this.fillSettings(cell, panel);
      });
      head.append(gear);
      box.append(panel);
    }
    // the calculator's line is cleared by x: what it asked goes into the history first
    if (this.line) {
      const x = el("button", { type: "button", class: "icon clear", "aria-label": "Clear the line", title: "Clear the line", onclick: () => this.line?.() });
      x.innerHTML = '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M6 6l12 12M18 6 6 18"/></svg>';
      head.append(x);
    }
    box.append(el("div", { class: "answer", "aria-live": "polite" }), el("div", { class: "also" }));
    this.shown.set(cell, { box, answer: null, above: 0, name: null });
    return box;
  }

  private fillSettings(cell: Cell, panel: HTMLElement) {
    panel.replaceChildren();
    if (cell.kind === "math")
      for (const key of presets[this.doc.preset].settings) {
        const now = setting(this.doc, cell, key);
        const group = el("span", { class: "choice", role: "group", "aria-label": CHOICES[key].title }, el("span", { text: CHOICES[key].title }));
        for (const [value, label] of CHOICES[key].values)
          group.append(el("button", { type: "button", class: "chip", "aria-pressed": String(value === now), text: label, onclick: () => {
            (cell.settings as Record<string, unknown>)[key] = value;
            this.fillSettings(cell, panel);
            // steps are the page's to show; numbers are the engine's, and ask again
            if (key === "steps") { const s = this.shown.get(cell); if (s) this.show(cell, s.answer, s.above); }
            this.host.changed();
          } }));
        panel.append(group);
      }
    if (this.doc.cells.length > 1)
      panel.append(el("button", { type: "button", class: "quiet small", text: "Remove", onclick: () => this.remove(cell) }));
  }

  private fill(cell: Cell, rows: HTMLElement, hint: HTMLElement) {
    if (cell.kind === "graph") {
      rows.append(el("span", { class: "of", text: `graph of ${cell.of}` }));
      return;
    }
    const tell = (word: string) => {
      const said = word ? this.host.hint(word) : null;
      hint.replaceChildren(...(said ? [el("b", { text: said[0] }), ` ${said[1]}. `, el("i", { text: said[2] })] : []));
    };
    if ("rows" in cell.written && this.host.editor) {
      const written = cell.written;
      rows.classList.toggle("many", written.rows.length > 1);
      written.rows.forEach((latex, k) => {
        const mf = document.createElement("math-field") as MathField;
        mf.mathVirtualKeyboardPolicy = "manual";
        mf.setAttribute("aria-label", written.rows.length > 1 ? `Row ${k + 1}` : "Formula");
        mf.value = latex;
        const look = () => tell(wordBefore(mf.getValue(0, mf.position, "latex")));
        mf.addEventListener("input", () => { written.rows[k] = mf.value; look(); this.host.changed(); });
        mf.addEventListener("selection-change", look);
        mf.addEventListener("focusin", () => this.host.focused(mf));
        mf.addEventListener("focusout", () => tell(""));
        mf.addEventListener("keydown", (e) => this.rowKeys(e, cell, k, mf.value), { capture: true });
        mf.addEventListener("paste", (e) => this.pasted(e, cell, k), { capture: true });
        rows.append(mf);
      });
      return;
    }
    // a text in the language, or rows of a formula while the editor is not in
    const t = el("textarea", { class: "text", rows: "1", spellcheck: "false", autocapitalize: "off", "aria-label": "Text in the language",
      placeholder: "2x^2 - 3x - 5 = 0" });
    t.value = "text" in cell.written ? cell.written.text : cell.written.rows.join("\n");
    const look = () => tell(wordBefore(t.value.slice(0, t.selectionStart ?? t.value.length)));
    t.addEventListener("input", () => { cell.written = { text: t.value }; look(); this.host.changed(); });
    t.addEventListener("keyup", look);
    t.addEventListener("blur", () => tell(""));
    t.addEventListener("keydown", (e) => {
      if (this.line && e.key === "Enter" && !e.shiftKey) { e.preventDefault(); this.line(); return; }
      // Enter on an empty last line, or with Ctrl, opens a card below; else it is a new line
      if (e.key !== "Enter" || e.shiftKey || !(e.ctrlKey || e.metaKey || t.value.endsWith("\n"))) return;
      e.preventDefault();
      cell.written = { text: t.value.trimEnd() };
      this.insertAfter(cell, [math({ text: "" })]);
    });
    rows.append(t);
  }

  // Enter writes a row below in the card; on an empty last row, or with Ctrl, it opens a card
  // below. Backspace in an empty row removes it, and in a card's only row the card.
  private rowKeys(e: KeyboardEvent, cell: Cell, k: number, value: string) {
    if (e.defaultPrevented || cell.kind !== "math" || !("rows" in cell.written)) return;
    const rows = cell.written.rows, empty = !value.trim();
    // the calculator: Enter is done, the line goes into the history; Shift+Enter writes a row below
    if (this.line && e.key === "Enter") {
      e.preventDefault();
      if (!e.shiftKey) { this.line(); return; }
      rows.splice(k + 1, 0, "");
      this.build(cell, k + 1);
      this.host.changed();
      return;
    }
    if (e.key === "Enter" && !e.shiftKey) {
      e.preventDefault();
      if (e.ctrlKey || e.metaKey || (empty && k === rows.length - 1 && rows.length > 1)) {
        if (empty && rows.length > 1) rows.splice(k, 1);
        this.insertAfter(cell, [math({ rows: [""] })]);
      } else {
        rows.splice(k + 1, 0, "");
        this.build(cell, k + 1);
        this.host.changed();
      }
    } else if (e.key === "Backspace" && value === "") {
      if (rows.length > 1) {
        e.preventDefault();
        rows.splice(k, 1);
        this.build(cell, Math.max(0, k - 1));
        this.host.changed();
      } else if (this.doc.cells.length > 1) {
        e.preventDefault();
        this.remove(cell);
      }
    }
  }

  // a pasted text of several lines -- a system, a page of definitions -- is rows of this card
  private pasted(e: ClipboardEvent, cell: Cell, k: number) {
    if (!presets[this.doc.preset].auto.includes("split") || cell.kind !== "math" || !("rows" in cell.written)) return;
    const lines = splitLines(e.clipboardData?.getData("text/plain") ?? "");
    if (lines.length < 2) return;
    e.preventDefault();
    e.stopPropagation();
    const rows = cell.written.rows;
    rows.splice(rows[k].trim() ? k + 1 : k, rows[k].trim() ? 0 : 1, ...lines);
    this.build(cell, rows.length - 1);
    this.host.changed();
  }

  private adder(): HTMLElement {
    const row = el("div", { class: "adder" });
    const add = (label: string, make: () => Cell) =>
      row.append(el("button", { type: "button", class: "quiet", text: label, onclick: () => {
        const cell = make();
        this.doc.cells.push(cell);
        this.build(cell);
        this.host.changed();
      } }));
    add("+ Formula", () => math({ rows: [""] }));
    if (presets[this.doc.preset].may.add === "any") {
      add("+ Text", () => math({ text: "" }));
      const last = [...this.shown.values()].map((s) => s.name).filter((n): n is string => Boolean(n)).at(-1);
      if (last) add(`+ Graph of ${last}`, () => graph(last));
    }
    return row;
  }

  // cells inserted after `cell`, the focus on the first; an empty last card is taken by them
  insertAfter(cell: Cell, cells: Cell[]) {
    const at = this.doc.cells.indexOf(cell);
    this.doc.cells.splice(at + 1, 0, ...cells);
    this.build(cells[0]);
    this.host.changed();
  }

  // cells added at the end, an empty last card taken
  append(cells: Cell[]) {
    const last = this.doc.cells.at(-1);
    if (last && !filled(last) && last.kind === "math") this.doc.cells.pop();
    this.doc.cells.push(...cells);
    this.build(cells.at(-1));
    this.host.changed();
  }

  private remove(cell: Cell) {
    const at = this.doc.cells.indexOf(cell);
    if (at < 0) return;
    this.doc.cells.splice(at, 1);
    this.build(this.doc.cells[Math.max(0, at - 1)]);
    this.host.changed();
  }

  focus(cell: Cell, row = 0) {
    const f = this.shown.get(cell)?.box.querySelectorAll<HTMLElement>(".rows > math-field, .rows > textarea")[row];
    if (f) requestAnimationFrame(() => f.focus());
  }
}
