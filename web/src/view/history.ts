// The calculator's history: what was asked, newest first, each a cell folded to its question and
// its result in one line. Opening one asks the core for its whole sheet and shows it under it, one
// open at a time; its question can be put back in the line, or the entry removed, with an undo.

import type { Answer } from "../core.js";
import type { Calc, Entry } from "../calc.js";
import { answerView } from "./answer.js";
import { el, math } from "./dom.js";

export interface HistoryHost {
  source: () => boolean;                          // formulas as LaTeX source
  sheet(entry: Entry): Promise<Answer | null>;    // the entry's whole answer, asked again
  use(entry: Entry): void;                        // its question in the line
  changed(calc: Calc): void;                      // the history changed: keep it
  ask: (text: string, rows?: string[]) => void;   // a question a sheet offers, as its text and rows
  undo(message: string, restore: () => void): void;
}

const ICON = {
  up: '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M12 19V5M5 12l7-7 7 7"/></svg>',
  x: '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M6 6l12 12M18 6 6 18"/></svg>',
};

export class History {
  readonly root = el("section", { class: "history", "aria-label": "History" });
  private open: Entry | null = null;

  constructor(private calc: Calc, private readonly host: HistoryHost) {}

  // the history drawn again; `opened` is shown open, as the question just asked is
  render(calc: Calc, opened: Entry | null = null) {
    this.calc = calc;
    this.open = opened;
    this.root.replaceChildren();
    if (!calc.history.length) return;
    const clear = el("button", { type: "button", class: "quiet small", text: "Clear history", onclick: () => {
      const was = this.calc;
      this.set({ ...this.calc, history: [] });
      this.host.undo("History cleared", () => this.set(was));
    } });
    this.root.append(el("div", { class: "history-head" }, el("h2", { text: "History" }), clear));
    for (const e of calc.history) this.root.append(this.cell(e));
  }

  private set(calc: Calc) {
    this.render(calc);
    this.host.changed(calc);
  }

  private cell(e: Entry): HTMLElement {
    const tex = (latex: string, cls = "math") => math(latex, { cls, source: this.host.source() });
    const rows = "rows" in e.written ? e.written.rows.filter((r) => r.trim()) : [];
    const question = rows.length ? el("span", { class: "q" }, ...rows.map((r) => tex(r))) : el("span", { class: "q mono", text: "text" in e.written ? e.written.text : "" });
    const icon = (name: keyof typeof ICON, label: string, act: () => void) => {
      const b = el("button", { type: "button", class: "icon dense", "aria-label": label, onclick: act });
      b.innerHTML = ICON[name];
      return b;
    };
    const ask = el("div", { class: "ask" }, el("span", { class: "prompt", "aria-hidden": "true", text: "›" }), question,
      el("span", { class: "acts" },
        icon("up", "Put in the line", () => this.host.use(e)),
        icon("x", "Remove from history", () => {
          const was = this.calc;
          this.set({ ...this.calc, history: this.calc.history.filter((x) => x !== e) });
          this.host.undo("Removed from history", () => this.set(was));
        })));
    const body = el("div", { class: "sheet-body" });
    const result = el("button", { type: "button", class: "result-row", "aria-expanded": "false" },
      e.result ? tex(e.result, "math res") : null,
      el("span", { class: e.ok ? "strength ok" : "strength", text: e.strength }), e.note ? el("span", { class: "note", text: e.note }) : null);
    const box = el("article", { class: "hcell cell" }, ask, result, body);
    const show = async () => {
      result.setAttribute("aria-expanded", "true");
      box.classList.add("open");
      body.replaceChildren(el("p", { class: "state body", text: "Asking again" }));
      const a = await this.host.sheet(e);
      if (this.open !== e) return;
      body.replaceChildren(a ? answerView(a, "", { steps: true, source: this.host.source(), ask: this.host.ask }) : el("p", { class: "state", text: "Not answered" }));
    };
    result.addEventListener("click", () => {
      const was = this.open;
      this.render(this.calc, was === e ? null : e);
    });
    if (this.open === e) void show();
    return box;
  }
}
