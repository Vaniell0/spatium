// An answer as the core writes it, shown: the sheet every answer carries (`sheet`, rsc/spatium/
// respond.hpp) -- how strongly it holds, the result in one line, and its sections in the core's
// order, each a list of steps, of facts or a table. Nothing here asks which route answered,
// computes, or writes a formula.

import type { Answer } from "../core.js";
import { el, math } from "./dom.js";

export interface Shown {
  steps: boolean;
  source: boolean;   // formulas as LaTeX source
}

interface Step { label: string; latex: string }
interface Fact { label: string; text?: string; latex?: string; ok?: boolean }
interface Table { columns: string[]; rows: string[][] }
interface Section { title: string; open: boolean; echo?: boolean; note?: string; rows?: string[]; steps?: Step[]; facts?: Fact[]; table?: Table }
export interface Sheet { outcome: string; guarantee: string; label: string; ok: boolean; result: string; sections: Section[] }

// what an answer is and how strongly it holds, in the core's words; a card asked in code and a
// scene, which carry no sheet, by what they are
export function kindOf(a: Answer): [ok: boolean, label: string] {
  const s = a.sheet as Sheet | undefined;
  if (s) return [s.ok, s.label];
  if (a.route === "scene") return [false, "A scene, as LaTeX"];
  if (a.refused) return [false, "Not answered"];
  const c = a.certificate as { error_below?: number } | undefined;
  return c ? [true, "Certified"] : [false, String(a.route ?? "")];
}

// a refusal: what was not answered and why, with the place in a line of text when the reason names a column
export function refusal(a: Answer, typed: string): HTMLElement {
  const why = String(a.refused);
  const box = el("div", { class: "refusal", role: "status" }, el("p", { class: "what", text: kindOf(a)[1] }), el("p", { class: "why", text: why }));
  const col = Number((/column (\d+)/.exec(why) ?? [])[1] ?? 0);
  if (typed && col > 0 && !typed.includes("\n")) {
    const at = Math.min(col - 1, typed.length);
    box.append(el("div", { class: "caret" }, typed.slice(0, at), el("mark", { text: typed.slice(at, at + 1) || " " }), typed.slice(at + 1)));
  }
  return box;
}

export function answerView(a: Answer, typed: string, shown: Shown): DocumentFragment {
  const out = document.createDocumentFragment();
  const tex = (latex: string, cls = "math", display = false) => math(latex, { display, cls, source: shown.source });
  const s = a.sheet as Sheet | undefined;
  if (!s) {
    if (a.refused) out.append(refusal(a, typed));
    else if (typeof a.latex === "string") out.append(tex(a.latex, "big", true));
    return out;
  }
  if (s.outcome === "refused") out.append(refusal(a, typed));
  else {
    out.append(el("p", { class: s.ok ? "kind ok" : "kind", text: s.label }));
    if (s.result) out.append(tex(s.result, "big accent", true));
  }
  for (const sec of s.sections) {
    // the question again, which the cell shows above its answer; steps where the settings hide them
    if (sec.echo || (sec.steps && !shown.steps)) continue;
    out.append(section(sec, tex));
  }
  return out;
}

function section(sec: Section, tex: (latex: string, cls?: string, display?: boolean) => HTMLElement) {
  const body = el("div", { class: "sec-body" });
  if (sec.rows) for (const r of sec.rows) body.append(tex(r, "math", true));
  if (sec.steps) {
    const list = el("ol", { class: "chain", "aria-label": sec.title });
    sec.steps.forEach((st, k) => list.append(el("li", { class: k === sec.steps!.length - 1 ? "pay" : "" },
      el("span", { class: "n", text: String(k + 1) }), el("span", { class: "lab", text: st.label }), tex(st.latex, "math", true))));
    body.append(list);
  }
  if (sec.facts) {
    const dl = el("dl", { class: "facts" });
    for (const f of sec.facts)
      dl.append(el("dt", { text: f.label }),
        el("dd", { class: f.ok === true ? "ok" : f.ok === false ? "no" : "" }, f.latex !== undefined ? tex(f.latex) : el("span", { text: f.text ?? "" })));
    body.append(dl);
  }
  if (sec.table) {
    const table = el("table", { class: "sheet-table" }, el("tr", {}, ...sec.table.columns.map((c) => el("th", { text: c }))),
      ...sec.table.rows.map((r) => el("tr", {}, ...r.map((c) => el("td", { text: c })))));
    body.append(el("div", { class: "scroll" }, table));
  }
  if (sec.note) body.append(el("p", { class: "note", text: sec.note }));
  const d = el("details", { class: "sec" }, el("summary", { text: sec.title }), body);
  if (sec.open) d.setAttribute("open", "");
  return d;
}
