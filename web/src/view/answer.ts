// An answer as the core writes it, shown: the sheet every answer carries (`sheet`, rsc/spatium/
// respond.hpp) -- how strongly it holds, the result in one line, and its sections in the core's
// order, each a list of steps, of facts or a table. Nothing here asks which route answered,
// computes, or writes a formula.

import type { Answer } from "../core.js";
import { el, math } from "./dom.js";
import type { Relation } from "../gen/language.js";

export interface Shown {
  steps: boolean;
  source: boolean;                // formulas as LaTeX source
  ask?: (text: string) => void;   // a question the sheet offers, asked
}

interface Step { label: string; latex: string; lead?: string; note?: string }
interface Fact { label: string; text?: string; latex?: string; ok?: boolean }
interface Table { columns: string[]; rows: string[][] }
interface Part { kind: "holds" | "fails" | "root" | "break" | "undecided"; lo: number; hi: number; text?: string }
interface Axis { relation: Relation; line: boolean; parts: Part[] }
interface Next { caption: string; text: string; latex: string }
interface Section { title: string; open: boolean; gist?: string; note?: string; rows?: string[]; steps?: Step[]; facts?: Fact[]; table?: Table; axis?: Axis; next?: Next[] }
interface RouteChain { route: string; cost?: number; taken: boolean; gist: string; steps: Step[] }
export interface Sheet { outcome: string; guarantee: string; verdict: { strength: string; notes: string[] }; label: string; ok: boolean; result: string;
  sections: Section[]; routes?: RouteChain[] }

// what a person chose before, kept in this browser only: which sections they keep open, and the
// route whose chain they read first; a page with no storage shows the core's defaults
function kept<T>(key: string): Record<string, T> {
  try { return JSON.parse(localStorage.getItem(key) ?? "{}") as Record<string, T>; } catch { return {}; }
}
function keep<T>(key: string, name: string, value: T) {
  try { localStorage.setItem(key, JSON.stringify({ ...kept<T>(key), [name]: value })); } catch { /* not kept */ }
}

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
    else if (typeof a.latex === "string") out.append(el("div", { class: "result" }, tex(a.latex, "math", true)));
    return out;
  }
  // the result, then how strongly it holds and how: "exact  two real roots  checked by two routes"
  if (s.outcome === "refused") out.append(el("div", { class: "result" }, refusal(a, typed)));
  else out.append(el("div", { class: "result" }, s.result ? tex(s.result, "math", true) : null,
    el("p", { class: s.ok ? "verdict ok" : "verdict" }, el("span", { class: "strength", text: s.verdict.strength }),
      ...s.verdict.notes.map((n) => el("span", { text: n })))));
  for (const sec of s.sections) {
    // the plot, which the page draws, goes where the sheet's order puts it: before how it was found
    if (sec.title === "How it was found") out.append(el("div", { class: "plot-slot" }));
    // steps where the settings hide them
    if (sec.steps && !shown.steps) continue;
    out.append(section(sec, tex, sec.title === "Steps" ? s.routes ?? [] : [], shown.ask));
  }
  return out;
}

type Tex = (latex: string, cls?: string, display?: boolean) => HTMLElement;

function chain(steps: Step[], title: string, tex: Tex) {
  // numbered by the list itself; what the step does above, its formula under it, the last one paid
  const list = el("ol", { class: "chain", "aria-label": title });
  steps.forEach((st, k) => list.append(el("li", { class: k === steps.length - 1 ? "step pay" : "step" },
    // what the step does: its lead in bold, the rest beside it; a step with neither says its label
    el("span", { class: "op" }, el("b", { text: st.lead ?? st.label.charAt(0).toUpperCase() + st.label.slice(1) }), st.note ? ` ${st.note}` : null),
    tex(st.latex, "math", true))));
  return list;
}

// a section of the sheet; its open state as the person left it, else the core's. Steps with more
// than one route answering carry a tab each, the route read before first, else the one taken.
export function section(sec: Section, tex: Tex, routes: RouteChain[] = [], ask?: (text: string) => void) {
  const body = el("div", { class: "body" });
  const gist = el("span", { class: "tag", text: sec.gist ?? "" });
  if (sec.rows) for (const r of sec.rows) body.append(tex(r, "math", true));
  if (sec.steps && routes.length > 1) {
    const before = kept<number>("spatium.route");
    const first = routes.reduce((b, r, k) => ((before[r.route] ?? 0) > (before[routes[b].route] ?? 0) ? k : b), Math.max(0, routes.findIndex((r) => r.taken)));
    const tabs = el("div", { class: "route-tabs", role: "tablist", "aria-label": "Routes" });
    const panel = el("div", { role: "tabpanel" });
    const pick = (k: number, remember: boolean) => {
      tabs.querySelectorAll("button").forEach((b, j) => b.setAttribute("aria-selected", String(j === k)));
      panel.replaceChildren(chain(routes[k].steps, routes[k].route, tex));
      gist.textContent = `${routes[k].route}, ${routes[k].gist}`;
      if (remember) keep("spatium.route", routes[k].route, Date.now());
    };
    routes.forEach((r, k) => tabs.append(el("button", { type: "button", role: "tab", onclick: () => pick(k, true) },
      el("span", { text: r.route }), r.cost !== undefined ? el("span", { class: "cost", text: String(r.cost) }) : null)));
    body.append(tabs, panel);
    pick(first, false);
  } else if (sec.steps) body.append(chain(sec.steps, sec.title, tex));
  if (sec.facts) {
    const dl = el("dl", { class: "facts" });
    for (const f of sec.facts)
      dl.append(el("dt", { text: f.label }),
        el("dd", { class: f.ok === true ? "agree" : f.ok === false ? "no" : "" }, f.latex !== undefined ? tex(f.latex) : el("span", { text: f.text ?? "" })));
    body.append(dl);
  }
  if (sec.table) {
    // the routes tried, the taken one marked; any other table is numbers, in columns
    const how = sec.title === "How it was found";
    const table = el("table", { class: how ? "" : "roots" }, el("tr", {}, ...sec.table.columns.map((c) => el("th", { text: c }))),
      ...sec.table.rows.map((r) => el("tr", {}, ...r.map((c, j) => el("td", { class: how ? (j === 0 && c.endsWith(", taken") ? "taken" : j === 1 ? "num" : "") : "", text: c })))));
    body.append(el("div", { class: how ? "scroll how" : "scroll" }, table));
  }
  if (sec.axis) body.append(axisView(sec.axis));
  // what else the question can become: a caption and the question as the library writes it
  if (sec.next) body.append(el("ul", { class: "next" }, ...sec.next.map((n) => {
    const inner = [el("span", { class: "caption", text: n.caption }), tex(n.latex, "math")];
    return el("li", {}, ask ? el("button", { type: "button", class: "ask", onclick: () => ask(n.text) }, ...inner) : el("div", {}, ...inner));
  })));
  if (sec.note) body.append(el("p", { class: "note", text: sec.note }));
  return foldable(sec.title, gist, sec.open, body);
}

// a titled part that folds, open as the person last left one of that title, else as `open`; `fill`
// builds its body the first time it opens, so what is never opened is never built
export function foldable(title: string, gist: HTMLElement | null, open: boolean, body: HTMLElement, fill?: () => void) {
  const summary = el("summary", {}, title, gist);
  const d = el("details", { class: "sec" }, summary, body);
  const fillOnce = () => { if (d.open && fill && !body.firstChild) fill(); };
  d.addEventListener("toggle", fillOnce);
  // only a person's choice is kept: a click (or Enter) on the title, not the core's default opening it
  summary.addEventListener("click", () => keep("spatium.open", title, !d.open));
  if (kept<boolean>("spatium.open")[title] ?? open) d.setAttribute("open", "");
  fillOnce();
  return d;
}

// An inequality's axis as the interval method draws it: its points evenly apart (a scheme, not a
// scale), the pieces where it holds shaded, the sign of its left side less its right over each
// piece, a root filled where it belongs to the set and hollow where not, a break hollow.
function axisView(a: Axis): HTMLElement {
  const NS = "http://www.w3.org/2000/svg";
  const svg = (tag: string, at: Record<string, string | number>, text?: string) => {
    const e = document.createElementNS(NS, tag);
    for (const [k, v] of Object.entries(at)) e.setAttribute(k, String(v));
    if (text !== undefined) e.textContent = text;
    return e;
  };
  const points = a.parts.filter((p) => p.kind !== "holds" && p.kind !== "fails");
  const step = 120, W = step * (points.length + 1), H = 72, y = 36;
  const root = svg("svg", { viewBox: `0 0 ${W} ${H}`, width: W, height: H, role: "img", "aria-label": "The axis: where the inequality holds" });
  const marks: Element[] = [];   // the points, drawn over the pieces
  root.append(svg("line", { x1: 4, y1: y, x2: W - 4, y2: y, class: "rule" }));
  const greater = a.relation.startsWith(">"), equal = a.relation.endsWith("=");
  // each piece between its neighbouring points, or the axis' ends
  let k = 0;
  a.parts.forEach((p) => {
    if (p.kind === "holds" || p.kind === "fails") {
      const x0 = k * step + (k ? 0 : 4), x1 = (k + 1) * step;
      if (p.kind === "holds") root.append(svg("rect", { x: x0, y: y - 6, width: x1 - x0, height: 12, class: "band" }));
      root.append(svg("text", { x: (x0 + x1) / 2, y: y - 14, class: "sign" }, (p.kind === "holds") === greater ? "+" : "\u2212"));
      return;
    }
    k += 1;
    const x = k * step;
    const filled = p.kind === "root" && equal;
    marks.push(svg("circle", { cx: x, cy: y, r: 6, class: filled ? "pt in" : p.kind === "undecided" ? "pt open" : "pt out" }));
    // on the scheme an end is short: its digits are in the answer and under "Ends"
    const short = (t: string) => (t.length > 8 && !t.includes("/") ? `\u2248${Number(Number(t).toPrecision(6))}` : t);
    const label = p.kind === "undecided" ? "?" : short(p.text ?? String((p.lo + p.hi) / 2));
    marks.push(svg("text", { x, y: y + 28, class: "at" }, label));
  });
  root.append(...marks);
  if (a.line) for (const x of [12, W - 12]) root.append(svg("text", { x, y: y + 28, class: "at" }, x < W / 2 ? "\u2212\u221e" : "\u221e"));
  return el("div", { class: "axis" }, root as unknown as Node);
}
