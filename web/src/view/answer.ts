// An answer as the core writes it, shown: what kind it is, the value or the roots, the chain of
// steps, the certificate, or the refusal and where it is. Nothing here computes or writes a formula.

import type { Answer } from "../core.js";
import { decade, decadeText, el, math, number } from "./dom.js";

export interface Shown {
  digits: number;
  steps: boolean;
  source: boolean;   // formulas as LaTeX source
}

interface Step { label: string; latex: string }
interface Certificate { strength: string; norm: string; points: number; scalar: string; computed: string; error_below: number }
interface Piece { lo: number; hi: number; form: string; latex: string; degree?: number; coefficients?: number[] }

const COMPUTED: Record<string, string> = { wasm: "in this browser", host: "on the host, when the gallery was written" };

const HOW: Record<string, string> = {
  evaluate: "Evaluated as written, in double precision, with nothing rewritten.",
  "closed form": "Solved in closed form. The coefficients are read from the equation's own values and the degree is checked at more points; the roots come from the quadratic formula in the form that keeps the small root, and Vieta's formulas check them by a second route.",
  synthesis: "Synthesised. The library computes the question's value at each parameter by its own solvers; the answer is the shortest expression in the parameter that matches it, searched bottom up, or a Chebyshev polynomial where none of a few terms does, and it is checked against those values on a dense grid.",
};

export function kindOf(a: Answer): [ok: boolean, label: string] {
  if (a.route === "scene") return [false, "A scene, as LaTeX"];
  if (a.refused) return [false, a.route ? `Not answered by ${a.route}` : "Not answered"];
  if (a.route === "evaluate") return [false, "Evaluated"];
  if (a.route === "expression") return [false, `A function of ${a.reads}, as written`];
  if (a.route === "closed form") {
    const check = a.check as { holds?: boolean } | undefined;
    return check?.holds ? [true, "Closed form, checked by Vieta's formulas"] : [false, "Closed form"];
  }
  const c = a.certificate as Certificate | undefined;
  if (!c) return [false, String(a.route ?? "")];
  const e = c.error_below;
  return [true, e ? `Certified, ${c.norm === "largest relative error" ? "relative " : ""}error below ${decadeText(e)}` : "Certified, exact on the grid"];
}

function facts(rows: [string, Node][]) {
  const dl = el("dl", { class: "facts" });
  for (const [k, v] of rows) dl.append(el("dt", { text: k }), el("dd", {}, v));
  return dl;
}

// the step a closed-form answer ends in: its roots, or the root of a linear equation
export function finalStep(steps: Step[]): Step {
  const named = steps.find((s) => s.label === "roots" || s.label === "root");
  if (named) return named;
  const solving = steps.filter((s) => !s.label.startsWith("check"));
  return solving[solving.length - 1] ?? steps[steps.length - 1];
}

function how(route: unknown) {
  const text = typeof route === "string" ? HOW[route] : undefined;
  return text ? el("details", { class: "more" }, el("summary", { text: "How it was found" }), el("p", { text })) : null;
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
  if (a.refused) {
    out.append(refusal(a, typed));
    return out;
  }
  const [ok, label] = kindOf(a);
  out.append(el("p", { class: ok ? "kind ok" : "kind", text: label }));
  if (a.route === "scene") {
    out.append(tex(String(a.latex), "big", true));
    return out;
  }
  // a parameter: the interval a later cell reads it on
  const last = (a.document as { statements?: { name: string; kind: string; lo?: string; hi?: string }[] } | undefined)?.statements?.at(-1);
  if (last?.kind === "parameter") {
    out.firstChild?.remove();
    out.append(el("p", { class: "kind", text: "Parameter" }), tex(`${last.name} \\in [${last.lo}, ${last.hi}]`, "big", true));
    return out;
  }
  if (a.route === "expression") {
    out.append(tex(String(a.latex), "big accent", true), facts([["Reads", el("span", { class: "mono", text: String(a.reads) })]]));
    return out;
  }
  if (a.route === "evaluate") {
    out.append(tex(String(a.latex), "big", true), facts([["Value", el("span", { class: "mono", text: number(a.value, shown.digits) })]]));
    const h = how(a.route);
    if (h) out.append(h);
    return out;
  }
  if (a.route === "closed form") {
    const steps = (a.steps ?? []) as Step[];
    const last = finalStep(steps);
    out.append(tex(last.latex, "big accent", true));
    if (shown.steps) {
      const chain = el("ol", { class: "chain", "aria-label": "Steps" });
      steps.forEach((s, k) => {
        const cls = s === last ? "pay" : s.label.startsWith("check") ? "check" : "";
        chain.append(el("li", { class: cls }, el("span", { class: "n", text: String(k + 1) }), el("span", { class: "lab", text: s.label }),
          tex(s.latex, "math", true)));
      });
      out.append(chain);
    }
    const check = a.check as { by: string; holds: boolean; error_below: number } | undefined;
    if (check) out.append(facts([
      ["Check", el("span", { class: check.holds ? "ok" : "", text: `${check.by} ${check.holds ? "hold" : "do not hold"}` })],
      [check.error_below ? "Relative error below" : "Error", decade(check.error_below)],
    ]));
    const h = how(a.route);
    if (h) out.append(h);
    return out;
  }
  if (a.route === "synthesis") {
    const pieces = (a.answer ?? []) as Piece[];
    const one = pieces.length === 1 && pieces[0].form !== "chebyshev";
    if (one) out.append(tex(pieces[0].latex, "big accent", true));
    else out.append(el("p", { class: "big", text: pieces.length === 1 ? `A Chebyshev series of degree ${pieces[0].degree}` : `${a.pieces} pieces` }));
    const c = a.certificate as Certificate | undefined;
    if (c) out.append(facts([
      ["Holds", el("span", { class: "ok", text: `${c.strength}, on ${c.points} points` })],
      [c.norm === "largest relative error" ? "Relative error below" : "Absolute error below", decade(c.error_below)],
      ["Scalar", el("span", { text: c.scalar })],
      ["Computed", el("span", { text: COMPUTED[c.computed] ?? c.computed })],
    ]));
    if (!one) {
      const list = el("div", { class: "pieces" });
      for (const p of pieces) list.append(pieceView(p, pieces.length > 1, shown, tex));
      out.append(list);
    }
    const h = how(a.route);
    if (h) out.append(h);
  }
  return out;
}

function pieceView(p: Piece, many: boolean, shown: Shown, tex: (latex: string) => HTMLElement) {
  const box = el("div", { class: "piece" }, el("p", { class: "label", text: `on [${number(p.lo, shown.digits)}, ${number(p.hi, shown.digits)}]` }));
  if (many && p.form !== "chebyshev") box.append(tex(p.latex));
  if (p.form === "chebyshev" && p.coefficients) {
    const list = el("ol", { class: "coeffs" });
    p.coefficients.forEach((c, k) => list.append(el("li", {}, el("span", { text: `c${k}` }), number(c, shown.digits))));
    box.append(el("details", {}, el("summary", { text: `Chebyshev degree ${p.degree}, ${p.coefficients.length} coefficients` }), list));
  }
  return box;
}
