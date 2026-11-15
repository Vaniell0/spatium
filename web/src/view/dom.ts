// Small helpers the views share: elements, formulas through KaTeX (loaded by the page), numbers.

declare const katex: { render(latex: string, into: HTMLElement, options: object): void } | undefined;

export const $ = <T extends Element = HTMLElement>(s: string, root: ParentNode = document) => root.querySelector<T>(s);

type Child = Node | string | null | undefined | false;

export function el<K extends keyof HTMLElementTagNameMap>(
  tag: K, props: Record<string, unknown> = {}, ...children: Child[]): HTMLElementTagNameMap[K] {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(props)) {
    if (v === undefined || v === null || v === false) continue;
    if (k === "text") e.textContent = String(v);
    else if (k === "class") e.className = String(v);
    else if (k.startsWith("on") && typeof v === "function") e.addEventListener(k.slice(2), v as EventListener);
    else e.setAttribute(k, String(v));
  }
  for (const c of children) if (c !== null && c !== undefined && c !== false) e.append(c);
  return e;
}

// The page's icons, one path each, drawn at 20 px by stroke: sized by attribute as well as by
// CSS, so one inside a shrinking flex item keeps its size.
const ICONS = {
  link: "M10 14a4.5 4.5 0 0 0 6.4 0l3-3a4.5 4.5 0 0 0-6.4-6.4l-1 1M14 10a4.5 4.5 0 0 0-6.4 0l-3 3a4.5 4.5 0 0 0 6.4 6.4l1-1",
  check: "M5 12.5l4.5 4.5L19 7.5",
  keyboard: "M4.5 6h15a2 2 0 0 1 2 2v8a2 2 0 0 1-2 2h-15a2 2 0 0 1-2-2V8a2 2 0 0 1 2-2ZM6 10h1M9.5 10h1M13 10h1M16.5 10h1M7 14h10",
  sliders: "M4 7h16M4 12h16M4 17h16M9 5v4M15 10v4M10 15v4",
  book: "M4 5.5A2.5 2.5 0 0 1 6.5 3H20v15H6.5A2.5 2.5 0 0 0 4 20.5ZM4 20.5A2.5 2.5 0 0 0 6.5 23H20v-5M8 7.5h8M8 11h6",
  x: "M6 6l12 12M18 6 6 18",
  up: "M12 19V5M5 12l7-7 7 7",
  plus: "M12 5v14M5 12h14",
};
export function icon(name: keyof typeof ICONS): SVGElement {
  const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  for (const [k, v] of Object.entries({ viewBox: "0 0 24 24", width: "20", height: "20", "aria-hidden": "true" })) svg.setAttribute(k, v);
  const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
  path.setAttribute("d", ICONS[name]);
  svg.append(path);
  return svg;
}

// A formula: rendered by KaTeX, or its LaTeX source when the page shows sources or KaTeX is not in.
export function math(latex: string, { display = false, cls = "math", source = false } = {}): HTMLElement {
  const box = el(display ? "div" : "span", { class: cls });
  if (source || typeof katex === "undefined") {
    box.append(el("code", { class: "src", text: latex }));
    return box;
  }
  try {
    katex.render(latex, box, { displayMode: display, throwOnError: false, output: "htmlAndMathml" });
  } catch {
    box.append(el("code", { class: "src", text: latex }));
  }
  return box;
}

// a number to `digits` significant digits, with the shortest form that keeps them
export function number(v: unknown, digits = 12): string {
  if (typeof v !== "number" || !isFinite(v)) return String(v);
  if (v === 0) return "0";
  return String(Number(v.toPrecision(digits)));
}

