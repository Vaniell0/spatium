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

