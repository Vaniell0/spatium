// A graph: the fragment shader the core writes for a value along one parameter, compiled in
// WebGL2. The shader computes the value at every pixel from the view it is given (`plot_view`)
// and the other parameters (`field_symbols`), so panning, zooming and sliding a parameter are
// new uniforms, never a new question to the core. Axes and their numbers are a 2D canvas under it.
// Zoomed past what fp32 tells apart, the shader's curve is a staircase; there the core gives the
// curve in doubles, a value a column (`curve`), drawn on the 2D canvas instead.

import type { Answer } from "../core.js";
import { el, number } from "./dom.js";

export type View = [number, number, number, number];   // x0 y0 x1 y1

interface Plot {
  vertex: string;
  fragment: string;
  interval: [number, number];
  range?: [number, number];
  symbols: { name: string; interval: [number, number] }[];
  weights: number[];
}

// the view a plot opens on: its interval, and the heights it reaches with a margin
export function fitted(p: Plot): View {
  const [x0, x1] = p.interval;
  let [y0, y1] = p.range ?? [-1, 1];
  if (!(y1 > y0)) { y0 -= 1; y1 += 1; }
  const pad = (y1 - y0) * 0.08;
  return [x0, y0 - pad, x1, y1 + pad];
}

export interface GraphOptions {
  view?: View;                      // as the cell keeps it; the fitted view otherwise
  height: number;
  onView?: (view: View) => void;    // after a pan or a zoom settles
  // the value at the centres of n columns of [x0, x1], in doubles, the sliders' values given
  curve?: (at: { x0: number; x1: number; n: number; slots: number[] }) => Promise<(number | null)[] | null>;
}

// fp32's step at the numbers a view spans: a view with fewer than a few of them a pixel shows them
const fp32Step = (a: number, b: number) => Math.max(Math.abs(a), Math.abs(b)) * 2 ** -23;
export const pastFp32 = ([x0, y0, x1, y1]: View, w: number, h: number) =>
  (x1 - x0) / w < 4 * fp32Step(x0, x1) || (y1 - y0) / h < 4 * fp32Step(y0, y1);

export function graphView(answer: Answer, options: GraphOptions): HTMLElement {
  const p = answer as unknown as Plot;
  const fig = el("figure", { class: "graph" });
  if (answer.refused) {
    fig.append(el("p", { class: "state", text: String(answer.refused) }));
    return fig;
  }
  const gl = el("canvas", { class: "gl", "aria-hidden": "true" });
  const axes = el("canvas", { class: "axes", "aria-hidden": "true" });
  const plane = el("div", { class: "plane", tabindex: "0", "aria-label": "Graph: drag to pan, wheel to zoom" }, axes, gl);
  const at = el("span", { class: "at mono" });
  const fit = el("button", { type: "button", class: "quiet small", text: "Fit" });
  const exact = el("span", { class: "at", text: "computed in double", hidden: true });
  fig.append(el("div", { class: "graph-bar" }, at, exact, fit), plane);

  const ctx = gl.getContext("webgl2", { premultipliedAlpha: false, antialias: false });
  if (!ctx) {
    fig.replaceChildren(el("p", { class: "state", text: "This browser has no WebGL2 to draw the graph with." }));
    return fig;
  }
  const shader = (type: number, text: string) => {
    const s = ctx.createShader(type)!;
    ctx.shaderSource(s, text);
    ctx.compileShader(s);
    return s;
  };
  const prog = ctx.createProgram()!;
  ctx.attachShader(prog, shader(ctx.VERTEX_SHADER, p.vertex));
  ctx.attachShader(prog, shader(ctx.FRAGMENT_SHADER, p.fragment));
  ctx.linkProgram(prog);
  if (!ctx.getProgramParameter(prog, ctx.LINK_STATUS)) {
    fig.replaceChildren(el("p", { class: "state", text: `The graph's shader did not link: ${ctx.getProgramInfoLog(prog)}` }));
    return fig;
  }

  let view: View = options.view ?? fitted(p);
  const held = p.symbols.map((s) => (s.interval[0] + s.interval[1]) / 2);
  // the curve in doubles, and the view and sliders it was asked for
  let curve: { key: string; values: (number | null)[] } | null = null;
  let asking: ReturnType<typeof setTimeout> | 0 = 0;

  // the parameters the graph holds: a slider each, a uniform each
  if (p.symbols.length) {
    const sliders = el("div", { class: "sliders" });
    p.symbols.forEach((s, k) => {
      const [lo, hi] = s.interval;
      const value = el("span", { class: "mono", text: number(held[k], 6) });
      const input = el("input", { type: "range", min: String(lo), max: String(hi), step: String((hi - lo) / 1000), value: String(held[k]), "aria-label": s.name });
      input.addEventListener("input", () => {
        held[k] = Number(input.value);
        value.textContent = number(held[k], 6);
        paint();
      });
      sliders.append(el("label", { class: "slider" }, el("span", { class: "name", text: s.name }), input, value));
    });
    fig.append(sliders);
  }

  function paint() {
    const dpr = devicePixelRatio || 1, w = plane.clientWidth || 480, h = options.height;
    for (const c of [gl, axes]) {
      c.width = Math.round(w * dpr);
      c.height = Math.round(h * dpr);
      c.style.height = `${h}px`;
    }
    const css = getComputedStyle(document.documentElement);
    const ink = rgb(css.getPropertyValue("--accent"));
    ctx!.viewport(0, 0, gl.width, gl.height);
    ctx!.useProgram(prog);
    const u = (n: string) => ctx!.getUniformLocation(prog, n);
    ctx!.uniform4fv(u("plot_view"), view);
    ctx!.uniform2f(u("plot_size"), gl.width, gl.height);
    ctx!.uniform1f(u("plot_width"), 1.1 * dpr);
    ctx!.uniform4f(u("plot_ink"), ink[0], ink[1], ink[2], 1);
    ctx!.uniform1i(u("plot_values"), 0);
    if (held.length) ctx!.uniform1fv(u("field_symbols"), new Float32Array(held));
    if (p.weights.length) ctx!.uniform4fv(u("field_weights"), new Float32Array(p.weights));
    ctx!.clearColor(0, 0, 0, 0);
    ctx!.clear(ctx!.COLOR_BUFFER_BIT);
    drawAxes(axes, view, dpr, css);
    // past fp32, the doubles' curve for this view once it is here, the shader's until then
    const fine = Boolean(options.curve) && pastFp32(view, w, h);
    const key = JSON.stringify([view, held, w]);
    exact.hidden = !(fine && curve?.key === key);
    if (fine && curve?.key === key) {
      drawCurve(axes, view, curve.values, dpr, css.getPropertyValue("--accent"));
      return;
    }
    ctx!.enable(ctx!.BLEND);
    ctx!.blendFunc(ctx!.SRC_ALPHA, ctx!.ONE_MINUS_SRC_ALPHA);
    ctx!.drawArrays(ctx!.TRIANGLES, 0, 3);
    if (!fine) return;
    clearTimeout(asking);
    asking = setTimeout(async () => {
      const values = await options.curve!({ x0: view[0], x1: view[2], n: Math.round(w), slots: [...held] });
      if (values && JSON.stringify([view, held, w]) === key) { curve = { key, values }; paint(); }
    }, 150);
  }

  // where a pixel of the plane is in the view
  const toView = (e: { clientX: number; clientY: number }): [number, number] => {
    const r = plane.getBoundingClientRect();
    const fx = (e.clientX - r.left) / r.width, fy = 1 - (e.clientY - r.top) / r.height;
    return [view[0] + fx * (view[2] - view[0]), view[1] + fy * (view[3] - view[1])];
  };
  let settle: ReturnType<typeof setTimeout> | 0 = 0;
  const moved = () => {
    paint();
    clearTimeout(settle);
    settle = setTimeout(() => options.onView?.([...view]), 300);
  };
  // no closer than the doubles tell apart, and no farther than 10^15 across
  const zoom = (cx: number, cy: number, k: number) => {
    const next: View = [cx + (view[0] - cx) * k, cy + (view[1] - cy) * k, cx + (view[2] - cx) * k, cy + (view[3] - cy) * k];
    const fine = (a: number, b: number) => b - a < 64 * Math.max(Math.abs(a), Math.abs(b), 1e-300) * Number.EPSILON;
    if (k < 1 && (fine(next[0], next[2]) || fine(next[1], next[3]))) return;
    if (k > 1 && (next[2] - next[0] > 1e15 || next[3] - next[1] > 1e15)) return;
    view = next;
    moved();
  };

  // drag with one pointer pans, two pointers pinch
  const pointers = new Map<number, { x: number; y: number }>();
  plane.addEventListener("pointerdown", (e) => {
    plane.setPointerCapture(e.pointerId);
    pointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
  });
  plane.addEventListener("pointermove", (e) => {
    const [vx, vy] = toView(e);
    at.textContent = `x ${number(vx, 4)}  y ${number(vy, 4)}`;
    const last = pointers.get(e.pointerId);
    if (!last) return;
    const r = plane.getBoundingClientRect();
    if (pointers.size === 1) {
      const dx = ((e.clientX - last.x) / r.width) * (view[2] - view[0]), dy = ((e.clientY - last.y) / r.height) * (view[3] - view[1]);
      view = [view[0] - dx, view[1] + dy, view[2] - dx, view[3] + dy];
      pointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
      moved();
    } else if (pointers.size === 2) {
      const other = [...pointers.entries()].find(([id]) => id !== e.pointerId)![1];
      const before = Math.hypot(last.x - other.x, last.y - other.y), after = Math.hypot(e.clientX - other.x, e.clientY - other.y);
      pointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
      if (before > 0 && after > 0) {
        const [cx, cy] = toView({ clientX: (e.clientX + other.x) / 2, clientY: (e.clientY + other.y) / 2 });
        zoom(cx, cy, before / after);
      }
    }
  });
  const up = (e: PointerEvent) => pointers.delete(e.pointerId);
  plane.addEventListener("pointerup", up);
  plane.addEventListener("pointercancel", up);
  plane.addEventListener("pointerleave", () => { at.textContent = ""; });
  plane.addEventListener("wheel", (e) => {
    e.preventDefault();
    const [cx, cy] = toView(e);
    zoom(cx, cy, Math.exp(e.deltaY * 0.0015));
  }, { passive: false });
  plane.addEventListener("dblclick", () => { view = fitted(p); moved(); });
  fit.addEventListener("click", () => { view = fitted(p); moved(); });

  requestAnimationFrame(paint);
  new ResizeObserver(() => paint()).observe(plane);
  return fig;
}

function rgb(color: string): [number, number, number] {
  const c = document.createElement("canvas").getContext("2d")!;
  c.fillStyle = color.trim();
  const hex = String(c.fillStyle).replace("#", "");
  return [0, 2, 4].map((k) => parseInt(hex.slice(k, k + 2), 16) / 255) as [number, number, number];
}

// round numbers along an axis, about one each 80 pixels
function ticks(a: number, b: number, pixels: number): number[] {
  const raw = (b - a) / Math.max(2, pixels / 80), mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 5, 10].map((m) => m * mag).find((s) => s >= raw)!;
  const out: number[] = [];
  for (let t = Math.ceil(a / step) * step; t <= b + step * 1e-9; t += step) out.push(Math.abs(t) < step * 1e-9 ? 0 : t);
  return out;
}

// a tick's number with as many digits as the step between ticks needs: at 0.99999 zoomed to 1e-7
// a tick is 0.9999901, not 0.99999 again
function tickText(t: number, step: number): string {
  const digits = Math.ceil(Math.log10(Math.max(Math.abs(t), step) / step)) + 1;
  return String(Number(t.toPrecision(Math.min(17, Math.max(1, digits)))));
}

// the curve a value a column, a gap where there is none
function drawCurve(canvas: HTMLCanvasElement, [, y0, , y1]: View, values: (number | null)[], dpr: number, ink: string) {
  const g = canvas.getContext("2d")!, W = canvas.width, H = canvas.height;
  g.strokeStyle = ink;
  g.lineWidth = 1.5 * dpr;
  g.beginPath();
  let pen = false;
  values.forEach((v, k) => {
    if (v === null) { pen = false; return; }
    const x = ((k + 0.5) / values.length) * W, y = H - ((v - y0) / (y1 - y0)) * H;
    if (pen) g.lineTo(x, y); else g.moveTo(x, y);
    pen = true;
  });
  g.stroke();
}

function drawAxes(canvas: HTMLCanvasElement, [x0, y0, x1, y1]: View, dpr: number, css: CSSStyleDeclaration) {
  const g = canvas.getContext("2d")!, W = canvas.width, H = canvas.height;
  const X = (x: number) => ((x - x0) / (x1 - x0)) * W, Y = (y: number) => H - ((y - y0) / (y1 - y0)) * H;
  g.clearRect(0, 0, W, H);
  g.font = `${12 * dpr}px ${css.getPropertyValue("--f-mono")}`;
  g.lineWidth = dpr;
  g.strokeStyle = css.getPropertyValue("--line-faint");
  g.fillStyle = css.getPropertyValue("--ink-3");
  const xs = ticks(x0, x1, W / dpr), ys = ticks(y0, y1, H / dpr);
  for (const t of xs) {
    g.beginPath(); g.moveTo(X(t), 0); g.lineTo(X(t), H); g.stroke();
    g.fillText(tickText(t, (x1 - x0) / Math.max(1, xs.length)), X(t) + 4 * dpr, H - 6 * dpr);
  }
  for (const t of ys) {
    g.beginPath(); g.moveTo(0, Y(t)); g.lineTo(W, Y(t)); g.stroke();
    g.fillText(tickText(t, (y1 - y0) / Math.max(1, ys.length)), 4 * dpr, Y(t) - 4 * dpr);
  }
  g.strokeStyle = css.getPropertyValue("--ink-3");
  if (y0 < 0 && y1 > 0) { g.beginPath(); g.moveTo(0, Y(0)); g.lineTo(W, Y(0)); g.stroke(); }
  if (x0 < 0 && x1 > 0) { g.beginPath(); g.moveTo(X(0), 0); g.lineTo(X(0), H); g.stroke(); }
}
