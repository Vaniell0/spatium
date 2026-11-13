// The Riemann zeta page: Hardy's Z(t) on the critical line with the zeros the core proved and the
// Gram points, zeta over a strip of the plane as colour, and the notes the library writes about the
// window. Every number is the core's (`zeta_request` in the worker); the page only draws.

import { Core, type Answer } from "./core.js";
import { $, el, math } from "./view/dom.js";

interface Zero { lo: number; hi: number; text: string; digits: number }
interface Zeros { zeros: Zero[]; counted: [number, number] | null; expected: number | null; complete: boolean }
interface Gram { n: number; t: number; good: boolean | null }
interface Note { title: string; text: string; strength: string; anchor: { t?: number; t0?: number; t1?: number } }

const PRESETS: [string, number, number][] = [["First zeros", 0, 50], ["Both sides of 0", -40, 40], ["Gram failure", 280, 284], ["Lehmer pair", 7004.5, 7005.5]];
const SVG = "http://www.w3.org/2000/svg";

const core = new Core({ glue: "../../wasm/spatium_core.js", wasm: "../../wasm/spatium_core.wasm", worker: "../../core-worker.js", cache: false });
const keys = { hardy: {}, zeros: {}, gram: {}, grid: {}, notes: {} };
const ask = async <T>(key: object, request: object): Promise<T | null> => {
  const a: Answer | null = await core.ask(key, "zeta", JSON.stringify(request));
  return a && !("refused" in a) ? (a as T) : null;
};

let win: [number, number] = [0, 50];
let asked = 0;

function svg<K extends keyof SVGElementTagNameMap>(tag: K, attrs: Record<string, string | number> = {}): SVGElementTagNameMap[K] {
  const e = document.createElementNS(SVG, tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, String(v));
  return e;
}

// ── Z(t): the curve, the zeros on the axis, the Gram points under it ──
function drawHardy(box: HTMLElement, t: number[], z: (number | null)[], zeros: Zeros | null, gram: Gram[]) {
  const w = box.clientWidth || 800, h = 240, pad = 32, axis = h - pad;
  const zs = z.filter((v): v is number => v !== null).map(Math.abs).sort((a, b) => a - b);
  const top = Math.max(1e-9, zs[Math.floor(zs.length * 0.98)] ?? 1) * 1.1;
  const x = (v: number) => ((v - win[0]) / (win[1] - win[0])) * w;
  const y = (v: number) => (axis / 2) - (v / top) * (axis / 2 - 8);
  const s = svg("svg", { viewBox: `0 0 ${w} ${h}`, width: w, height: h, class: "zplot", role: "img", "aria-label": "Hardy's Z(t)" });
  s.append(svg("line", { x1: 0, x2: w, y1: y(0), y2: y(0), class: "zaxis" }));
  let d = "";
  t.forEach((tv, i) => { const v = z[i]; if (v !== null) d += `${d ? "L" : "M"}${x(tv).toFixed(1)},${y(Math.max(-top, Math.min(top, v))).toFixed(1)}`; });
  s.append(svg("path", { d, class: "zcurve" }));
  for (const g of gram) {
    const gx = x(g.t);
    const tick = svg("line", { x1: gx, x2: gx, y1: axis + 4, y2: axis + 12, class: g.good === false ? "gram bad" : "gram" });
    tick.append(Object.assign(svg("title"), { textContent: `g${g.n} ≈ ${g.t.toFixed(6)}${g.good === false ? ", Gram's law fails" : ""}` }));
    s.append(tick);
    if (g.good === false) s.append(Object.assign(svg("text", { x: gx, y: axis + 26, class: "gram-label" }), { textContent: `g${g.n}` }));
  }
  for (const r of zeros?.zeros ?? []) {
    const zx = x((r.lo + r.hi) / 2), zy = y(0);
    const m = svg("path", { d: `M${zx},${zy - 7}L${zx + 6},${zy + 4}L${zx - 6},${zy + 4}Z`, class: "zzero" });
    m.append(Object.assign(svg("title"), { textContent: `t ∈ [${r.lo}, ${r.hi}]${zeros?.complete ? ", simple" : ""}` }));
    s.append(m);
  }
  box.replaceChildren(s);
}

// ── zeta over sigma in [-1, 2] and the window: hue the argument, contours where |zeta| is a power of 2 ──
function oklch(l: number, c: number, hDeg: number): [number, number, number] {
  const h = (hDeg * Math.PI) / 180, a = c * Math.cos(h), b = c * Math.sin(h);
  const l_ = (l + 0.3963377774 * a + 0.2158037573 * b) ** 3, m_ = (l - 0.1055613458 * a - 0.0638541728 * b) ** 3,
    s_ = (l - 0.0894841775 * a - 1.291485548 * b) ** 3;
  const lin = [4.0767416621 * l_ - 3.3077115913 * m_ + 0.2309699292 * s_, -1.2684380046 * l_ + 2.6097574011 * m_ - 0.3413193965 * s_,
    -0.0041960863 * l_ - 0.7034186147 * m_ + 1.707614701 * s_];
  return lin.map((v) => 255 * Math.max(0, Math.min(1, v <= 0.0031308 ? 12.92 * v : 1.055 * v ** (1 / 2.4) - 0.055))) as [number, number, number];
}

function drawMap(canvas: HTMLCanvasElement, w: number, h: number, re: (number | null)[], im: (number | null)[]) {
  const img = new ImageData(w, h);
  for (let j = 0; j < h; ++j)
    for (let i = 0; i < w; ++i) {
      const k = j * w + i, p = ((h - 1 - j) * w + i) * 4, a = re[k], b = im[k];
      if (a === null || b === null) { img.data.set([0, 0, 0, 255], p); continue; }
      const mod = Math.hypot(a, b), band = ((Math.log2(mod) % 1) + 1) % 1;
      const l = 0.72 * (band < 0.08 ? 0.94 : 1) * (mod / (mod + 0.08)) ** 0.5;
      const [r, g, bl] = oklch(l, 0.11, (Math.atan2(b, a) * 180) / Math.PI);
      img.data.set([r, g, bl, 255], p);
    }
  canvas.width = w, canvas.height = h;
  canvas.getContext("2d")!.putImageData(img, 0, 0);
}

function mapOverlay(box: HTMLElement, zeros: Zeros | null) {
  const s = svg("svg", { viewBox: "0 0 300 300", preserveAspectRatio: "none", class: "zoverlay", "aria-hidden": "true" });
  const x = (sigma: number) => ((sigma + 1) / 3) * 300, y = (t: number) => 300 - ((t - win[0]) / (win[1] - win[0])) * 300;
  s.append(svg("line", { x1: x(0.5), x2: x(0.5), y1: 0, y2: 300, class: "crit" }), svg("line", { x1: x(1), x2: x(1), y1: 0, y2: 300, class: "edge" }));
  if (win[0] <= 0 && 0 <= win[1]) s.append(svg("path", { d: `M${x(1) - 4},${y(0) - 4}l8,8m0,-8l-8,8`, class: "pole" }));
  for (const r of zeros?.zeros ?? []) s.append(svg("circle", { cx: x(0.5), cy: y((r.lo + r.hi) / 2), r: 3, class: "ring" }));
  box.querySelector(".zoverlay")?.remove();
  box.append(s);
}

// a note's text, its \( \) pieces as formulas
function prose(text: string): HTMLElement {
  const p = el("p", { class: "note-text" });
  text.split(/\\\((.*?)\\\)/g).forEach((piece, i) => p.append(i % 2 ? math(piece) : piece));
  return p;
}

async function show() {
  const mine = ++asked;
  const [t0, t1] = win;
  location.replace(`#t=${t0},${t1}`);
  for (const i of document.querySelectorAll<HTMLInputElement>(".zwin input")) i.value = String(i.name === "t0" ? t0 : t1);
  const plot = $("#hardy")!, map = $("#map")!, notes = $("#notes")!;
  const gist = (id: string, text: string) => { $(`#${id} .gist`)!.textContent = text; };
  gist("v-hardy", "Computing Z(t)");
  gist("v-map", "Computing ζ");
  notes.replaceChildren(el("li", { class: "state", text: "The library is writing the notes" }));
  const n = Math.min(2048, Math.max(256, Math.round((plot.clientWidth || 800) * 1.5)));
  const [hardy, gram] = await Promise.all([ask<{ t: number[]; z: (number | null)[] }>(keys.hardy, { what: "hardy", t0, t1, n }),
    ask<{ gram: Gram[] }>(keys.gram, { what: "gram", t0, t1 })]);
  if (mine !== asked || !hardy) return;
  drawHardy(plot, hardy.t, hardy.z, null, gram?.gram ?? []);
  gist("v-hardy", "Proving the zeros");
  const zeros = await ask<Zeros>(keys.zeros, { what: "zeros", t0, t1 });
  if (mine !== asked) return;
  drawHardy(plot, hardy.t, hardy.z, zeros, gram?.gram ?? []);
  const count = zeros?.zeros.length ?? 0;
  gist("v-hardy", zeros?.complete ? `${count} zero${count === 1 ? "" : "s"} proved on [${t0}, ${t1}], all there are` :
    `${count} zero${count === 1 ? "" : "s"} found on [${t0}, ${t1}]; their count not proved`);
  const canvas = map.querySelector("canvas")!;
  for (const size of [64, 256]) {
    const g = await ask<{ re: (number | null)[]; im: (number | null)[] }>(keys.grid, { what: "grid", s0: [-1, t0], s1: [2, t1], w: size, h: size });
    if (mine !== asked) return;
    if (g) drawMap(canvas, size, size, g.re, g.im);
    if (g) mapOverlay(map, zeros);
    if (size === 64) gist("v-map", "Refining");
  }
  gist("v-map", "s = σ + it: σ = Re s across, from −1 to 2, t = Im s up; hue is arg ζ, rings where |ζ| is a power of 2");
  const list = await ask<{ notes: Note[] }>(keys.notes, { what: "notes", t0, t1 });
  if (mine !== asked) return;
  notes.replaceChildren(...(list?.notes ?? []).map((nt, i) => el("li", { class: "note" },
    el("span", { class: "badge", text: String(i + 1) }),
    el("div", {}, el("h3", { text: nt.title }), prose(nt.text),
      el("span", { class: nt.strength === "proved here" ? "strength ok" : "strength", text: nt.strength })))));
  const zl = $("#zero-list")!;
  zl.replaceChildren(...(zeros?.zeros ?? []).map((r) => el("li", { class: "mono", text: `t ≈ ${r.text}` })));
}

function main() {
  const m = /#t=(-?[\d.]+),(-?[\d.]+)/.exec(location.hash);
  if (m && +m[1] < +m[2]) win = [+m[1], +m[2]];
  const form = $<HTMLFormElement>(".zwin")!;
  form.addEventListener("submit", (e) => {
    e.preventDefault();
    const a = +form.t0.value, b = +form.t1.value;
    if (a >= -100000 && a < b && b <= 100000) { win = [a, b]; void show(); }
  });
  $(".zpresets")!.append(...PRESETS.map(([name, a, b]) => el("button", { type: "button", class: "quiet dense", text: name,
    onclick: () => { win = [a, b]; void show(); } })));
  void core.load().then(() => show());
}

main();
