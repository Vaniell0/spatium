// The page: the calculator -- one line and the history of what was asked -- or the canvas, a
// document of cells, as the bar picks, answered live by the core in a worker. It keeps what is
// written in its address (#w=) and in the browser, and shows the JSON the core and the library's
// files hold; it composes nothing.

import { Answers, notation, type Parse } from "./answers.js";
import { type Calc, type Entry, blank, commit, empty, load as loadCalc, save as saveCalc } from "./calc.js";
import { Core, type State } from "./core.js";
import { layout } from "./keys.js";
import { type Doc, type Preset, fromJSON, fromLink, math, toJSON, toLink } from "./doc.js";
import { Cards } from "./view/cells.js";
import { $, el, math as tex } from "./view/dom.js";
import type { Sheet } from "./view/answer.js";
import { History } from "./view/history.js";

// ── settings of the page (a cell's own are in the document) ─────────────

const DEFAULTS = { theme: "system", scale: "1", formulas: "rendered", keyboard: "on", cache: "on" };
type PageSettings = typeof DEFAULTS;
const settings: PageSettings = { ...DEFAULTS };
try { Object.assign(settings, JSON.parse(localStorage.getItem("spatium.settings") || "{}")); } catch {}
const saveSettings = () => { try { localStorage.setItem("spatium.settings", JSON.stringify(settings)); } catch {} };

// the one transient message: what was just undone-able, under the bar, for five seconds
let undoTimer: ReturnType<typeof setTimeout> | 0 = 0, undoing: (() => void) | null = null;
function undoStrip(message: string, restore: () => void) {
  const strip = $("#undo")!;
  clearTimeout(undoTimer);
  undoing = restore;
  strip.replaceChildren(el("span", { text: message }), el("button", { type: "button", class: "quiet small", text: "Undo", onclick: () => {
    undoing?.();
    undoing = null;
    strip.hidden = true;
  } }));
  strip.hidden = false;
  undoTimer = setTimeout(() => { strip.hidden = true; undoing = null; }, 5000);
}

let toastTimer: ReturnType<typeof setTimeout> | 0 = 0;
function toast(text: string, bad = false) {
  const t = $("#toast")!;
  clearTimeout(toastTimer);
  t.textContent = text;
  t.hidden = false;
  t.classList.toggle("bad", bad);
  if (!bad) toastTimer = setTimeout(() => { t.hidden = true; }, 3000);
}

function dot(state: State, title: string) {
  const d = $("#core-dot")!;
  d.className = `core-dot ${state === "loading" ? "" : state}`;
  d.title = title;
}

// ── the core, the editor, the library's files ───────────────────────────

const core = new Core({
  glue: "wasm/spatium_core.js", wasm: "wasm/spatium_core.wasm", worker: "core-worker.js",
  cache: settings.cache === "on", onState: dot,
});

const editor = {
  ready: false,
  parse: (() => ["Error", "'the formula editor did not load'"]) as Parse,
  toLatex: (text: string) => text,   // a line of the language as a formula (MathLive's ascii-math reader)
};

async function loadEditor() {
  try {
    const [ml, ls] = await Promise.all([
      import(new URL("../../vendor/mathlive/mathlive.min.mjs", import.meta.url).href),
      import(new URL("../../vendor/latex-syntax/latex-syntax.min.js", import.meta.url).href),
    ]);
    ml.MathfieldElement.fontsDirectory = new URL("vendor/katex/fonts/", location.href).href;
    ml.MathfieldElement.soundsDirectory = null;
    ml.MathfieldElement.computeEngine = null;
    const kb = (window as unknown as { mathVirtualKeyboard: { container: Element; layouts: unknown } }).mathVirtualKeyboard;
    kb.container = $("#keys-body")!;
    kb.layouts = [layout()];
    editor.parse = (latex: string) => ls.parse(notation(latex));
    editor.toLatex = ml.convertAsciiMathToLatex;
    editor.ready = true;
  } catch (e) {
    toast(`The formula editor did not load (${e instanceof Error ? e.message : e}); cells are written as text.`, true);
  }
}

interface Reference {
  functions: { name: string; arity: number; caption: string; domain: string; example?: string }[];
  operators: { written: string; caption: string; example?: string }[];
  constants: { name: string; caption: string; example?: string }[];
  questions: { form: string; caption: string; example?: string; answer?: { question_latex?: string[] } }[];
}
const data: { reference: Reference | null } = { reference: null };

async function loadData() {
  data.reference = await fetch("data/reference.json").then((x) => x.json()).catch(() => null);
}

// ── the document: from the address, else the browser, else a question to start from ─

const KEPT = "spatium.document";

function starting(): Doc {
  return { preset: "calculator", cells: [math({ rows: ["2x^2-3x-5=0"] })] };
}

async function fromAddress(): Promise<Doc | null> {
  const hash = location.hash.slice(1);
  if (hash.startsWith("w=")) {
    const doc = await fromLink(hash.slice(2));
    if (!doc) toast("The link to a document could not be read", true);
    return doc;
  }
  const [path, query] = hash.split("?");
  if (path === "solve" || path === "") {
    const q = new URLSearchParams(query || "").get("q");
    if (q !== null) return { preset: "calculator", cells: [math({ text: q })] };
  }
  return null;
}

function fromBrowser(): Doc | null {
  try {
    const kept = JSON.parse(localStorage.getItem(KEPT) || "null");
    return kept ? fromJSON(kept) : null;
  } catch { return null; }
}

// ── asking and keeping ──────────────────────────────────────────────────

let cards: Cards;
let answers: Answers;
let calc: Calc;              // the calculator's line and its history
let canvas: Doc;             // the canvas's document, kept while the calculator is shown
let past: History;
type HistoryHostSheet = ConstructorParameters<typeof History>[1]["sheet"];
let runTimer: ReturnType<typeof setTimeout> | 0 = 0, saveTimer: ReturnType<typeof setTimeout> | 0 = 0;

function changed() {
  clearTimeout(runTimer);
  runTimer = setTimeout(run, 120);
  clearTimeout(saveTimer);
  saveTimer = setTimeout(save, 300);
}

async function run() {
  if (!core.ready) {
    if (core.failed) for (const c of cards.document.cells) cards.show(c, { refused: `the core did not load (${core.failed})` });
    return;
  }
  cards.clearDerived();
  cards.pending(() => core.stop());
  await answers.run(cards.document, (cell, answer, above, from) => cards.show(cell, answer, above, from));
}

async function save() {
  const doc = cards.document;
  if (doc.preset === "calculator") {
    calc.line = doc.cells[0]?.kind === "math" ? doc.cells[0].written : empty();
    saveCalc(calc);
  } else try { localStorage.setItem(KEPT, JSON.stringify(toJSON(doc))); } catch {}
  try { history.replaceState(null, "", `#w=${await toLink(doc)}`); } catch {}
}

async function copyLink() {
  clearTimeout(saveTimer);
  await save();
  // said in place: the icon turns into a check for a moment, and a screen reader hears it
  const button = $("#copy-link")!, was = button.innerHTML;
  try {
    await navigator.clipboard.writeText(location.href);
    button.innerHTML = '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M5 12.5l4.5 4.5L19 7.5"/></svg>';
    $("#said")!.textContent = "Link copied";
    setTimeout(() => { button.innerHTML = was; $("#said")!.textContent = ""; }, 1500);
  } catch {
    toast("The address bar holds the link to this document");
  }
}

// the calculator's line as a document of one cell, its rows the line's own
const lineDoc = (): Doc => ({ preset: "calculator", cells: [math(calc.line)] });

function usePreset(preset: Preset) {
  const now = cards.document;
  if (now.preset === "canvas") canvas = now;
  else if (now.cells[0]?.kind === "math") calc.line = now.cells[0].written;
  for (const b of document.querySelectorAll<HTMLButtonElement>("[data-preset]")) b.setAttribute("aria-pressed", String(b.dataset.preset === preset));
  document.body.dataset.preset = preset;
  past.root.hidden = preset !== "calculator";
  if (preset !== now.preset) cards.setDoc(preset === "calculator" ? lineDoc() : canvas);
  else cards.build();
  changed();
}

// the line is done: what it asked goes on top of the history, shown open when `open`, and the
// line is empty again
function commitLine(open: boolean) {
  const cell = cards.document.cells[0];
  if (cards.document.preset !== "calculator" || cell?.kind !== "math") return;
  const sheet = cards.answerOf(cell)?.sheet as Sheet | undefined;
  let opened: Entry | null = null;
  if (sheet && !blank(cell.written)) {
    calc = commit(calc, { written: structuredClone(cell.written), result: sheet.result, strength: sheet.verdict?.strength ?? sheet.label,
      note: sheet.verdict?.notes[0] ?? "", ok: sheet.ok, t: Date.now() });
    opened = open ? calc.history[0] : null;
  }
  calc.line = empty();
  cards.setDoc(lineDoc());
  past.render(calc, opened);
  changed();
}

// a question put in the line, what was there committed first
function intoLine(written: Calc["line"]) {
  commitLine(false);
  calc.line = structuredClone(written);
  cards.setDoc(lineDoc());
  changed();
}

// ── the reference panel: each line of it a question to put in the line ─

// a question: its rows as the library writes them (`question_latex`), or a line of text -- on the
// calculator, which takes formulas only, that line as a formula -- put in the calculator's line,
// or added as a card at the end of the canvas
const addRows = (rows: string[] | undefined, text: string) => {
  const written = rows?.length ? { rows } : cards.document.preset === "calculator" && editor.ready ? { rows: [editor.toLatex(text)] } : { text };
  if (cards.document.preset === "calculator") intoLine(written);
  else cards.append([math(written)]);
};
const addText = (text: string) => addRows(undefined, text);
// a question shown as its rows, or as its text where the library wrote none
const shownAs = (rows: string[] | undefined, text: string, cls: string) =>
  rows?.length ? el("span", { class: cls }, ...rows.map((r) => tex(r, { cls: "math" }))) : el("span", { class: cls, text });

function panel(which: "ref" | null) {
  const box = $("#panel")!;
  for (const b of document.querySelectorAll<HTMLButtonElement>("[data-panel]")) b.setAttribute("aria-pressed", String(b.dataset.panel === which));
  if (!which) { box.hidden = true; return; }
  box.hidden = false;
  const body = el("div", { class: "panel-body" });
  {
    const r = data.reference;
    const filter = el("input", { class: "filter", type: "search", placeholder: "Find a function, an operator, a form", "aria-label": "Find in the reference" });
    const list = el("div", { class: "ref" });
    const entry = (name: string, e: { caption: string; example?: string }, domain?: string) =>
      el("div", { class: "entry", "data-find": `${name} ${e.caption}`.toLowerCase() },
        el("div", { class: "name mono", text: name }),
        el("div", {}, el("p", { class: "caption", text: e.caption + (domain ? `. Defined for ${domain}.` : "") }),
          e.example ? el("button", { type: "button", class: "ex", text: e.example, onclick: () => addText(e.example!) }) : null));
    if (!r) list.append(el("p", { class: "state", text: "The reference did not load." }));
    else list.append(
      el("h3", { text: "Questions" }), ...r.questions.map((q) =>
        el("div", { class: "entry", "data-find": `${q.form} ${q.caption}`.toLowerCase() },
          el("div", { class: "name" }, tex(q.form, { cls: "math" })),
          el("div", {}, el("p", { class: "caption", text: q.caption }),
            q.example ? el("button", { type: "button", class: "ex", onclick: () => addRows(q.answer?.question_latex, q.example!) },
              shownAs(q.answer?.question_latex, q.example, "asked")) : null))),
      el("h3", { text: "Functions" }), ...r.functions.map((f) => entry(`${f.name}(${f.arity === 2 ? "x, y" : "x"})`, f, f.domain)),
      el("h3", { text: "Operators" }), ...r.operators.map((o) => entry(o.written, o)),
      el("h3", { text: "Constants" }), ...r.constants.map((c) => entry(c.name, c)));
    filter.addEventListener("input", () => {
      const q = filter.value.trim().toLowerCase();
      for (const e of list.querySelectorAll<HTMLElement>(".entry")) e.hidden = Boolean(q) && !e.dataset.find!.includes(q);
    });
    body.append(filter, list);
  }
  box.replaceChildren(body);
}

// ── a hint for the word under the caret, from the reference ─────────────

// what the reference says of a function or a constant the word names, or begins alone
function hint(word: string): [string, string, string] | null {
  const r = data.reference;
  if (!r) return null;
  const words = [...r.functions.map((f) => ({ name: f.name, shown: `${f.name}(${f.arity === 2 ? "x, y" : "x"})`, e: f, domain: f.domain })),
    ...r.constants.map((c) => ({ name: c.name, shown: c.name, e: c, domain: "" }))];
  const begun = words.filter((w) => w.name.startsWith(word));
  const w = words.find((x) => x.name === word) ?? (begun.length === 1 ? begun[0] : undefined);
  return w ? [w.shown, w.e.caption + (w.domain ? `, for ${w.domain}` : ""), w.e.example ?? ""] : null;
}

// ── the keyboard: MathLive's own, in a window near the row being written ─

let activeField: HTMLElement | null = null;

function setKeys(open: boolean) {
  const keys = $("#keys")!;
  keys.hidden = !open;
  $("#open-keys")!.setAttribute("aria-expanded", String(open));
  const kb = (window as unknown as { mathVirtualKeyboard?: { show(): void; hide(): void } }).mathVirtualKeyboard;
  if (!editor.ready || !kb) return;
  if (open) {
    if (!activeField || !document.body.contains(activeField)) activeField = document.querySelector("math-field");
    activeField?.focus();
    kb.show();
  } else kb.hide();
  placeKeys();
}

// the keyboard is a window of its own over the page: it takes no room from the page, is moved by
// its title bar, and comes back where it was left (kept in this browser), inside the viewport
function placeKeys() {
  const keys = $("#keys")!;
  if (keys.hidden) return;
  let at: { x: number; y: number } | null = null;
  try { at = JSON.parse(localStorage.getItem("spatium.keys") ?? "null"); } catch { /* where CSS puts it */ }
  if (!at) { keys.style.left = keys.style.top = keys.style.bottom = keys.style.transform = ""; return; }
  const x = Math.min(Math.max(0, at.x), innerWidth - keys.offsetWidth), y = Math.min(Math.max(0, at.y), innerHeight - keys.offsetHeight);
  Object.assign(keys.style, { left: `${x}px`, top: `${y}px`, bottom: "auto", transform: "none" });
}

function dragKeys() {
  const keys = $("#keys")!, head = keys.querySelector<HTMLElement>(".keys-head")!;
  head.addEventListener("pointerdown", (e) => {
    if ((e.target as Element).closest("button")) return;
    const r = keys.getBoundingClientRect(), dx = e.clientX - r.left, dy = e.clientY - r.top;
    head.setPointerCapture(e.pointerId);
    const move = (m: PointerEvent) => {
      try { localStorage.setItem("spatium.keys", JSON.stringify({ x: m.clientX - dx, y: m.clientY - dy })); } catch { /* not kept */ }
      placeKeys();
    };
    head.addEventListener("pointermove", move);
    head.addEventListener("pointerup", () => head.removeEventListener("pointermove", move), { once: true });
  });
}

// ── settings of the page ────────────────────────────────────────────────

function applySettings() {
  const root = document.documentElement;
  if (settings.theme === "system") delete root.dataset.theme;
  else root.dataset.theme = settings.theme;
  root.style.setProperty("--scale", settings.scale);
  $("#open-keys")!.hidden = settings.keyboard !== "on" || !editor.ready;
  if (settings.keyboard !== "on") setKeys(false);
  for (const seg of document.querySelectorAll<HTMLElement>(".seg[data-setting]"))
    for (const b of seg.querySelectorAll<HTMLButtonElement>("button"))
      b.setAttribute("aria-pressed", String(b.value === settings[seg.dataset.setting as keyof PageSettings]));
}

function showCoreState() {
  const out = $("#core-state")!;
  if (core.failed) out.textContent = `Not loaded: ${core.failed}`;
  else if (!core.ready) out.textContent = "Loading";
  else out.textContent = `Core ${core.version || "loaded"}, ${Math.round(core.size / 1024)} KB, from the ${core.from}`;
}

function wire() {
  for (const b of document.querySelectorAll<HTMLButtonElement>("[data-preset]"))
    b.addEventListener("click", () => usePreset(b.dataset.preset as Preset));
  for (const b of document.querySelectorAll<HTMLButtonElement>("[data-panel]"))
    b.addEventListener("click", () => panel(b.getAttribute("aria-pressed") === "true" ? null : (b.dataset.panel as "ref")));
  $("#copy-link")!.addEventListener("click", copyLink);
  $("#open-keys")!.addEventListener("pointerdown", (e) => e.preventDefault());
  $("#open-keys")!.addEventListener("click", () => setKeys($("#keys")!.hidden === true));
  $("#close-keys")!.addEventListener("click", () => setKeys(false));
  // Esc closes the keyboard; with it closed, Esc twice within a second clears the line
  let escAt = 0;
  document.addEventListener("keydown", (e) => {
    if (e.key !== "Escape") return;
    if (!$("#keys")!.hidden) { setKeys(false); return; }
    if (cards.document.preset === "calculator" && Date.now() - escAt < 1000) { escAt = 0; commitLine(false); }
    else escAt = Date.now();
  });
  // forgetting asks once, in place: the button becomes Forget and Cancel
  const forgetRow = () => $("#forget-row")!.replaceChildren(el("button", { type: "button", class: "quiet warn", text: "Forget history and canvas",
    onclick: () => $("#forget-row")!.replaceChildren(
      el("button", { type: "button", class: "quiet warn", text: "Forget", onclick: () => {
        try { for (const k of ["spatium.calc", KEPT]) localStorage.removeItem(k); } catch { /* nothing kept */ }
        location.replace(location.pathname);
      } }),
      el("button", { type: "button", class: "quiet", text: "Cancel", onclick: forgetRow })) }));
  forgetRow();
  addEventListener("resize", placeKeys);
  dragKeys();
  const dialog = $<HTMLDialogElement>("#settings")!;
  $("#open-settings")!.addEventListener("click", () => { showCoreState(); dialog.showModal(); });
  dialog.addEventListener("click", (e) => { if (e.target === dialog) dialog.close(); });
  for (const seg of document.querySelectorAll<HTMLElement>(".seg[data-setting]"))
    seg.addEventListener("click", (e) => {
      const b = (e.target as Element).closest("button");
      if (!b) return;
      settings[seg.dataset.setting as keyof PageSettings] = b.value;
      saveSettings();
      applySettings();
      if (seg.dataset.setting === "formulas") cards.build();
      changed();
    });
  $("#check-update")!.addEventListener("click", async () => {
    $("#core-state")!.textContent = "Checking";
    try {
      if (await core.update()) { toast(`Core updated to ${core.version}`); changed(); }
      showCoreState();
    } catch (e) {
      $("#core-state")!.textContent = `Could not check: ${e instanceof Error ? e.message : e}`;
    }
  });
  addEventListener("hashchange", async () => {
    const doc = await fromAddress();
    if (doc) { cards.setDoc(doc); usePreset(doc.preset); }
  });
}

// ── start ───────────────────────────────────────────────────────────────

wire();
applySettings();
const coreLoaded = core.load();
await Promise.all([loadData(), loadEditor()]);
applySettings();
answers = new Answers(core, (latex) => editor.parse(latex));
const asked = await fromAddress();
calc = loadCalc();
canvas = fromBrowser() ?? starting();
canvas.preset = "canvas";
if (asked?.preset === "canvas") canvas = asked;
else if (asked?.cells[0]?.kind === "math") calc.line = asked.cells[0].written;
const doc = asked?.preset === "canvas" ? canvas : lineDoc();
past = new History(calc, {
  source: () => settings.formulas === "source",
  sheet: async (e) => {
    const one: Doc = { preset: "calculator", cells: [math(structuredClone(e.written))] };
    let got: Awaited<ReturnType<HistoryHostSheet>> = null;
    await answers.run(one, (cell, answer, _above, from) => { if (!from && cell === one.cells[0]) got = answer; });
    return got;
  },
  use: (e) => intoLine(e.written),
  changed: (c) => { calc = c; saveCalc(calc); },
  ask: (text) => addText(text),
  undo: undoStrip,
});
past.render(calc);
cards = new Cards(doc, {
  editor: editor.ready,
  source: () => settings.formulas === "source",
  changed,
  focused: (f) => { activeField = f; placeKeys(); },
  hint,
  curve: (cell) => (at) => answers.curve(cell, at),
  ask: (text) => addText(text),
  commit: () => commitLine(true),
});
$("#sheet")!.append(cards.root, past.root);
usePreset(doc.preset);
// the gallery's bar links here by name: the canvas, or the reference open
if (location.hash === "#canvas") usePreset("canvas");
if (location.hash === "#ref") panel("ref");
await coreLoaded;
showCoreState();
if (core.failed) toast(`The core did not load: ${core.failed}`, true);
changed();
