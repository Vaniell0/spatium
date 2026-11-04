// The page: one document of cells under the policy the bar picks -- the calculator or the canvas
// -- answered live by the core in a worker. It keeps the document in its address (#w=) and in
// the browser, and shows the JSON the core and the library's files hold; it composes nothing.

import { Answers, notation, type Parse } from "./answers.js";
import { Core, type State } from "./core.js";
import { type Doc, type Preset, fromJSON, fromLink, math, toJSON, toLink } from "./doc.js";
import { Cards } from "./view/cells.js";
import { $, el, math as tex } from "./view/dom.js";
import { kindOf } from "./view/answer.js";

// ── settings of the page (a cell's own are in the document) ─────────────

const DEFAULTS = { theme: "system", scale: "1", formulas: "rendered", keyboard: "on", cache: "on" };
type PageSettings = typeof DEFAULTS;
const settings: PageSettings = { ...DEFAULTS };
try { Object.assign(settings, JSON.parse(localStorage.getItem("spatium.settings") || "{}")); } catch {}
const saveSettings = () => { try { localStorage.setItem("spatium.settings", JSON.stringify(settings)); } catch {} };

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
    (window as unknown as { mathVirtualKeyboard: { container: Element } }).mathVirtualKeyboard.container = $("#keys-body")!;
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
interface Card { id: string; question?: string; question_latex?: string[]; typed?: boolean; route?: string; refused?: string; sheet?: { guarantee: string } }
const data: { gallery: Card[]; groups: { guarantee: string; title: string }[]; reference: Reference | null } = { gallery: [], groups: [], reference: null };

async function loadData() {
  const [g, r] = await Promise.all([
    fetch("gallery/gallery.json").then((x) => x.json()).catch(() => ({ cards: [] })),
    fetch("data/reference.json").then((x) => x.json()).catch(() => null),
  ]);
  data.gallery = g.cards;
  data.groups = g.groups ?? [];
  data.reference = r;
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
  await answers.run(cards.document, (cell, answer, above, from) => cards.show(cell, answer, above, from));
}

async function save() {
  const doc = cards.document;
  try { localStorage.setItem(KEPT, JSON.stringify(toJSON(doc))); } catch {}
  try { history.replaceState(null, "", `#w=${await toLink(doc)}`); } catch {}
}

async function copyLink() {
  clearTimeout(saveTimer);
  await save();
  try {
    await navigator.clipboard.writeText(location.href);
    toast("Link to this document copied");
  } catch {
    toast("The address bar holds the link to this document");
  }
}

function usePreset(preset: Preset) {
  cards.document.preset = preset;
  for (const b of document.querySelectorAll<HTMLButtonElement>("[data-preset]")) b.setAttribute("aria-pressed", String(b.dataset.preset === preset));
  document.body.dataset.preset = preset;
  cards.build();
  changed();
}

// ── panels: the reference and the gallery, each line of them a cell to add ─

// a question added as a card at the end: its rows as the library writes them (`question_latex`),
// or a line of text -- on the calculator, which takes formulas only, that line as a formula
const addRows = (rows: string[] | undefined, text: string) =>
  cards.append([math(rows?.length ? { rows } : cards.document.preset === "calculator" && editor.ready ? { rows: [editor.toLatex(text)] } : { text })]);
const addText = (text: string) => addRows(undefined, text);
// a question shown as its rows, or as its text where the library wrote none
const shownAs = (rows: string[] | undefined, text: string, cls: string) =>
  rows?.length ? el("span", { class: cls }, ...rows.map((r) => tex(r, { cls: "math" }))) : el("span", { class: cls, text });

function panel(which: "ref" | "gallery" | null) {
  const box = $("#panel")!;
  for (const b of document.querySelectorAll<HTMLButtonElement>("[data-panel]")) b.setAttribute("aria-pressed", String(b.dataset.panel === which));
  if (!which) { box.hidden = true; return; }
  box.hidden = false;
  const body = el("div", { class: "panel-body" });
  if (which === "gallery") {
    // by how strongly each answer holds, in the gallery's own order of groups
    for (const g of data.groups) {
      const picks = data.gallery.filter((c) => c.question && c.typed !== false && c.sheet?.guarantee === g.guarantee);
      if (!picks.length) continue;
      body.append(el("h3", { text: g.title }));
      for (const c of picks) {
        const [ok, label] = kindOf(c as unknown as Record<string, unknown>);
        body.append(el("button", { type: "button", class: "pick", onclick: () => addRows(c.question_latex, c.question!) },
          shownAs(c.question_latex, c.question!, "q asked"), el("span", { class: ok ? "kind ok" : "kind", text: label })));
      }
    }
  } else {
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

// the panel stays at the bottom edge; the page gets room below its last row, and the row being
// written is scrolled above the panel when the panel would cover it
function placeKeys() {
  const keys = $("#keys")!;
  document.body.style.paddingBottom = keys.hidden ? "" : `${keys.offsetHeight + 32}px`;
  if (keys.hidden || !activeField) return;
  const r = activeField.getBoundingClientRect(), top = keys.getBoundingClientRect().top;
  if (r.bottom + 16 > top) scrollBy({ top: r.bottom + 16 - top, behavior: "smooth" });
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
    b.addEventListener("click", () => panel(b.getAttribute("aria-pressed") === "true" ? null : (b.dataset.panel as "ref" | "gallery")));
  $("#copy-link")!.addEventListener("click", copyLink);
  $("#open-keys")!.addEventListener("pointerdown", (e) => e.preventDefault());
  $("#open-keys")!.addEventListener("click", () => setKeys($("#keys")!.hidden === true));
  $("#close-keys")!.addEventListener("click", () => setKeys(false));
  document.addEventListener("keydown", (e) => { if (e.key === "Escape") setKeys(false); });
  addEventListener("resize", placeKeys);
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
const doc = (await fromAddress()) ?? fromBrowser() ?? starting();
cards = new Cards(doc, {
  editor: editor.ready,
  source: () => settings.formulas === "source",
  changed,
  focused: (f) => { activeField = f; placeKeys(); },
  hint,
  curve: (cell) => (at) => answers.curve(cell, at),
});
$("#sheet")!.append(cards.root);
usePreset(doc.preset);
await coreLoaded;
showCoreState();
if (core.failed) toast(`The core did not load: ${core.failed}`, true);
changed();
