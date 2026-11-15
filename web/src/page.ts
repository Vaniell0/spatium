// The page: the calculator -- one line and the history of what was asked -- or the canvas, a
// document of cells, as the bar picks, answered live by the core in a worker. It keeps what is
// written in its address (#w=) and in the browser, and shows the JSON the core and the library's
// files hold; it composes nothing.

import { Answers, notation, type Parse } from "./answers.js";
import { type Calc, type Entry, blank, commit, empty, load as loadCalc, save as saveCalc } from "./calc.js";
import { Core, type State } from "./core.js";
import { type Doc, type Preset, fromJSON, fromLink, math, toJSON, toLink } from "./doc.js";
import { Cards } from "./view/cells.js";
import { $, el, icon, math as tex } from "./view/dom.js";
import type { Sheet } from "./view/answer.js";
import { History } from "./view/history.js";
import { Keys } from "./view/keys.js";

// ── settings of the page (a cell's own are in the document) ─────────────

const DEFAULTS = { theme: "system", scale: "1", formulas: "rendered", keyboard: "auto", cache: "on" };
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

// the core's state is the dot's colour; its words are in the settings dialog
function dot(state: State) {
  $("#core-dot")!.className = `core-dot ${state === "loading" ? "" : state}`;
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
    editor.parse = (latex: string) => ls.parse(notation(latex));
    editor.toLatex = ml.convertAsciiMathToLatex;
    editor.ready = true;
  } catch (e) {
    toast(`The formula editor did not load (${e instanceof Error ? e.message : e}); cells are written as text.`, true);
  }
}

// the reference, as site_reference writes it from the library's tables: every example answered
// by the library (`answer`), every formula the library's or read back by the editor's parser in CI
interface Answered { route?: string; refused?: string; question_latex?: string[]; label?: string }
interface Reference {
  writing: { typed: string; reads: string; caption: string; answer: Answered }[];
  questions: { form: string; caption: string; example: string; answer: Answered }[];
  functions: { form?: string; reads?: string; name: string; arity?: number; caption: string; domain: string; rows?: string[]; example: string; answer: Answered }[];
  operators: { written: string; caption: string; example: string; answer: Answered }[];
  constants: { name: string; caption: string }[];
  language: { form: string; caption: string; example: string; latex: string[] }[];
  shapes: { name: string; arguments: string; caption: string }[];
  told: { caption: string; example: string; told: string }[];
  guarantees: { name: string; title: string; promise: string }[];
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
  const button = $("#copy-link")!;
  try {
    await navigator.clipboard.writeText(location.href);
    button.replaceChildren(icon("check"));
    $("#said")!.textContent = "Link copied";
    setTimeout(() => { button.replaceChildren(icon("link")); $("#said")!.textContent = ""; }, 1500);
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
  showExamples();
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
  showExamples();
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

// The reference: a drawer of sections, as a page of docs is -- how to write, the questions, the
// functions, the operators, the constants, the canvas's language, what an answer promises -- each
// entry in both notations, the page's own first (LaTeX on the calculator, the language on the
// canvas), and its example a button that puts it in the line or on the canvas.
const REF = "spatium.ref";   // the drawer's scroll, for the session
const sections = ["Writing", "Questions", "Functions", "Operators", "Constants", "Language", "Answers"] as const;

function panel(which: "ref" | null) {
  const box = $("#panel")!;
  for (const b of document.querySelectorAll<HTMLButtonElement>("[data-panel]")) b.setAttribute("aria-pressed", String(b.dataset.panel === which));
  $("#scrim")!.hidden = !which;
  if (!which) { box.hidden = true; return; }
  box.hidden = false;
  // under the bar, which is two rows on a phone
  const top = `${Math.round($(".bar")!.getBoundingClientRect().bottom)}px`;
  box.style.top = top, $("#scrim")!.style.top = top;
  if (box.dataset.preset === cards.document.preset && box.childElementCount) return;   // drawn already, for this page
  box.dataset.preset = cards.document.preset;
  const r = data.reference;
  const dslFirst = cards.document.preset === "canvas";
  const close = el("button", { type: "button", class: "icon", "aria-label": "Close the reference", onclick: () => panel(null) }, icon("x"));
  const filter = el("input", { class: "filter", type: "search", placeholder: "Find in the reference", "aria-label": "Find in the reference" });
  const index = el("nav", { class: "ref-index", "aria-label": "Sections" });
  const doc = el("div", { class: "ref" });
  const rows = (latex: string[]) => el("span", { class: "rows" }, ...latex.map((l) => tex(l, { cls: "math" })));
  const mono = (text: string) => el("code", { class: "mono", text });
  // the notations of one thing, the page's own first
  const both = (latex: HTMLElement | null, dsl: HTMLElement | null) =>
    el("div", { class: "forms" }, ...(dslFirst ? [dsl, latex] : [latex, dsl]));
  const example = (latex: string[] | undefined, text: string, label?: string) =>
    el("button", { type: "button", class: "ex", "aria-label": `Ask ${text}`, onclick: () => addRows(latex, text) },
      dslFirst || !latex?.length ? mono(text) : rows(latex), label ? el("span", { class: "said", text: label }) : null);
  const entry = (find: string, head: HTMLElement, ...rest: (HTMLElement | null)[]) =>
    el("div", { class: "entry", "data-find": find.toLowerCase() }, head, el("div", { class: "what" }, ...rest));
  const caption = (text: string) => el("p", { class: "caption", text });
  const part = (title: (typeof sections)[number], intro: string | null, ...entries: HTMLElement[]) =>
    el("section", { class: "ref-sec", id: `ref-${title.toLowerCase()}` }, el("h3", { text: title }), intro ? caption(intro) : null, ...entries);
  if (!r) doc.append(el("p", { class: "state", text: "The reference did not load." }));
  else doc.append(
    part("Writing", "What is typed in the formula editor, the formula it becomes, and the line of the language it reads as.",
      ...r.writing.map((w) => entry(`${w.typed} ${w.reads} ${w.caption}`,
        el("kbd", { text: w.typed }), example(w.answer.question_latex, w.reads), dslFirst ? null : mono(w.reads), caption(w.caption)))),
    part("Questions", null, ...r.questions.map((q) => entry(`${q.form} ${q.caption} ${q.example}`,
      el("div", { class: "name" }, tex(q.form, { cls: "math" })), caption(q.caption),
      example(q.answer.question_latex, q.example, q.answer.label)))),
    part("Functions", null, ...r.functions.map((f) => entry(`${f.name} ${f.form ?? ""} ${f.caption} ${f.example}`,
      both(f.form ? tex(f.form, { cls: "math" }) : null, mono(f.name ? `${f.name}(${f.arity === 2 ? "x, y" : "x"})` : f.reads ?? "")),
      caption(`${f.caption}; defined for ${f.domain}`), example(f.rows ?? f.answer.question_latex, f.example)))),
    part("Operators", "Tightest last.", ...r.operators.map((o) => entry(`${o.written} ${o.caption}`,
      mono(o.written), caption(o.caption), example(o.answer.question_latex, o.example)))),
    part("Constants", null, ...r.constants.map((c) => entry(`${c.name} ${c.caption}`, mono(c.name), caption(c.caption)))),
    part("Language", "The canvas's text: a line a definition, reading the names the lines above it define.",
      ...r.language.map((l) => entry(`${l.form} ${l.caption} ${l.example}`,
        mono(l.form), caption(l.caption), example(l.latex.length ? l.latex : undefined, l.example))),
      entry("shape shapes " + r.shapes.map((s) => s.name).join(" "), el("span", { class: "name", text: "Shapes" }),
        el("dl", { class: "shapes" }, ...r.shapes.flatMap((s) => [el("dt", {}, mono(`${s.name}(${s.arguments})`)), el("dd", { text: s.caption })]))),
      ...r.told.map((t) => entry(`${t.caption} ${t.told}`, el("pre", { class: "mono", text: t.example }),
        caption(t.caption), el("p", { class: "told", text: t.told })))),
    part("Answers", "How strongly an answer holds, strongest first.", ...r.guarantees.map((g) => entry(`${g.title} ${g.promise}`,
      el("span", { class: "name", text: g.title }), caption(g.promise)))));
  for (const title of sections)
    index.append(el("button", { type: "button", class: "chip", text: title,
      onclick: () => doc.querySelector(`#ref-${title.toLowerCase()}`)?.scrollIntoView({ block: "start" }) }));
  filter.addEventListener("input", () => {
    const q = filter.value.trim().toLowerCase();
    for (const e of doc.querySelectorAll<HTMLElement>(".entry")) e.hidden = Boolean(q) && !e.dataset.find!.includes(q);
    for (const sec of doc.querySelectorAll<HTMLElement>(".ref-sec")) sec.hidden = !sec.querySelector(".entry:not([hidden])");
  });
  doc.addEventListener("scroll", () => { try { sessionStorage.setItem(REF, String(doc.scrollTop)); } catch { /* not kept */ } }, { passive: true });
  box.replaceChildren(el("div", { class: "drawer-head" }, filter, close), index, doc);
  try { doc.scrollTop = Number(sessionStorage.getItem(REF) ?? 0); } catch { /* from the top */ }
}

// ── a hint for the word under the caret, from the reference ─────────────

// what the reference says of a function or a constant the word names, or begins alone
function hint(word: string): [string, string, string] | null {
  const r = data.reference;
  if (!r) return null;
  const words = [...r.functions.filter((f) => f.name).map((f) => ({ name: f.name, shown: `${f.name}(${f.arity === 2 ? "x, y" : "x"})`, caption: f.caption, domain: f.domain, example: f.example })),
    ...r.constants.map((c) => ({ name: c.name, shown: c.name, caption: c.caption, domain: "", example: "" }))];
  const begun = words.filter((w) => w.name.startsWith(word));
  const w = words.find((x) => x.name === word) ?? (begun.length === 1 ? begun[0] : undefined);
  return w ? [w.shown, w.caption + (w.domain ? `, for ${w.domain}` : ""), w.example] : null;
}

// ── the keyboard: the page's own panel over the line ────────────────────

let activeField: HTMLElement | null = null;
let keys: Keys;
const coarse = () => matchMedia("(pointer: coarse)").matches;

// the row being written, or the line's first row
const field = () => {
  if (!activeField || !document.body.contains(activeField)) activeField = document.querySelector("math-field");
  return activeField as (HTMLElement & { executeCommand(c: unknown): void }) | null;
};

// the panel floats over the page and adds nothing to it: the live cell's ask row turns sticky so
// the line stays above the keys. On touch the hidden textarea takes no OS keyboard while it is open.
function setKeys(open: boolean) {
  keys.root.hidden = !open;
  document.documentElement.style.setProperty("--bar-bottom", `${Math.round($(".bar")!.getBoundingClientRect().bottom)}px`);
  document.body.classList.toggle("keys-open", open);
  $("#open-keys")!.setAttribute("aria-expanded", String(open));
  const mf = field();
  mf?.shadowRoot?.querySelector("textarea, [contenteditable]")?.setAttribute("inputmode", open && coarse() ? "none" : "");
  if (open) mf?.focus();
}

// a row took the focus: the keyboard opens as the setting says -- always, or on touch
function rowFocused(f: HTMLElement) {
  activeField = f;
  if (settings.keyboard === "on" || (settings.keyboard === "auto" && coarse())) setKeys(true);
}

// ── settings of the page ────────────────────────────────────────────────

function applySettings() {
  const root = document.documentElement;
  if (settings.theme === "system") delete root.dataset.theme;
  else root.dataset.theme = settings.theme;
  root.style.setProperty("--scale", settings.scale);
  $("#open-keys")!.hidden = settings.keyboard === "off" || !editor.ready;
  if (settings.keyboard === "off") setKeys(false);
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
  $("#scrim")!.addEventListener("click", () => panel(null));
  $("#copy-link")!.addEventListener("click", copyLink);
  $("#open-keys")!.addEventListener("pointerdown", (e) => e.preventDefault());
  $("#open-keys")!.addEventListener("click", () => setKeys(keys.root.hidden));
  // Esc closes the keyboard; with it closed, Esc twice within a second clears the line
  let escAt = 0;
  document.addEventListener("keydown", (e) => {
    if (e.key !== "Escape") return;
    if (!keys.root.hidden) { setKeys(false); return; }
    if (!$("#panel")!.hidden) { panel(null); return; }
    if (cards.document.preset === "calculator" && Date.now() - escAt < 1000) {
      escAt = 0;
      const was = structuredClone(calc);
      if (blank(calc.line)) return;
      commitLine(false);
      undoStrip("Line cleared", () => { calc = was; cards.setDoc(lineDoc()); past.render(calc); changed(); });
    } else escAt = Date.now();
  });
  // Ctrl+Z undoes what the strip offers to; inside a formula that has something in it, the editor's
  // own undo goes first
  document.addEventListener("keydown", (e) => {
    if (!(e.ctrlKey || e.metaKey) || e.shiftKey || e.key.toLowerCase() !== "z" || !undoing) return;
    const field = document.activeElement?.closest("math-field") as (HTMLElement & { value: string }) | null;
    if (field?.value || (document.activeElement instanceof HTMLTextAreaElement && document.activeElement.value)) return;
    e.preventDefault();
    undoing();
    undoing = null;
    $("#undo")!.hidden = true;
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

// ── examples: three of the gallery's questions, shown while the history is empty ──

const examples = el("section", { class: "examples", "aria-label": "Examples", hidden: true });

async function loadExamples() {
  interface Card { question_latex?: string[]; sheet?: { guarantee: string; result: string; verdict?: { strength: string } } }
  const g = await fetch("gallery/gallery.json").then((x) => x.json() as Promise<{ groups: { guarantee: string }[]; cards: Card[] }>).catch(() => null);
  // one of each of the first three guarantees, the shortest question of each: the range, not the depth
  const short = (c: Card) => c.question_latex!.join("").length;
  const picks = (g?.groups ?? []).map((grp) => (g?.cards ?? []).filter((c) => c.question_latex?.length && c.sheet?.result && c.sheet.guarantee === grp.guarantee)
    .sort((a, b) => short(a) - short(b))[0]).filter((c): c is Card => Boolean(c)).slice(0, 3);
  if (!picks.length) return;
  examples.replaceChildren(el("div", { class: "history-head" }, el("h2", { text: "Examples" })), ...picks.map((c) => {
    const rows = c.question_latex!;
    return el("article", { class: "cell hcell" },
      el("div", { class: "ask" }, el("span", { class: "prompt", "aria-hidden": "true", text: "›" }), el("span", { class: "q" }, ...rows.map((r) => tex(r, { cls: "math" })))),
      el("button", { type: "button", class: "result-row", "aria-label": "Put in the line", onclick: () => intoLine({ rows: [...rows] }) },
        tex(c.sheet!.result, { cls: "math res" }), c.sheet!.verdict ? el("span", { class: "strength ok", text: c.sheet!.verdict.strength }) : null));
  }));
  showExamples();
}

const showExamples = () => { if (cards) examples.hidden = cards.document.preset !== "calculator" || calc.history.length > 0 || !examples.childElementCount; };

// ── start ───────────────────────────────────────────────────────────────

keys = new Keys({ field, commit: () => commitLine(true), close: () => setKeys(false) });
document.body.append(keys.root);
wire();
applySettings();
const coreLoaded = core.load();
await Promise.all([loadData(), loadEditor(), loadExamples()]);
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
  changed: (c) => { calc = c; saveCalc(calc); showExamples(); },
  ask: (text, rows) => addRows(rows, text),
  undo: undoStrip,
});
past.render(calc);
cards = new Cards(doc, {
  editor: editor.ready,
  source: () => settings.formulas === "source",
  changed,
  focused: rowFocused,
  hint,
  curve: (cell) => (at) => answers.curve(cell, at),
  ask: (text, rows) => addRows(rows, text),
  commit: () => commitLine(true),
});
$("#sheet")!.append(cards.root, examples, past.root);
showExamples();
usePreset(doc.preset);
// the gallery's bar links here by name: the canvas, or the reference open
if (location.hash === "#canvas") usePreset("canvas");
if (location.hash === "#ref") panel("ref");
await coreLoaded;
showCoreState();
if (core.failed) toast(`The core did not load: ${core.failed}`, true);
changed();
