// The page: windows on a workspace, each a view of what the library answers. Nothing here
// computes an answer or writes a formula. A question is typed as rows of formulas in a LaTeX
// editor (MathLive); the editor's parser (the Compute Engine's latex-syntax) writes each row as
// MathJSON, and the core (wasm/spatium_core.wasm, the library built to WebAssembly, run in
// core-worker.js) reads the rows as a document and answers with its steps and its certificate as
// JSON. A question may also be typed as text in the language. The gallery and the reference are
// files the library wrote. The page chooses where each part goes.

const CORE_JS = "wasm/spatium_core.js";
const CORE_WASM = "wasm/spatium_core.wasm";
const CACHE = "spatium-core";

// ── settings ──────────────────────────────────────────────────────────

const DEFAULTS = { theme: "system", scale: "1", digits: "12", formulas: "rendered", keyboard: "on", cache: "on" };
const settings = { ...DEFAULTS };
try { Object.assign(settings, JSON.parse(localStorage.getItem("spatium.settings") || "{}")); } catch {}

function saveSettings() {
  try { localStorage.setItem("spatium.settings", JSON.stringify(settings)); } catch {}
}

function applySettings() {
  const root = document.documentElement;
  if (settings.theme === "system") delete root.dataset.theme;
  else root.dataset.theme = settings.theme;
  root.style.setProperty("--scale", settings.scale);
  $("#open-keys").hidden = settings.keyboard !== "on" || !editor.ready;
  if (settings.keyboard !== "on") setKeys(false);
  for (const seg of document.querySelectorAll(".seg[data-setting]"))
    for (const b of seg.querySelectorAll("button"))
      b.setAttribute("aria-pressed", String(b.value === settings[seg.dataset.setting]));
}

// ── small helpers ─────────────────────────────────────────────────────

const $ = (s, root = document) => root.querySelector(s);

function el(tag, props = {}, ...children) {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(props)) {
    if (k === "class") e.className = v;
    else if (k === "text") e.textContent = v;
    else if (k.startsWith("on")) e.addEventListener(k.slice(2), v);
    else e.setAttribute(k, v);
  }
  for (const c of children) if (c != null && c !== false) e.append(c);
  return e;
}

// a formula as the library wrote it, rendered, or its LaTeX source when that is asked for
function math(latex, { display = false, cls = "math" } = {}) {
  const box = el("div", { class: cls });
  if (settings.formulas === "source" || !window.katex) {
    box.append(el("code", { class: "source", text: latex }));
    return box;
  }
  window.katex.render(latex, box, { throwOnError: false, displayMode: display });
  return box;
}

function number(v) {
  if (typeof v !== "number" || !Number.isFinite(v)) return String(v);
  return String(Number(v.toPrecision(Number(settings.digits))));
}

// an error bound written as the power of ten above it, as the library gives it
function decade(e) {
  if (!e) return el("span", { text: "exact, to the last digit" });
  return math(`10^{${Math.round(Math.log10(e))}}`);
}

// 10^-12 as a person reads it in a line of text
function decadeText(e) {
  const sup = { "-": "⁻", 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹" };
  return "10" + [...String(Math.round(Math.log10(e)))].map((d) => sup[d]).join("");
}

let toastTimer;
function toast(text, bad = false) {
  const t = $("#toast");
  clearTimeout(toastTimer);
  t.textContent = text;
  t.hidden = false;
  t.classList.toggle("bad", bad);
  t.classList.remove("gone");
  if (!bad) toastTimer = setTimeout(() => {
    t.classList.add("gone");
    toastTimer = setTimeout(() => { t.hidden = true; }, 500);
  }, 3000);
}

// ── the core, in a worker ─────────────────────────────────────────────

const core = { ready: false, failed: null, version: "", size: 0, from: "", glue: "", wasm: null, worker: null };

async function fetchBytes(path, fresh) {
  const useCache = settings.cache === "on" && "caches" in window;
  if (useCache && !fresh) {
    try {
      const hit = await (await caches.open(CACHE)).match(path);
      if (hit) return { bytes: await hit.arrayBuffer(), from: "cache" };
    } catch {}
  }
  const res = await fetch(path, { cache: fresh ? "no-store" : "default" });
  if (!res.ok) throw new Error(`${path} answered ${res.status}`);
  return { bytes: await res.arrayBuffer(), from: "network" };
}

async function storeBytes(path, bytes, type) {
  if (settings.cache !== "on" || !("caches" in window)) return;
  try { await (await caches.open(CACHE)).put(path, new Response(bytes.slice(0), { headers: { "content-type": type } })); } catch {}
}

async function digest(bytes) {
  if (!crypto.subtle) return "";
  const h = new Uint8Array(await crypto.subtle.digest("SHA-256", bytes));
  return [...h.slice(0, 4)].map((b) => b.toString(16).padStart(2, "0")).join("");
}

function dot(state, title) {
  const d = $("#core-dot");
  d.className = `core-dot ${state}`;
  d.title = title;
}

// the glue and the module are kept and replaced together, so a cached module never meets a newer glue
async function useCore(glue, wasm) {
  if (core.glue) URL.revokeObjectURL(core.glue);
  core.glue = URL.createObjectURL(new Blob([glue], { type: "text/javascript" }));
  core.wasm = wasm;
  core.version = await digest(wasm);
  core.size = wasm.byteLength;
  await spawn();
}

// a worker with the core in it; a new one replaces one that is ended mid-question
function spawn() {
  if (core.worker) core.worker.terminate();
  core.ready = false;
  const w = new Worker("core-worker.js");
  core.worker = w;
  return new Promise((done) => {
    w.onmessage = (e) => {
      const m = e.data;
      if (m.type === "ready") { core.ready = true; core.failed = null; dot("ready", "Core ready"); done(); pump(); }
      else if (m.type === "failed") { core.failed = m.why; dot("failed", `Core failed: ${m.why}`); done(); }
      else if (m.type === "answer") answered(m);
    };
    w.onerror = (e) => { core.failed = e.message || "the worker did not start"; dot("failed", core.failed); done(); };
    w.postMessage({ type: "start", glue: core.glue, wasm: core.wasm.slice(0) });
  });
}

async function loadCore() {
  dot("", "Core loading");
  try {
    const [glue, wasm] = await Promise.all([fetchBytes(CORE_JS), fetchBytes(CORE_WASM)]);
    await useCore(glue.bytes, wasm.bytes);
    core.from = wasm.from;
    if (core.failed) throw new Error(core.failed);
    if (wasm.from === "network") {
      await storeBytes(CORE_JS, glue.bytes, "text/javascript");
      await storeBytes(CORE_WASM, wasm.bytes, "application/wasm");
    }
  } catch (e) {
    core.failed = e.message || String(e);
    dot("failed", `Core did not load: ${core.failed}`);
    toast(`The core did not load: ${core.failed}`, true);
    for (const item of work.items) if (item.kind === "ask") recompute(item);
  }
  showCoreState();
}

async function checkUpdate() {
  const out = $("#core-state");
  out.textContent = "Checking";
  try {
    const [glue, wasm] = await Promise.all([fetchBytes(CORE_JS, true), fetchBytes(CORE_WASM, true)]);
    const version = await digest(wasm.bytes);
    if (core.ready && version && version === core.version) {
      out.textContent = `Up to date, core ${version}`;
      return;
    }
    await storeBytes(CORE_JS, glue.bytes, "text/javascript");
    await storeBytes(CORE_WASM, wasm.bytes, "application/wasm");
    await useCore(glue.bytes, wasm.bytes);
    core.from = "network";
    showCoreState();
    toast(`Core updated to ${core.version}`);
    for (const item of work.items) if (item.kind === "ask") recompute(item);
  } catch (e) {
    out.textContent = `Could not check: ${e.message || e}`;
  }
}

function showCoreState() {
  const out = $("#core-state");
  if (core.failed) out.textContent = `Not loaded: ${core.failed}`;
  else if (!core.ready) out.textContent = "Loading";
  else out.textContent = `Core ${core.version || "loaded"}, ${Math.round(core.size / 1024)} KB, from the ${core.from}`;
}

// One question at a time goes to the worker; each window keeps only its newest. A question made
// stale while the worker is on it -- a synthesis typed past -- is given a moment, then the
// worker is ended and a new one started: that is how a synthesis is cancelled.
const STALE_MS = 400;
const waiting = new Map();   // window id -> { kind, payload, resolve }
let inflight = null;         // { item, id, resolve, timer }
let asked = 0;

function ask(item, kind, payload) {
  return new Promise((resolve) => {
    const old = waiting.get(item);
    if (old) old.resolve(null);
    waiting.set(item, { kind, payload, resolve });
    if (inflight && inflight.item === item && !inflight.timer)
      inflight.timer = setTimeout(() => {
        if (!inflight || inflight.item !== item) return;
        inflight.resolve(null);
        inflight = null;
        spawn();
      }, STALE_MS);
    pump();
  });
}

function pump() {
  if (!core.ready || inflight || waiting.size === 0) return;
  const [item, q] = waiting.entries().next().value;
  waiting.delete(item);
  inflight = { item, id: ++asked, resolve: q.resolve, timer: 0 };
  dot("busy", "Core working");
  core.worker.postMessage({ type: "ask", id: inflight.id, kind: q.kind, payload: q.payload });
}

function answered(m) {
  if (!inflight || m.id !== inflight.id) return;
  clearTimeout(inflight.timer);
  const done = inflight.resolve;
  inflight = null;
  dot("ready", "Core ready");
  done(m.answer);
  pump();
}

// ── data the library wrote ────────────────────────────────────────────

const data = { gallery: [], reference: null };

async function loadData() {
  const [g, r] = await Promise.all([
    fetch("gallery/gallery.json").then((x) => x.json()).catch(() => ({ cards: [] })),
    fetch("data/reference.json").then((x) => x.json()).catch(() => null),
  ]);
  data.gallery = g.cards;
  data.reference = r;
}

// ── the editor: MathLive for the rows, latex-syntax for their MathJSON ─

const editor = { ready: false, parse: null };

async function loadEditor() {
  try {
    const [ml, ls] = await Promise.all([
      import("./vendor/mathlive/mathlive.min.mjs"),
      import("./vendor/latex-syntax/latex-syntax.min.js"),
    ]);
    ml.MathfieldElement.fontsDirectory = new URL("vendor/katex/fonts/", location.href).href;
    ml.MathfieldElement.soundsDirectory = null;
    ml.MathfieldElement.computeEngine = null;
    window.mathVirtualKeyboard.container = $("#keys-body");
    editor.parse = ls.parse;
    editor.ready = true;
  } catch (e) {
    toast(`The formula editor did not load (${e.message || e}); questions are typed as text.`, true);
  }
}

// ── an answer ─────────────────────────────────────────────────────────

const HOW = {
  evaluate: "Evaluated as written, in double precision, with nothing rewritten.",
  "closed form": "Solved in closed form. The coefficients are read from the equation's own values and the degree is checked at more points; the roots come from the quadratic formula in the form that keeps the small root, and Vieta's formulas check them by a second route.",
  synthesis: "Synthesised. The library computes the question's value at each parameter by its own solvers; the answer is the shortest expression in the parameter that matches it, searched bottom up, or a Chebyshev polynomial where none of a few terms does, and it is checked against those values on a dense grid.",
};
const COMPUTED = { wasm: "in this browser", host: "on the host, when the gallery was written" };

function kindOf(a) {
  if (a.route === "scene") return ["", "A scene, as LaTeX"];
  if (a.refused) return ["", a.route ? `Not answered by ${a.route}` : "Not answered"];
  if (a.route === "evaluate") return ["", "Evaluated"];
  if (a.route === "expression") return ["", `A function of ${a.reads}, as written`];
  if (a.route === "closed form") return a.check && a.check.holds ? ["ok", "Closed form, checked by Vieta's formulas"] : ["", "Closed form"];
  const c = a.certificate;
  if (!c) return ["", a.route || ""];
  const e = c.error_below;
  return ["ok", e ? `Certified, ${c.norm === "largest relative error" ? "relative " : ""}error below ${decadeText(e)}` : "Certified, exact on the grid"];
}

function facts(rows) {
  const dl = el("dl", { class: "facts" });
  for (const [k, v] of rows) dl.append(el("dt", { text: k }), el("dd", {}, v));
  return dl;
}

// the step a closed-form answer ends in: its roots, or the root of a linear equation
function finalStep(steps) {
  const named = steps.find((s) => s.label === "roots" || s.label === "root");
  if (named) return named;
  const solving = steps.filter((s) => !s.label.startsWith("check"));
  return solving[solving.length - 1];
}

function how(route) {
  return HOW[route] ? el("details", { class: "more" }, el("summary", { text: "How it was found" }), el("p", { text: HOW[route] })) : null;
}

// a refusal: what was not answered and why, with the place in a line of text when the reason names a column
function refusal(a, typed) {
  const box = el("div", { class: "refusal", role: "status" },
    el("p", { class: "what", text: kindOf(a)[1] }),
    el("p", { class: "why", text: a.refused }));
  const col = Number((/column (\d+)/.exec(a.refused) || [])[1] || 0);
  if (typed && col > 0 && !typed.includes("\n")) {
    const at = Math.min(col - 1, typed.length);
    box.append(el("div", { class: "caret" }, typed.slice(0, at), el("mark", { text: typed.slice(at, at + 1) || " " }), typed.slice(at + 1)));
  }
  return box;
}

function answerView(a, typed) {
  const out = document.createDocumentFragment();
  if (a.refused) {
    out.append(refusal(a, typed));
    return out;
  }
  const [ok, label] = kindOf(a);
  out.append(el("p", { class: ok ? "kind ok" : "kind", text: label }));
  if (a.route === "scene") {
    out.append(math(a.latex, { display: true, cls: "big" }));
    return out;
  }
  // a value of one parameter: its own expression, as written
  if (a.route === "expression") {
    out.append(math(a.latex, { display: true, cls: "big accent" }));
    out.append(facts([["Reads", el("span", { class: "mono", text: a.reads })]]));
    return out;
  }
  if (a.route === "evaluate") {
    out.append(math(a.latex, { display: true, cls: "big" }));
    out.append(facts([["Value", el("span", { class: "mono", text: number(a.value) })]]), how(a.route));
    return out;
  }
  if (a.route === "closed form") {
    const last = finalStep(a.steps);
    out.append(math(last.latex, { display: true, cls: "big accent" }));
    const chain = el("ol", { class: "chain", "aria-label": "Steps" });
    a.steps.forEach((s, k) => {
      const cls = s === last ? "pay" : s.label.startsWith("check") ? "check" : "";
      chain.append(el("li", { class: cls }, el("span", { class: "n", text: String(k + 1) }), el("span", { class: "lab", text: s.label }),
        math(s.latex, { display: true })));
    });
    out.append(chain);
    if (a.check) out.append(facts([
      ["Check", el("span", { class: a.check.holds ? "ok" : "", text: `${a.check.by} ${a.check.holds ? "hold" : "do not hold"}` })],
      [a.check.error_below ? "Relative error below" : "Error", decade(a.check.error_below)],
    ]));
    out.append(how(a.route));
    return out;
  }
  if (a.route === "synthesis") {
    const pieces = a.answer || [];
    if (pieces.length === 1 && pieces[0].form !== "chebyshev") out.append(math(pieces[0].latex, { display: true, cls: "big accent" }));
    else out.append(el("p", { class: "big", text: pieces.length === 1 ? `A Chebyshev series of degree ${pieces[0].degree}` : `${a.pieces} pieces` }));
    const c = a.certificate;
    if (c) out.append(facts([
      ["Holds", el("span", { class: "ok", text: `${c.strength}, on ${c.points} points` })],
      [c.norm === "largest relative error" ? "Relative error below" : "Absolute error below", decade(c.error_below)],
      ["Scalar", el("span", { text: c.scalar })],
      ["Computed", el("span", { text: COMPUTED[c.computed] || c.computed })],
    ]));
    if (pieces.length > 1 || (pieces[0] && pieces[0].form === "chebyshev")) {
      const list = el("div", { class: "pieces" });
      for (const p of pieces) list.append(pieceView(p, pieces.length > 1));
      out.append(list);
    }
    out.append(how(a.route));
    return out;
  }
  return out;
}

function pieceView(p, many) {
  const box = el("div", { class: "piece" });
  box.append(el("p", { class: "label", text: `on [${number(p.lo)}, ${number(p.hi)}]` }));
  if (many && p.form !== "chebyshev") box.append(math(p.latex));
  if (p.form === "chebyshev" && p.coefficients) {
    const list = el("ol", { class: "coeffs" });
    p.coefficients.forEach((c, k) => list.append(el("li", {}, el("span", { text: `c${k}` }), number(c))));
    box.append(el("details", {}, el("summary", { text: `Chebyshev degree ${p.degree}, ${p.coefficients.length} coefficients` }), list));
  }
  return box;
}

// ── the workspace: windows in an auto grid ────────────────────────────

const work = { items: [], next: 1 };

function addWindow(item, { focus = false } = {}) {
  item.id = work.next++;
  work.items.push(item);
  const win = item.kind === "ask" ? askWindow(item) : item.kind === "ref" ? refWindow(item) : galleryWindow(item);
  item.win = win;
  if (item.wide) win.classList.add("wide");
  $("#work").append(win);
  saveSoon();
  if (focus) {
    win.scrollIntoView({ block: "nearest", behavior: "smooth" });
    focusFirst(item);
  }
  return item;
}

function closeWindow(item) {
  work.items = work.items.filter((x) => x !== item);
  waiting.delete(item);
  item.win.remove();
  saveSoon();
  if (editor.active && !document.body.contains(editor.active)) { editor.active = null; setKeys(false); }
}

function head(item, title, ...extra) {
  return el("div", { class: "head" },
    el("h2", { class: "title", text: title }),
    ...extra,
    el("button", { class: "icon small", type: "button", "aria-label": `Close ${title}`, onclick: () => closeWindow(item) },
      svgClose()));
}

function svgClose() {
  const s = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  s.setAttribute("viewBox", "0 0 24 24");
  s.setAttribute("aria-hidden", "true");
  s.innerHTML = '<path d="M6 6l12 12M18 6 6 18"/>';
  return s;
}

function focusFirst(item) {
  if (item.kind === "ask") {
    if (item.mode === "formula") { const mf = item.win.querySelector("math-field"); if (mf) mf.focus(); }
    else if (item.mode === "text") item.win.querySelector(".text").focus();
  } else if (item.kind === "ref") item.win.querySelector(".filter").focus();
}

// ── an ask window ─────────────────────────────────────────────────────
// item: { kind: "ask", mode: "formula" | "text" | "fixed", rows: [latex], text, card, answer }

function askWindow(item) {
  if (item.mode === "formula" && !editor.ready) item.mode = "text";
  const modes = el("div", { class: "seg", role: "group", "aria-label": "Written as" });
  for (const [m, label] of [["formula", "Formula"], ["text", "Text"]]) {
    const b = el("button", { type: "button", value: m, text: label, onclick: () => setMode(item, m) });
    if (m === "formula" && !editor.ready) b.disabled = true;
    modes.append(b);
  }
  const win = el("section", { class: "win", "aria-label": "Question" },
    head(item, item.mode === "fixed" ? "From the gallery" : "Ask", item.mode === "fixed" ? null : modes),
    el("div", { class: "body" }, el("div", { class: "input" }), el("p", { class: "hint", "aria-live": "polite" }),
      el("div", { class: "answer", "aria-live": "polite" }), el("p", { class: "read-as" })));
  item.modes = modes;
  item.win = win;
  buildInput(item);
  return win;
}

function setMode(item, mode) {
  if (item.mode === mode) return;
  // the text a formula was read as is the core's, so a question carries over to text as the language writes it
  if (mode === "text" && !item.text && item.answer && item.answer.question) item.text = item.answer.question;
  item.mode = mode;
  buildInput(item);
  focusFirst(item);
  saveSoon();
}

function buildInput(item) {
  const box = item.win.querySelector(".input");
  for (const b of item.modes.querySelectorAll("button")) b.setAttribute("aria-pressed", String(b.value === item.mode));
  box.replaceChildren();
  item.win.querySelector(".hint").replaceChildren();
  if (item.mode === "formula") {
    if (!item.rows || !item.rows.length) item.rows = [""];
    const rows = el("div", { class: "rows" });
    box.append(rows);
    item.rows.forEach((latex, k) => rows.append(rowView(item, k, latex)));
  } else if (item.mode === "text") {
    const t = el("textarea", { class: "text", rows: "2", spellcheck: "false", autocapitalize: "off", "aria-label": "Question as text",
      placeholder: "solve x^2 - p = 0 for x in [0, 10], p in [0.5, 4]" });
    t.value = item.text || "";
    t.addEventListener("input", () => { item.text = t.value; complete.update(item, t); schedule(item); });
    t.addEventListener("keydown", (e) => complete.key(e, item, t));
    for (const ev of ["click", "keyup", "focus"]) t.addEventListener(ev, () => hint(item, t));
    t.addEventListener("blur", () => setTimeout(() => complete.close(), 120));
    box.append(t);
  } else {
    const typed = item.card.question && item.card.typed !== false;
    box.append(el("p", { class: typed ? "fixed" : "fixed said", text: item.card.question || "A torus whose shell thickness is the root of an equation in its radii" }));
  }
  recompute(item);
}

function rowView(item, k, latex) {
  const mf = document.createElement("math-field");
  mf.setAttribute("aria-label", `Row ${k + 1}`);
  mf.mathVirtualKeyboardPolicy = "manual";
  if (k === 0) mf.setAttribute("placeholder", "2x^2-3x-5=0");
  mf.value = latex;
  const row = el("div", { class: "row" }, el("span", { class: "n", text: String(k + 1), "aria-hidden": "true" }), mf);
  mf.addEventListener("input", () => { item.rows[indexOf(row)] = mf.value; schedule(item); });
  mf.addEventListener("focusin", () => { editor.active = mf; placeKeys(); });
  // Enter opens a row below; Backspace in an empty row closes it
  mf.addEventListener("keydown", (e) => {
    const at = indexOf(row);
    if (e.key === "Enter" && !e.shiftKey) {
      e.preventDefault();
      item.rows.splice(at + 1, 0, "");
      const next = rowView(item, at + 1, "");
      row.after(next);
      renumber(item);
      next.querySelector("math-field").focus();
      saveSoon();
    } else if (e.key === "Backspace" && mf.value === "" && item.rows.length > 1) {
      e.preventDefault();
      const prev = row.previousElementSibling || row.nextElementSibling;
      item.rows.splice(at, 1);
      row.remove();
      renumber(item);
      const pm = prev.querySelector("math-field");
      pm.focus();
      pm.executeCommand("moveToMathfieldEnd");
      schedule(item);
    }
  }, { capture: true });
  return row;
}

function indexOf(row) { return [...row.parentNode.children].indexOf(row); }

function renumber(item) {
  item.win.querySelectorAll(".row").forEach((r, k) => {
    r.querySelector(".n").textContent = String(k + 1);
    r.querySelector("math-field").setAttribute("aria-label", `Row ${k + 1}`);
  });
}

// live: each change asks the core again, a moment after the last keystroke
function schedule(item) {
  saveSoon();
  clearTimeout(item.timer);
  item.timer = setTimeout(() => recompute(item), 120);
}

async function recompute(item) {
  const out = item.win.querySelector(".answer");
  const readAs = item.win.querySelector(".read-as");
  for (const r of item.win.querySelectorAll(".row.refused")) r.classList.remove("refused");
  let kind, payload, typed = "", shown = [];
  if (item.mode === "fixed") {
    show(item, item.card, "");
    return;
  }
  if (item.mode === "formula") {
    const rows = [];
    item.rows.forEach((latex, k) => {
      if (!latex.trim()) return;
      rows.push(editor.parse(latex) ?? ["Error", "'missing'"]);
      shown.push(k);
    });
    if (!rows.length) { out.replaceChildren(); readAs.replaceChildren(); item.answer = null; return; }
    kind = "rows", payload = JSON.stringify({ rows });
  } else {
    typed = (item.text || "").trim();
    if (!typed) { out.replaceChildren(); readAs.replaceChildren(); item.answer = null; return; }
    kind = "text", payload = typed;
  }
  if (!core.ready) {
    // the gallery's card at once, and the live answer once the core is in
    const card = kind === "text" && data.gallery.find((c) => c.question === typed);
    if (card) show(item, card, typed);
    else if (core.failed) show(item, { refused: `the core did not load (${core.failed}), so only the gallery's questions are answered here` }, typed);
    else out.replaceChildren(el("p", { class: "state", text: "The core is loading; the answer follows." }));
    if (core.failed) return;
  }
  const a = await ask(item, kind, payload);
  if (a === null || !document.body.contains(item.win)) return;
  show(item, a, typed, shown);
}

function show(item, a, typed, shown = []) {
  item.answer = a;
  item.typed = typed;
  item.shown = shown;
  item.win.querySelector(".answer").replaceChildren(answerView(a, typed));
  if (!a.graph) item.win.querySelector(".graph")?.remove();
  // in formula rows, the text they were read as: the second view of the same question
  const readAs = item.win.querySelector(".read-as");
  readAs.replaceChildren();
  if (item.mode === "formula" && a.question) readAs.append(el("span", { class: "sr", text: "Read as: " }), a.question);
  if (a.graph) drawGraph(item, a.graph);
  if (a.row && shown[a.row - 1] !== undefined) {
    const r = item.win.querySelectorAll(".row")[shown[a.row - 1]];
    if (r) r.classList.add("refused");
  }
}

// ── a graph: the shader the core writes for a value of one parameter ──
// The page asks the core for the graph the answer names (`graph`, a request for spatium_plot),
// compiles the fragment shader it returns and draws it with WebGL2; the axes and their numbers
// are drawn on a 2D canvas over it. The y axis spans the heights the core reports (`range`).

const graphKeys = new WeakMap();   // a window -> the key its graph is asked under

async function drawGraph(item, request) {
  if (!core.ready) return;
  if (!graphKeys.has(item)) graphKeys.set(item, {});
  const p = await ask(graphKeys.get(item), "plot", JSON.stringify(request));
  if (!p || p.refused || item.answer?.graph !== request || !document.body.contains(item.win)) return;
  const fig = el("figure", { class: "graph" });
  const gl = el("canvas", { class: "gl", "aria-hidden": "true" }), axes = el("canvas", { class: "axes", "aria-hidden": "true" });
  const grow = el("button", { class: "icon small grow", type: "button", "aria-label": "Expand the graph", "aria-pressed": String(item.win.classList.contains("wide")) });
  grow.innerHTML = '<svg viewBox="0 0 24 24" aria-hidden="true"><path d="M4 9V4h5M15 4h5v5M20 15v5h-5M9 20H4v-5"/></svg>';
  grow.addEventListener("click", () => {
    const wide = item.win.classList.toggle("wide");
    item.wide = wide;
    saveSoon();
    grow.setAttribute("aria-pressed", String(wide));
    paint();
  });
  // the axes under the shader's curve; a caption only for the parameters held at their middles
  fig.append(el("div", { class: "plane" }, axes, gl), grow);
  if (p.symbols.length) fig.append(el("figcaption", { class: "label", text: p.symbols.map((s) => `${s.name} = ${number((s.interval[0] + s.interval[1]) / 2)}`).join(", ") }));
  const old = item.win.querySelector(".graph");
  if (old) old.replaceWith(fig); else item.win.querySelector(".answer").after(fig);
  const [x0, x1] = p.interval;
  let [y0, y1] = p.range || [-1, 1];
  if (!(y1 > y0)) { y0 -= 1; y1 += 1; }
  const pad = (y1 - y0) * 0.08;
  const view = [x0, y0 - pad, x1, y1 + pad];
  const ctx = gl.getContext("webgl2", { premultipliedAlpha: false, antialias: false });
  if (!ctx) { fig.replaceChildren(el("p", { class: "state", text: "This browser has no WebGL2 to draw the graph with." })); return; }
  const shader = (type, text) => { const s = ctx.createShader(type); ctx.shaderSource(s, text); ctx.compileShader(s); return s; };
  const prog = ctx.createProgram();
  ctx.attachShader(prog, shader(ctx.VERTEX_SHADER, p.vertex));
  ctx.attachShader(prog, shader(ctx.FRAGMENT_SHADER, p.fragment));
  ctx.linkProgram(prog);
  if (!ctx.getProgramParameter(prog, ctx.LINK_STATUS)) { fig.replaceChildren(el("p", { class: "state", text: `The graph's shader did not link: ${ctx.getProgramInfoLog(prog)}` })); return; }
  function paint() {
    const dpr = devicePixelRatio || 1, w = fig.querySelector(".plane").clientWidth, h = item.win.classList.contains("wide") ? 420 : 220;
    for (const c of [gl, axes]) { c.width = Math.round(w * dpr); c.height = Math.round(h * dpr); c.style.height = `${h}px`; }
    const css = getComputedStyle(document.documentElement);
    const ink = rgb(css.getPropertyValue("--accent"));
    ctx.viewport(0, 0, gl.width, gl.height);
    ctx.useProgram(prog);
    const u = (n) => ctx.getUniformLocation(prog, n);
    ctx.uniform4fv(u("plot_view"), view);
    ctx.uniform2f(u("plot_size"), gl.width, gl.height);
    ctx.uniform1f(u("plot_width"), 1.1 * dpr);
    ctx.uniform4f(u("plot_ink"), ink[0], ink[1], ink[2], 1);
    ctx.uniform1i(u("plot_values"), 0);
    if (p.symbols.length) ctx.uniform1fv(u("field_symbols"), new Float32Array(p.symbols.map((s) => (s.interval[0] + s.interval[1]) / 2)));
    if (p.weights.length) ctx.uniform4fv(u("field_weights"), new Float32Array(p.weights));
    ctx.clearColor(0, 0, 0, 0);
    ctx.clear(ctx.COLOR_BUFFER_BIT);
    ctx.enable(ctx.BLEND);
    ctx.blendFunc(ctx.SRC_ALPHA, ctx.ONE_MINUS_SRC_ALPHA);
    ctx.drawArrays(ctx.TRIANGLES, 0, 3);
    drawAxes(axes, view, dpr, css);
  }
  paint();
  new ResizeObserver(() => paint()).observe(fig.querySelector(".plane"));
}

function rgb(color) {
  const c = document.createElement("canvas").getContext("2d");
  c.fillStyle = color.trim();
  const hex = c.fillStyle.replace("#", "");
  return [0, 2, 4].map((k) => parseInt(hex.slice(k, k + 2), 16) / 255);
}

// round numbers along an axis, about one each 80 pixels
function ticks(a, b, pixels) {
  const raw = (b - a) / Math.max(2, pixels / 80), mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 5, 10].map((m) => m * mag).find((s) => s >= raw);
  const out = [];
  for (let t = Math.ceil(a / step) * step; t <= b + step * 1e-9; t += step) out.push(Math.abs(t) < step * 1e-9 ? 0 : t);
  return out;
}

function drawAxes(canvas, [x0, y0, x1, y1], dpr, css) {
  const g = canvas.getContext("2d"), W = canvas.width, H = canvas.height;
  const X = (x) => ((x - x0) / (x1 - x0)) * W, Y = (y) => H - ((y - y0) / (y1 - y0)) * H;
  g.clearRect(0, 0, W, H);
  g.font = `${12 * dpr}px ${css.getPropertyValue("--f-mono")}`;
  g.lineWidth = dpr;
  g.strokeStyle = css.getPropertyValue("--line-faint");
  g.fillStyle = css.getPropertyValue("--ink-3");
  for (const t of ticks(x0, x1, W / dpr)) {
    g.beginPath(); g.moveTo(X(t), 0); g.lineTo(X(t), H); g.stroke();
    g.fillText(String(Number(t.toPrecision(6))), X(t) + 4 * dpr, H - 6 * dpr);
  }
  for (const t of ticks(y0, y1, H / dpr)) {
    g.beginPath(); g.moveTo(0, Y(t)); g.lineTo(W, Y(t)); g.stroke();
    g.fillText(String(Number(t.toPrecision(6))), 4 * dpr, Y(t) - 4 * dpr);
  }
  g.strokeStyle = css.getPropertyValue("--ink-3");
  if (y0 < 0 && y1 > 0) { g.beginPath(); g.moveTo(0, Y(0)); g.lineTo(W, Y(0)); g.stroke(); }
  if (x0 < 0 && x1 > 0) { g.beginPath(); g.moveTo(X(0), 0); g.lineTo(X(0), H); g.stroke(); }
}

// the word under the caret in a line of text, and what the reference says of it
function wordAt(t) {
  const text = t.value, at = t.selectionStart ?? text.length;
  let a = at, b = at;
  while (a > 0 && /[A-Za-z_]/.test(text[a - 1])) a--;
  while (b < text.length && /[A-Za-z_]/.test(text[b])) b++;
  return { word: text.slice(a, b), from: a, to: b, before: text.slice(a, at) };
}

function hint(item, t) {
  const box = item.win.querySelector(".hint");
  box.replaceChildren();
  const r = data.reference;
  const w = wordAt(t).word;
  if (!r || !w) return;
  const f = r.functions.find((x) => x.name === w);
  if (f) { box.append(el("code", { text: `${f.name}(${f.arity === 2 ? "x, y" : "x"})` }), ` ${f.caption}. Defined for ${f.domain}.`); return; }
  const c = r.constants.find((x) => x.name === w);
  if (c) { box.append(el("code", { text: c.name }), ` ${c.caption}.`); return; }
  const q = r.questions.find((x) => x.form.split(" ")[0] === w);
  if (q) box.append(el("code", { text: q.form }), ` ${q.caption}.`);
}

// completion at the caret from the reference: Tab or Enter takes the word, Escape closes
const complete = {
  list: null, items: [], at: 0, item: null, t: null,
  words() {
    const r = data.reference;
    if (!r) return [];
    const out = [];
    for (const f of r.functions) out.push({ word: f.name, insert: `${f.name}(`, caption: f.caption });
    for (const c of r.constants) out.push({ word: c.name, insert: c.name, caption: c.caption });
    for (const w of ["solve", "integrate", "for", "in", "on", "integral"]) {
      const q = r.questions.find((x) => x.form.split(" ")[0] === w);
      out.push({ word: w, insert: `${w} `, caption: q ? q.caption : "" });
    }
    return out;
  },
  update(item, t) {
    const w = wordAt(t);
    const found = w.before.length >= 1 ? this.words().filter((x) => x.word.startsWith(w.before) && x.word !== w.word) : [];
    if (!found.length) { this.close(); return; }
    this.item = item, this.t = t, this.items = found.slice(0, 6), this.at = 0;
    if (!this.list) {
      this.list = el("ul", { class: "complete", role: "listbox", "aria-label": "Completions" });
      document.body.append(this.list);
    }
    this.list.replaceChildren(...this.items.map((x, k) => {
      const li = el("li", { role: "option", "aria-selected": String(k === 0) }, el("code", { text: x.word }), x.caption ? el("span", { text: x.caption }) : null);
      li.addEventListener("pointerdown", (e) => { e.preventDefault(); this.at = k; this.take(); });
      return li;
    }));
    const r = t.getBoundingClientRect();
    this.list.style.left = `${Math.max(8, Math.min(r.left, innerWidth - 368))}px`;
    this.list.style.top = `${r.bottom + 4}px`;
  },
  key(e, item, t) {
    if (!this.list || this.t !== t) return;
    if (e.key === "ArrowDown" || e.key === "ArrowUp") {
      e.preventDefault();
      this.at = (this.at + (e.key === "ArrowDown" ? 1 : this.items.length - 1)) % this.items.length;
      [...this.list.children].forEach((li, k) => li.setAttribute("aria-selected", String(k === this.at)));
    } else if (e.key === "Tab" || e.key === "Enter") {
      e.preventDefault();
      this.take();
    } else if (e.key === "Escape") this.close();
  },
  take() {
    const t = this.t, x = this.items[this.at], w = wordAt(t);
    t.setRangeText(x.insert, w.from, w.to, "end");
    this.item.text = t.value;
    this.close();
    hint(this.item, t);
    schedule(this.item);
  },
  close() { if (this.list) { this.list.remove(); this.list = null; } },
};

// ── a reference window ───────────────────────────────────────────────

function refWindow(item) {
  const r = data.reference;
  const list = el("div", { class: "ref" });
  const filter = el("input", { class: "filter", type: "search", placeholder: "Find a function, an operator, a form", "aria-label": "Find in the reference", spellcheck: "false" });
  const entry = (name, e, domain) => el("div", { class: "entry", "data-find": `${name} ${e.caption}`.toLowerCase() },
    el("div", { class: "name", text: name }),
    el("div", {},
      el("p", { class: "caption", text: e.caption + (domain ? `. Defined for ${domain}.` : "") }),
      e.example ? el("button", { class: "ex", type: "button", text: e.example, onclick: () => openText(e.example) }) : null));
  if (!r) list.append(el("p", { class: "state", text: "The reference did not load." }));
  else list.append(
    el("h3", { text: "Questions" }), ...r.questions.map((q) => entry(q.form, q)),
    el("h3", { text: "Functions" }), ...r.functions.map((f) => entry(`${f.name}(${f.arity === 2 ? "x, y" : "x"})`, f, f.domain)),
    el("h3", { text: "Operators" }), ...r.operators.map((o) => entry(o.written, o)),
    el("h3", { text: "Constants" }), ...r.constants.map((c) => entry(c.name, c)));
  filter.addEventListener("input", () => {
    const q = filter.value.trim().toLowerCase();
    for (const e of list.querySelectorAll(".entry")) e.hidden = Boolean(q) && !e.dataset.find.includes(q);
  });
  return el("section", { class: "win", "aria-label": "Reference" }, head(item, "Reference"), el("div", { class: "body" }, filter, list));
}

function openText(text) {
  addWindow({ kind: "ask", mode: "text", text }, { focus: true });
}

// ── a gallery window ─────────────────────────────────────────────────

function cardAnswer(c) {
  if (c.route === "scene") return null;
  if (c.refused) return null;
  if (c.route === "evaluate") return math(c.latex);
  if (c.route === "closed form") return math(finalStep(c.steps).latex);
  if (c.answer && c.answer.length === 1 && c.answer[0].form !== "chebyshev") return math(c.answer[0].latex);
  return null;
}

function galleryWindow(item) {
  const cards = el("div", { class: "cards" });
  for (const c of data.gallery) {
    const typed = c.question && c.typed !== false;
    const [ok, label] = kindOf(c);
    cards.append(el("button", { class: "card", type: "button", onclick: () => {
        if (typed) openText(c.question);
        else addWindow({ kind: "ask", mode: "fixed", card: c }, { focus: true });
      } },
      el("p", { class: typed ? "q" : "q said", text: c.question || "A torus whose shell thickness is the root of an equation in its radii" }),
      el("div", { class: "a" }, cardAnswer(c), el("p", { class: ok ? "kind ok" : "kind", text: label }))));
  }
  return el("section", { class: "win", "aria-label": "Gallery" }, head(item, "Gallery"), el("div", { class: "body" }, cards));
}

// ── the keyboard: MathLive's own, in a window near the row being written ─

function setKeys(open) {
  const keys = $("#keys");
  if (!keys) return;
  keys.hidden = !open;
  $("#open-keys").setAttribute("aria-expanded", String(open));
  if (!editor.ready) return;
  if (open) {
    if (!editor.active || !document.body.contains(editor.active)) {
      let item = work.items.find((x) => x.kind === "ask" && x.mode === "formula");
      if (!item) item = addWindow({ kind: "ask", mode: "formula", rows: [""] });
      editor.active = item.win.querySelector("math-field");
    }
    editor.active.focus();
    window.mathVirtualKeyboard.show();
    placeKeys();
  } else window.mathVirtualKeyboard.hide();
}

// below the active row when there is room, above it when not; on a phone, along the bottom (style.css)
function placeKeys() {
  const keys = $("#keys");
  if (keys.hidden || !editor.active) return;
  const r = editor.active.getBoundingClientRect(), h = keys.offsetHeight, w = keys.offsetWidth;
  const below = r.bottom + 8 + h <= innerHeight;
  keys.style.top = `${below ? r.bottom + 8 : Math.max(8, r.top - h - 8)}px`;
  keys.style.left = `${Math.max(16, Math.min(r.left, innerWidth - w - 16))}px`;
}

// ── the workspace as a link: #w= is its windows, kept in the address and in this browser ─
// A window is written as what it holds, an array led by a letter: ["f", rows...] a formula,
// ["t", text] a text, ["g", id] a gallery card, ["r"] the reference, ["G"] the gallery; "+" after
// the letter when its graph is expanded. JSON, deflate-raw, base64url. #solve?q=… opens one question as text.

const WORKSPACE = "spatium.workspace";

function snapshot() {
  const windows = work.items.map((item) => {
    const wide = item.wide ? "+" : "";
    if (item.kind === "ref") return ["r" + wide];
    if (item.kind === "gallery") return ["G" + wide];
    if (item.mode === "formula") return ["f" + wide, ...item.rows];
    if (item.mode === "text") return ["t" + wide, item.text || ""];
    return ["g" + wide, item.card.id];
  });
  return [1, ...windows];
}

// a window as the page opens it, or null for one it cannot (an unknown letter, a card no longer in the gallery)
function windowOf(w) {
  if (!Array.isArray(w) || typeof w[0] !== "string" || !w.slice(1).every((x) => typeof x === "string")) return null;
  const [letter, ...held] = w, wide = letter.endsWith("+");
  switch (letter.replace(/\+$/, "")) {
    case "r": return { kind: "ref", wide };
    case "G": return { kind: "gallery", wide };
    case "f": return { kind: "ask", mode: "formula", rows: held.length ? held : [""], wide };
    case "t": return { kind: "ask", mode: "text", text: held[0] ?? "", wide };
    case "g": {
      const card = data.gallery.find((c) => c.id === held[0]);
      return card ? { kind: "ask", mode: "fixed", card, wide } : null;
    }
  }
  return null;
}

function restore(state) {
  if (!Array.isArray(state) || state[0] !== 1) return false;
  const windows = state.slice(1).map(windowOf).filter(Boolean);
  for (const item of [...work.items]) closeWindow(item);
  for (const w of windows) addWindow(w);
  return true;
}

const base64url = (bytes) => btoa(String.fromCharCode(...bytes)).replace(/\+/g, "-").replace(/\//g, "_").replace(/=+$/, "");
const unbase64url = (s) => Uint8Array.from(atob(s.replace(/-/g, "+").replace(/_/g, "/")), (c) => c.charCodeAt(0));
const through = (bytes, stream) => new Response(new Blob([bytes]).stream().pipeThrough(stream)).arrayBuffer();

async function linkOf(state) {
  if (!("CompressionStream" in window)) return null;
  return base64url(new Uint8Array(await through(JSON.stringify(state), new CompressionStream("deflate-raw"))));
}

async function stateOf(code) {
  return JSON.parse(new TextDecoder().decode(await through(unbase64url(code), new DecompressionStream("deflate-raw"))));
}

// a moment after the last change: the address names the workspace, and this browser keeps it
let saveTimer, restoring = false;
function saveSoon() {
  if (restoring) return;
  clearTimeout(saveTimer);
  saveTimer = setTimeout(save, 300);
}

async function save() {
  const state = snapshot();
  try { localStorage.setItem(WORKSPACE, JSON.stringify(state)); } catch {}
  const code = await linkOf(state).catch(() => null);
  if (code !== null) history.replaceState(null, "", `#w=${code}`);
}

async function copyLink() {
  clearTimeout(saveTimer);
  await save();
  try {
    await navigator.clipboard.writeText(location.href);
    toast("Link to this workspace copied");
  } catch {
    toast("The address bar holds the link to this workspace");
  }
}

// what the address asks for: true when it opened something
async function fromHash() {
  const hash = location.hash.slice(1);
  if (hash.startsWith("w=")) {
    let state = null;
    try { state = await stateOf(hash.slice(2)); } catch {}
    if (opened(state)) return true;
    toast("The link to a workspace could not be read", true);
    return false;
  }
  const [path, query] = hash.split("?");
  if (path === "solve" || path === "") {
    const q = new URLSearchParams(query || "").get("q");
    if (q !== null) { openText(q); return true; }
  }
  return false;
}

function opened(state) {
  restoring = true;
  try { return restore(state); } finally { restoring = false; saveSoon(); }
}

function fromBrowser() {
  try { return opened(JSON.parse(localStorage.getItem(WORKSPACE) || "null")); } catch { return false; }
}

function rerender() {
  for (const item of work.items)
    if (item.kind === "ask" && item.answer) show(item, item.answer, item.typed, item.shown);
}

// ── wiring ────────────────────────────────────────────────────────────

function wire() {
  for (const b of document.querySelectorAll("[data-open]"))
    b.addEventListener("click", () => {
      const kind = b.dataset.open;
      addWindow(kind === "ask" ? { kind, mode: "formula", rows: [""] } : { kind }, { focus: true });
    });
  $("#open-keys").addEventListener("pointerdown", (e) => e.preventDefault());
  $("#open-keys").addEventListener("click", () => setKeys($("#keys").hidden));
  $("#close-keys").addEventListener("click", () => setKeys(false));
  document.addEventListener("keydown", (e) => { if (e.key === "Escape" && !$("#keys").hidden) setKeys(false); });
  addEventListener("resize", placeKeys);
  addEventListener("scroll", placeKeys, { passive: true });
  $("#open-settings").addEventListener("click", () => { showCoreState(); $("#settings").showModal(); });
  $("#settings").addEventListener("click", (e) => { if (e.target === $("#settings")) $("#settings").close(); });
  for (const seg of document.querySelectorAll(".seg[data-setting]"))
    seg.addEventListener("click", (e) => {
      const b = e.target.closest("button");
      if (!b) return;
      settings[seg.dataset.setting] = b.value;
      saveSettings();
      applySettings();
      if (["digits", "formulas"].includes(seg.dataset.setting)) rerender();
    });
  $("#check-update").addEventListener("click", checkUpdate);
  $("#copy-link").addEventListener("click", copyLink);
  addEventListener("hashchange", fromHash);
}

wire();
applySettings();
loadCore();
await Promise.all([loadData(), loadEditor()]);
applySettings();
// a link opens what it names, else the workspace this browser kept; otherwise two questions to
// start from: an equation solved with its steps, and a root as a function of a parameter,
// synthesised and certified
if (!(await fromHash()) && !fromBrowser()) {
  addWindow({ kind: "ask", mode: "formula", rows: ["2x^2-3x-5=0"] });
  addWindow({ kind: "ask", mode: "formula", rows: ["p\\in[0.5,4]", "x^2-p=0,\\ x\\in[0,10]"] });
}
