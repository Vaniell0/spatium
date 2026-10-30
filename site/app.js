// The page: it reads what the library answers and shows it. Nothing here computes an answer or
// writes a formula: the answer, its steps and its certificate come from the core
// (wasm/spatium_core.wasm, the library built to WebAssembly) as JSON, and the gallery and the
// reference are files the library wrote. The page chooses where each part goes.

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
  $("#open-keys").hidden = settings.keyboard !== "on";
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
function math(latex, { display = false, small = false } = {}) {
  const box = el("div", { class: small ? "math small" : "math" });
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
  return math(`10^{${Math.round(Math.log10(e))}}`, { small: true });
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
  }, 3500);
}

// ── the core ──────────────────────────────────────────────────────────

const core = { ready: false, failed: null, respond: null, preview: null, version: "", size: 0, from: "" };

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

// the glue and the module are kept and replaced together, so a cached module never meets a newer glue
async function startCore(glue, wasm) {
  const url = URL.createObjectURL(new Blob([glue], { type: "text/javascript" }));
  await new Promise((done, fail) => {
    const s = el("script", { src: url });
    s.onload = done;
    s.onerror = () => fail(new Error("the core's script did not run"));
    document.head.append(s);
  });
  URL.revokeObjectURL(url);
  const m = await window.SpatiumCore({ wasmBinary: new Uint8Array(wasm) });
  core.respond = m.cwrap("spatium_respond", "string", ["string", "number"]);
  core.preview = m.cwrap("spatium_preview", "string", ["string"]);
  core.version = await digest(wasm);
  core.size = wasm.byteLength;
  core.ready = true;
  core.failed = null;
}

async function loadCore() {
  try {
    const [glue, wasm] = await Promise.all([fetchBytes(CORE_JS), fetchBytes(CORE_WASM)]);
    await startCore(glue.bytes, wasm.bytes);
    core.from = wasm.from;
    if (wasm.from === "network") {
      await storeBytes(CORE_JS, glue.bytes, "text/javascript");
      await storeBytes(CORE_WASM, wasm.bytes, "application/wasm");
    }
    toast(`Core loaded, ${Math.round(core.size / 1024)} KB`);
  } catch (e) {
    core.failed = e.message || String(e);
    toast(`The core did not load: ${core.failed}`, true);
  }
  showCoreState();
  onCoreReady();
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
    await startCore(glue.bytes, wasm.bytes);
    core.from = "network";
    showCoreState();
    toast(`Core updated to ${core.version}`);
    onCoreReady();
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

// ── the answer sheet ──────────────────────────────────────────────────

const HOW = {
  evaluate: "Evaluated as written, in double precision, with nothing rewritten.",
  "closed form": "Solved in closed form. The coefficients are read from the equation's own values and the degree is checked at more points; the roots come from the quadratic formula in the form that keeps the small root, and Vieta's formulas check them by a second route.",
  synthesis: "Synthesised. The library computes the question's value at each parameter by its own solvers; the answer is the shortest expression in the parameter that matches it, searched bottom up, or a Chebyshev polynomial where none of a few terms does, and it is checked against those values on a dense grid.",
};

const COMPUTED = { wasm: "in this browser", host: "on the host, when the gallery was written" };

function resultRow(content, question) {
  const row = el("div", { class: "result" }, content);
  if (question) {
    const copy = el("button", { class: "quiet", type: "button", text: "Copy link" });
    copy.addEventListener("click", async () => {
      const url = `${location.origin}${location.pathname}#solve?q=${encodeURIComponent(question)}`;
      try { await navigator.clipboard.writeText(url); copy.textContent = "Copied"; }
      catch { copy.textContent = "Copy failed"; }
      setTimeout(() => { copy.textContent = "Copy link"; }, 1600);
    });
    row.append(copy);
  }
  return row;
}

function part(title, ...children) {
  return el("section", { class: "part" }, el("h2", { text: title }), ...children);
}

function facts(rows) {
  const dl = el("dl", { class: "facts" });
  for (const [k, v] of rows) dl.append(el("dt", { text: k }), el("dd", {}, v));
  return dl;
}

function refusal(reason) {
  return el("div", { class: "refusal" }, el("p", { text: "Refused" }), el("p", { text: reason }));
}

// the step a closed-form answer ends in: its roots, or the root of a linear equation
function finalStep(steps) {
  const named = steps.find((s) => s.label === "roots" || s.label === "root");
  if (named) return named;
  const solving = steps.filter((s) => !s.label.startsWith("check"));
  return solving[solving.length - 1];
}

function sheet(a, question) {
  const out = document.createDocumentFragment();
  const routeNote = HOW[a.route];

  if (a.refused) {
    out.append(refusal(a.refused));
    if (question) out.append(resultRow(el("span"), question));
    if (routeNote) out.append(part("How it was tried", el("p", { text: routeNote })));
    return out;
  }

  if (a.route === "evaluate") {
    out.append(resultRow(math(a.latex, { display: true }), question));
    out.append(facts([["Value", el("span", { class: "mono", text: number(a.value) })]]));
    out.append(part("How found", el("p", { text: routeNote })));
    return out;
  }

  if (a.route === "closed form") {
    const last = finalStep(a.steps);
    out.append(resultRow(math(last.latex, { display: true }), question));
    const chain = el("ol", { class: "chain" });
    if (a.check) chain.classList.add(a.check.holds ? "holds" : "fails");
    for (const s of a.steps) chain.append(el("li", {}, el("div", { class: "label", text: s.label }), math(s.latex, { display: true })));
    out.append(part("Steps", chain));
    if (a.check) out.append(facts([
      ["Check", el("span", { class: a.check.holds ? "holds-text" : "fails-text",
                             text: `${a.check.by} ${a.check.holds ? "hold" : "do not hold"}` })],
      [a.check.error_below ? "Relative error below" : "Error", decade(a.check.error_below)],
    ]));
    out.append(part("How found", el("p", { text: routeNote })));
    return out;
  }

  if (a.route === "synthesis") {
    const pieces = a.answer || [];
    if (pieces.length === 1) out.append(resultRow(math(pieces[0].latex, { display: true }), question));
    else out.append(resultRow(el("p", { class: "value", text: `${a.pieces} pieces` }), question));
    const c = a.certificate;
    if (c) out.append(part("Certificate", facts([
      ["Holds", el("span", { class: "holds-text", text: `${c.strength}, on ${c.points} points` })],
      [c.norm === "largest relative error" ? "Relative error below" : "Absolute error below", decade(c.error_below)],
      ["Scalar", el("span", { text: c.scalar })],
      ["Computed", el("span", { text: COMPUTED[c.computed] || c.computed })],
    ])));
    if (pieces.length) {
      const list = el("div");
      for (const p of pieces) list.append(pieceView(p, pieces.length > 1));
      out.append(part(pieces.length > 1 ? "Pieces" : "The answer", list));
    } else if (a.pieces) {
      out.append(part("Pieces", el("p", { text: `${a.pieces} pieces, too many to write out here.` })));
    }
    out.append(part("How found", el("p", { text: routeNote })));
    return out;
  }
  return out;
}

function pieceView(p, many) {
  const box = el("div", { class: "piece" });
  box.append(el("div", { class: "on", text: `on [${number(p.lo)}, ${number(p.hi)}]` }));
  if (many) box.append(math(p.latex));
  if (p.form === "expression") {
    box.append(facts([
      ["Terms", el("span", { text: String(p.terms) })],
      ["As text", el("code", { text: p.text })],
    ]));
  } else if (p.form === "chebyshev") {
    box.append(facts([["Chebyshev degree", el("span", { text: String(p.degree) })]]));
    if (p.coefficients) {
      const list = el("ol", { class: "coeffs" });
      p.coefficients.forEach((c, k) => list.append(el("li", {}, el("span", { text: `c${k}` }), number(c))));
      box.append(el("details", {}, el("summary", { text: `${p.coefficients.length} coefficients` }), list));
    }
  }
  return box;
}

// ── solve ─────────────────────────────────────────────────────────────

const input = () => $("#question");
let pending = null; // a question asked before the core loaded

function previewNow() {
  const text = input().value;
  const box = $("#preview");
  box.replaceChildren();
  hintNow();
  if (!text.trim() || !core.ready) return;
  const p = JSON.parse(core.preview(text));
  if (p.latex !== undefined) {
    box.append(math(p.latex));
    return;
  }
  const col = Number((/column (\d+)/.exec(p.error) || [])[1] || 0);
  const typed = el("div", { class: "typed" });
  if (col > 0) {
    const at = Math.min(col - 1, text.length);
    typed.append(text.slice(0, at), el("mark", { text: text.slice(at, at + 1) || " " }), text.slice(at + 1));
  } else typed.append(text);
  box.append(typed, el("p", { class: "why", text: p.error }));
}

function wordAtCaret() {
  const q = input();
  const text = q.value, at = q.selectionStart ?? text.length;
  let a = at, b = at;
  while (a > 0 && /[A-Za-z_]/.test(text[a - 1])) a--;
  while (b < text.length && /[A-Za-z_]/.test(text[b])) b++;
  return text.slice(a, b);
}

function hintNow() {
  const hint = $("#hint");
  hint.replaceChildren();
  const r = data.reference;
  const w = wordAtCaret();
  if (!r || !w) return;
  const f = r.functions.find((x) => x.name === w);
  if (f) {
    const args = f.arity === 2 ? "x, y" : "x";
    hint.append(el("code", { text: `${f.name}(${args})` }), ` ${f.caption}. Defined for ${f.domain}.`);
    return;
  }
  const c = r.constants.find((x) => x.name === w);
  if (c) {
    hint.append(el("code", { text: c.name }), ` ${c.caption}.`);
    return;
  }
  const form = r.questions.find((x) => x.form.split(" ")[0] === w);
  if (form) hint.append(el("code", { text: form.form }), ` ${form.caption}.`);
}

function ask(text, { focusResult = false } = {}) {
  const result = $("#result");
  const q = text.trim();
  if (!q) { result.replaceChildren(); return; }
  if (!core.ready) {
    const card = data.gallery.find((c) => c.question === q);
    // the gallery's card at once, and the live answer once the core is in
    if (card) { result.replaceChildren(sheet(card, q)); if (!core.failed) pending = q; return; }
    if (core.failed) {
      result.replaceChildren(refusal(`The core did not load (${core.failed}), so only the gallery's questions are answered here.`));
      return;
    }
    pending = q;
    result.replaceChildren(el("p", { class: "state", text: "The core is loading; the answer follows." }));
    return;
  }
  result.replaceChildren(el("p", { class: "state", text: "Working" }));
  $("#ask .go").disabled = true;
  // two frames so the working state is drawn before the core holds the thread
  requestAnimationFrame(() => requestAnimationFrame(() => {
    let a;
    try { a = JSON.parse(core.respond(q, 0)); }
    catch (e) { a = { refused: `the core failed: ${e.message || e}` }; }
    result.replaceChildren(sheet(a, q));
    $("#ask .go").disabled = false;
    if (focusResult) result.focus({ preventScroll: true });
  }));
}

function onCoreReady() {
  previewNow();
  if (pending) { const q = pending; pending = null; ask(q); }
}

// ── gallery ───────────────────────────────────────────────────────────

// 10^-12 as a person reads it in a line of text
function decadeText(e) {
  const sup = { "-": "\u207b", 0: "\u2070", 1: "\u00b9", 2: "\u00b2", 3: "\u00b3", 4: "\u2074", 5: "\u2075", 6: "\u2076", 7: "\u2077", 8: "\u2078", 9: "\u2079" };
  return "10" + [...String(Math.round(Math.log10(e)))].map((d) => sup[d]).join("");
}

function kindOf(c) {
  if (c.route === "scene") return ["scene", "A scene, as LaTeX"];
  if (c.refused) return ["refused", "Refused"];
  if (c.route === "evaluate") return ["evaluated", "Evaluated"];
  if (c.route === "closed form") return c.check && c.check.holds ? ["checked", "Closed form, checked by Vieta's formulas"] : ["solved", "Closed form"];
  const e = c.certificate.error_below;
  return ["certified", e ? `Certified, ${c.certificate.norm === "largest relative error" ? "relative " : ""}error below ${decadeText(e)}` : "Certified, exact on the grid"];
}

function cardAnswer(c) {
  if (c.route === "scene") return null;
  if (c.refused) return el("p", { class: "fails-text", text: c.refused });
  if (c.route === "evaluate") return math(c.latex, { small: true });
  if (c.route === "closed form") return math(finalStep(c.steps).latex, { small: true });
  if (c.answer && c.answer.length === 1) {
    const p = c.answer[0];
    if (p.form === "chebyshev") return el("p", { text: `A Chebyshev series of degree ${p.degree}` });
    return math(p.latex, { small: true });
  }
  return el("p", { text: `${c.pieces} pieces` });
}

function cardHref(c) {
  if (c.route === "scene") return "#scene";
  return c.typed === false ? `#gallery/${c.id}` : `#solve?q=${encodeURIComponent(c.question)}`;
}

function renderGallery(id) {
  const cards = $("#cards"), detail = $("#card-sheet");
  const card = id && data.gallery.find((c) => c.id === id);
  cards.hidden = Boolean(card);
  detail.hidden = !card;
  if (card) {
    detail.replaceChildren(
      el("a", { class: "back", href: "#gallery", text: "All cards" }),
      el("p", { class: card.typed === false ? "asked said" : "asked mono", text: card.question }),
      sheet(card, null),
    );
    return;
  }
  cards.replaceChildren(...data.gallery.map((c) => {
    const [kind, label] = kindOf(c);
    const typed = c.question && c.typed !== false;
    return el("a", { class: "card", href: cardHref(c), "data-kind": kind },
      el("p", { class: typed ? "q" : "q said", text: c.question || "A torus whose shell thickness is the root of an equation in its radii" }),
      el("div", { class: "a" }, cardAnswer(c), el("p", { class: "kind", text: label })));
  }));
}

// ── reference ─────────────────────────────────────────────────────────

function example(e) {
  if (!e.example) return null;
  const box = el("div", { class: "example" }, el("a", { href: `#solve?q=${encodeURIComponent(e.example)}`, text: e.example }));
  const a = e.answer;
  if (a && a.route === "evaluate") box.append(math(`${a.latex} = ${a.value}`, { small: true }));
  return box;
}

function entry(name, e, domain) {
  return el("div", { class: "entry" },
    el("div", { class: "name", text: name }),
    el("div", {},
      el("p", { class: "caption", text: e.caption }),
      domain ? el("p", { class: "domain", text: `Defined for ${domain}` }) : null,
      example(e)));
}

function renderReference() {
  const r = data.reference;
  const box = $("#reference");
  if (!r) { box.replaceChildren(el("p", { class: "state", text: "The reference did not load." })); return; }
  box.replaceChildren(
    el("div", { class: "about" },
      el("p", {}, el("strong", { text: "Closed-form solvers with certified steps." }),
        " An expression is evaluated; an equation of degree at most two is solved by a fixed method, every step shown, and the roots checked by a second route."),
      el("p", {}, el("strong", { text: "Synthesis, which is not a general symbolic engine." }),
        " Nothing is rewritten symbolically. A question with a parameter is answered by the library's solvers number by number, and the answer is an expression found by search and certified against them, or refused with the reason.")),
    el("section", { class: "ref" }, el("h2", { text: "Questions" }),
      ...r.questions.map((q) => entry(q.form, q))),
    el("section", { class: "ref" }, el("h2", { text: "Functions" }),
      ...r.functions.map((f) => entry(`${f.name}(${f.arity === 2 ? "x, y" : "x"})`, f, f.domain))),
    el("section", { class: "ref" }, el("h2", { text: "Operators" }),
      ...r.operators.map((o) => entry(o.written, o))),
    el("section", { class: "ref" }, el("h2", { text: "Constants" }),
      ...r.constants.map((c) => entry(c.name, c))),
  );
}

// ── scene ─────────────────────────────────────────────────────────────

function renderScene() {
  const scene = data.gallery.find((c) => c.route === "scene");
  $("#scene").replaceChildren(
    scene ? math(scene.latex, { display: true }) : el("p", { class: "state", text: "No scene in the gallery." }),
    el("div", { class: "window", text: "The scene renderer arrives with the scene editor." }));
}

// ── keyboard ──────────────────────────────────────────────────────────

const KEYS = [
  ["sin", "sin("], ["cos", "cos("], ["log", "log("], ["exp", "exp("], ["sqrt", "sqrt("],
  ["min", "min("], ["max", "max("], ["less", "less("], ["pi", "pi"],
  ["x", "x"], ["p", "p"], ["t", "t"], ["^", "^"], ["(", "("], [")", ")"],
  ["*", " * "], ["/", " / "], ["+", " + "], ["-", " - "], ["=", " = "],
  ["[", "["], ["]", "]"], [",", ", "],
  ["solve", "solve "], ["integrate", "integrate "], ["for", " for "], ["in", " in "],
  ["⌫", null],
];

function insert(text) {
  const q = input();
  const a = q.selectionStart ?? q.value.length, b = q.selectionEnd ?? a;
  if (text === null) {
    const from = a === b ? Math.max(0, a - 1) : a;
    q.setRangeText("", from, b, "end");
  } else q.setRangeText(text, a, b, "end");
  previewNow();
}

function setKeys(open) {
  const keys = $("#keys");
  keys.hidden = !open;
  document.body.classList.toggle("keys-open", open);
  $("#open-keys").setAttribute("aria-expanded", String(open));
}

function buildKeys() {
  $("#keys").replaceChildren(...KEYS.map(([label, text]) => {
    const b = el("button", { type: "button", text: label, "aria-label": text === null ? "Delete" : label });
    // pointerdown keeps the caret in the question line
    b.addEventListener("pointerdown", (e) => e.preventDefault());
    b.addEventListener("click", () => { input().focus(); insert(text); });
    return b;
  }));
}

// ── routes: #solve?q=…, #gallery/<id>, #reference, #scene ─────────────

function route() {
  const [path, query] = location.hash.slice(1).split("?");
  const [tab, id] = path.split("/");
  const name = ["solve", "gallery", "reference", "scene"].includes(tab) ? tab : "solve";
  for (const s of document.querySelectorAll("[data-view]")) s.hidden = s.dataset.view !== name;
  for (const a of document.querySelectorAll(".tabs a"))
    if (a.dataset.tab === name) a.setAttribute("aria-current", "page"); else a.removeAttribute("aria-current");
  if (name !== "solve") document.body.classList.remove("focus");
  if (name === "solve") {
    const q = new URLSearchParams(query || "").get("q");
    if (q !== null && q !== input().value) {
      input().value = q;
      previewNow();
      ask(q);
    }
  } else if (name === "gallery") renderGallery(id);
  else if (name === "reference") renderReference();
  else renderScene();
}

function rerender() {
  previewNow();
  const q = input().value.trim();
  if (q && $("#result").childElementCount) ask(q);
  route();
}

// ── wiring ────────────────────────────────────────────────────────────

function wire() {
  $("#ask").addEventListener("submit", (e) => {
    e.preventDefault();
    const q = input().value.trim();
    history.replaceState(null, "", q ? `#solve?q=${encodeURIComponent(q)}` : "#solve");
    ask(q, { focusResult: true });
  });
  input().addEventListener("input", previewNow);
  for (const ev of ["keyup", "click", "focus"]) input().addEventListener(ev, hintNow);
  $("#open-keys").addEventListener("click", () => setKeys($("#keys").hidden));
  $("#focus").addEventListener("click", () => {
    const on = document.body.classList.toggle("focus");
    $("#focus").setAttribute("aria-pressed", String(on));
  });
  document.addEventListener("keydown", (e) => {
    if (e.key !== "Escape") return;
    if (document.body.classList.contains("focus")) {
      document.body.classList.remove("focus");
      $("#focus").setAttribute("aria-pressed", "false");
    } else if (!$("#keys").hidden) setKeys(false);
  });
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
  window.addEventListener("hashchange", route);
}

buildKeys();
applySettings();
wire();
loadData().then(route);
loadCore();
