// The calculator page held against the design document's criteria in a headless browser: one
// scenario per criterion, PASS or FAIL each, and the acceptance shots at 400 and 1280 wide in
// both themes. Node alone drives Chromium over the DevTools protocol; nothing is installed.
//   node site/check/criteria.mjs <bundle dir> <shots dir>
// where <bundle dir> is what site/bundle.sh laid out. Exit code 1 when a criterion fails.
// Each criterion names the section of the design document it measures.

import { spawn } from "node:child_process";
import { createServer } from "node:http";
import { existsSync, mkdirSync, mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { extname, join } from "node:path";
import { inflateSync } from "node:zlib";

const [bundle, shots] = process.argv.slice(2);
if (!bundle || !shots) { console.error("usage: node site/check/criteria.mjs <bundle dir> <shots dir>"); process.exit(2); }
mkdirSync(shots, { recursive: true });

// ── the bundle, served ──────────────────────────────────────────────────

const TYPES = { ".html": "text/html", ".js": "text/javascript", ".mjs": "text/javascript", ".css": "text/css", ".json": "application/json",
  ".wasm": "application/wasm", ".woff2": "font/woff2", ".woff": "font/woff", ".ttf": "font/ttf", ".svg": "image/svg+xml", ".png": "image/png" };
const server = createServer((req, res) => {
  let path = decodeURIComponent(new URL(req.url, "http://x").pathname);
  if (path.endsWith("/")) path += "index.html";
  const file = join(bundle, path);
  if (!file.startsWith(bundle) || !existsSync(file)) { res.writeHead(404); res.end(); return; }
  res.writeHead(200, { "content-type": TYPES[extname(file)] ?? "application/octet-stream" });
  res.end(readFileSync(file));
});
await new Promise((r) => server.listen(0, "127.0.0.1", r));
const origin = `http://127.0.0.1:${server.address().port}/`;

// ── Chromium over the DevTools protocol ─────────────────────────────────

const profile = mkdtempSync(join(shots, "profile-"));
const chrome = spawn("chromium", ["--headless=new", "--disable-gpu", "--no-sandbox", "--hide-scrollbars", "--window-size=1280,900",
  "--remote-debugging-port=0", `--user-data-dir=${profile}`, "about:blank"], { stdio: "ignore" });
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
let ws, id = 0;
const waiting = new Map();
for (let k = 0; k < 100 && !ws; ++k) {
  await sleep(200);
  const at = join(profile, "DevToolsActivePort");
  if (!existsSync(at)) continue;
  const port = readFileSync(at, "utf8").split("\n")[0];
  try {
    const targets = await (await fetch(`http://127.0.0.1:${port}/json`)).json();
    const page = targets.find((t) => t.type === "page");
    if (page) ws = new WebSocket(page.webSocketDebuggerUrl);
  } catch { /* not up yet */ }
}
if (!ws) { console.error("chromium did not come up"); process.exit(2); }
await new Promise((r) => (ws.onopen = r));
const thrown = [];   // the page's own uncaught errors
ws.onmessage = (e) => {
  const m = JSON.parse(e.data);
  if (m.id && waiting.has(m.id)) { waiting.get(m.id)(m); waiting.delete(m.id); }
  else if (m.method === "Runtime.exceptionThrown") thrown.push((m.params.exceptionDetails.exception?.description ?? m.params.exceptionDetails.text).split("\n")[0]);
};
const send = (method, params = {}) => new Promise((r) => { const i = ++id; waiting.set(i, r); ws.send(JSON.stringify({ id: i, method, params })); });
const js = async (expr) => {
  const r = await send("Runtime.evaluate", { expression: expr, awaitPromise: true, returnByValue: true });
  if (r.result?.exceptionDetails) throw new Error(r.result.exceptionDetails.exception?.description ?? "evaluation failed");
  return r.result?.result?.value;
};
const key = async (k, mods = 0) => {
  const code = { Escape: 27, Enter: 13, Backspace: 8, ArrowUp: 38 }[k];
  const at = { key: k, code: k, windowsVirtualKeyCode: code, nativeVirtualKeyCode: code, modifiers: mods };
  await send("Input.dispatchKeyEvent", { type: "keyDown", ...at });
  await send("Input.dispatchKeyEvent", { type: "keyUp", ...at });
};
const click = async (x, y) => {
  await send("Input.dispatchMouseEvent", { type: "mouseMoved", x, y });
  for (const type of ["mousePressed", "mouseReleased"])
    await send("Input.dispatchMouseEvent", { type, x, y, button: "left", buttons: type === "mousePressed" ? 1 : 0, clickCount: 1 });
};
const clickOn = async (selector) => {
  const r = await js(`(() => { const e = document.querySelector(${JSON.stringify(selector)}); if (!e) return null; e.scrollIntoView({ block: "center" }); const b = e.getBoundingClientRect(); return [b.x + b.width / 2, b.y + b.height / 2]; })()`);
  if (!r) throw new Error(`nothing at ${selector}`);
  await click(r[0], r[1]);
};
const viewport = (width, height = 900) => send("Emulation.setDeviceMetricsOverride", { width, height, deviceScaleFactor: 1, mobile: width < 600 });
const theme = (which) => send("Emulation.setEmulatedMedia", { features: [{ name: "prefers-color-scheme", value: which }] });

await send("Page.enable");
await send("Runtime.enable");
// the worker's questions can be held back, so a slow answer and a stale one can be staged
await send("Page.addScriptToEvaluateOnNewDocument", { source: `
  window.__hold = [];   // ms to hold each next question, consumed in order
  const W = window.Worker;
  window.Worker = class extends W {
    postMessage(m, ...rest) {
      const hold = m && m.type === "ask" ? window.__hold.shift() : undefined;
      if (hold) setTimeout(() => super.postMessage(m, ...rest), hold);
      else super.postMessage(m, ...rest);
    }
  };` });

// ── the page: opened, written in, asked ─────────────────────────────────

async function open(hash = "", clean = true) {
  if (clean) {
    await send("Page.navigate", { url: origin + "index.html" });
    await sleep(300);
    await js("localStorage.clear(); sessionStorage.clear()");
  }
  await send("Page.navigate", { url: origin + "index.html" + hash });
  // the core and the editor loaded: the first answer drawn, or the line ready
  for (let k = 0; k < 60; ++k) {
    await sleep(250);
    if (await js("Boolean(document.querySelector('math-field')) && document.querySelector('#core-dot')?.classList.contains('ready')")) break;
  }
  await sleep(300);
}

// the line written as the rows given, and its answer waited for
async function write(latex, { settle = 4000 } = {}) {
  await js(`(() => { const mf = document.querySelector("math-field"); mf.focus(); mf.value = ${JSON.stringify(latex)}; mf.dispatchEvent(new Event("input", { bubbles: true })); })()`);
  for (let k = 0; k < settle / 100; ++k) {
    await sleep(100);
    if (await js(`!document.querySelector(".card.pending") && Boolean(document.querySelector(".card .answer > *"))`)) break;
  }
  await sleep(150);
}

const line = () => js(`document.querySelector("math-field")?.value ?? null`);
const answerText = () => js(`document.querySelector(".card .answer")?.textContent ?? ""`);
const shot = async (name) => {
  const r = await send("Page.captureScreenshot", { format: "png" });
  writeFileSync(join(shots, `${name}.png`), Buffer.from(r.result.data, "base64"));
  return Buffer.from(r.result.data, "base64");
};

// ── a PNG read back: the share of pixels near a colour ──────────────────

function pixels(png) {
  let at = 8, width = 0, height = 0;
  const idat = [];
  while (at < png.length) {
    const len = png.readUInt32BE(at), type = png.toString("ascii", at + 4, at + 8);
    if (type === "IHDR") { width = png.readUInt32BE(at + 8); height = png.readUInt32BE(at + 12); }
    if (type === "IDAT") idat.push(png.subarray(at + 8, at + 8 + len));
    at += 12 + len;
  }
  const raw = inflateSync(Buffer.concat(idat)), bpp = 4, stride = width * bpp, out = Buffer.alloc(stride * height);
  for (let y = 0, src = 0; y < height; ++y) {
    const filter = raw[src++], row = y * stride, up = row - stride;
    for (let x = 0; x < stride; ++x, ++src) {
      const a = x >= bpp ? out[row + x - bpp] : 0, b = y ? out[up + x] : 0, c = x >= bpp && y ? out[up + x - bpp] : 0;
      let v = raw[src];
      if (filter === 1) v += a; else if (filter === 2) v += b; else if (filter === 3) v += (a + b) >> 1;
      else if (filter === 4) { const p = a + b - c, pa = Math.abs(p - a), pb = Math.abs(p - b), pc = Math.abs(p - c); v += pa <= pb && pa <= pc ? a : pb <= pc ? b : c; }
      out[row + x] = v & 255;
    }
  }
  return { width, height, data: out };
}
function share(png, [r, g, b], within = 40) {
  const { width, height, data } = pixels(png);
  let near = 0;
  for (let k = 0; k < data.length; k += 4)
    if (Math.abs(data[k] - r) + Math.abs(data[k + 1] - g) + Math.abs(data[k + 2] - b) <= within) ++near;
  return near / (width * height);
}
const rgb = (hex) => [1, 3, 5].map((k) => parseInt(hex.slice(k, k + 2), 16));

// ── contrast, as WCAG counts it ─────────────────────────────────────────

const luminance = ([r, g, b]) => [r, g, b].map((c) => { c /= 255; return c <= .03928 ? c / 12.92 : ((c + .055) / 1.055) ** 2.4; })
  .reduce((s, c, k) => s + c * [.2126, .7152, .0722][k], 0);
const contrast = (a, b) => { const [x, y] = [luminance(a), luminance(b)].sort((p, q) => q - p); return (x + .05) / (y + .05); };
const token = async (name) => js(`getComputedStyle(document.documentElement).getPropertyValue(${JSON.stringify(name)}).trim()`);

// every element with its own text, with the smallest font it is set in; hidden and empty ones
// skipped, and the inside of a rendered formula, whose sub- and superscripts are the formula's
// own: the formula itself is measured at its root
const SMALLEST = `(() => { let out = []; for (const e of document.querySelectorAll("body *")) {
  if (!e.checkVisibility || !e.checkVisibility()) continue;
  if (e.closest(".katex") && !e.classList.contains("katex")) continue;
  if (![...e.childNodes].some((n) => n.nodeType === 3 && n.textContent.trim())) continue;
  const px = parseFloat(getComputedStyle(e).fontSize); if (px < 14) out.push([e.tagName.toLowerCase() + "." + e.className, px]); }
  return out.slice(0, 8); })()`;

// ── the criteria ────────────────────────────────────────────────────────

const results = [];
async function criterion(name, run) {
  let why = "";
  try { why = (await run()) ?? ""; } catch (e) { why = `threw: ${e.message}`; }
  results.push([name, why === ""]);
  console.log(`${why === "" ? "PASS" : "FAIL"}  ${name}${why ? `  -- ${why}` : ""}`);
}

const QUAD = "2x^2-3x-5=0";

await viewport(1280);
await theme("light");
await open();

// The page loads: its script threw nothing, the editor and the core are in
await criterion("page loads without an error", async () => {
  if (thrown.length) return thrown[0];
  return (await js(`Boolean(document.querySelector("math-field")) && document.querySelector("#core-dot")?.classList.contains("ready")`)) ? "" : "the line or the core is missing";
});

// Acceptance shots and measured checks: accent pixels on the first view at or below 2 %
await criterion("first view: accent pixels at or below 2 %", async () => {
  const png = await shot("first-1280-light");
  const p = share(png, rgb(await token("--accent")));
  return p <= .02 ? "" : `${(p * 100).toFixed(1)} % of the first view is accent`;
});

// Calculator, states: the first visit is an empty line with a placeholder and three examples
await criterion("first visit: empty line, placeholder, examples", async () => {
  const v = await line();
  if (v !== "") return `the line starts with ${JSON.stringify(v)}`;
  if (!(await js(`document.querySelector("math-field")?.getAttribute("placeholder")`))) return "no placeholder";
  const n = await js(`document.querySelectorAll(".examples .cell").length`);
  return n === 3 ? "" : `${n} example cells, not 3`;
});

// Five rules, 1: one question. The calculator has one line and no card stack, adder or derived cell
await criterion("one question: one line, no adder, no second card", async () => {
  const [cards, adder, fields] = await js(`[document.querySelectorAll(".card").length, document.querySelectorAll(".adder").length, document.querySelectorAll("math-field").length]`);
  return cards === 1 && adder === 0 && fields === 1 ? "" : `${cards} cards, ${adder} adders, ${fields} fields`;
});

// Complaint 11: per-card settings and the gear are gone
await criterion("no per-card settings, no gear", async () => {
  const n = await js(`document.querySelectorAll(".gear, .options").length`);
  return n === 0 ? "" : `${n} gear or options elements`;
});

// Icons: no title tooltips anywhere; every icon svg carries width and height attributes
await criterion("no title tooltips", async () => {
  const n = await js(`[...document.querySelectorAll("[title]")].map((e) => e.id || e.className).slice(0, 5)`);
  return n.length === 0 ? "" : `title on ${n.join(", ")}`;
});
await criterion("icon svgs sized by attribute, 20 x 20", async () => {
  const bad = await js(`[...document.querySelectorAll(".icon svg")].filter((s) => s.getAttribute("width") !== "20" || s.getAttribute("height") !== "20").length`);
  return bad === 0 ? "" : `${bad} icon svgs without width=20 height=20`;
});

// Icons: Settings uses the sliders icon, three lines with knobs, not the gear path
await criterion("settings icon is sliders", async () => {
  const d = await js(`document.querySelector("#open-settings svg")?.innerHTML ?? ""`);
  return d.includes("sliders") || /M4 (6|7)h/.test(d) ? "" : "the settings button draws the gear";
});

// Type: the line's formula at 1.25 x body, controls 40 tall
await criterion("line formula at 1.25 x body (21.25 px), controls 40 tall", async () => {
  const [mf, icons] = await js(`[parseFloat(getComputedStyle(document.querySelector("math-field")).fontSize),
    [...document.querySelectorAll(".bar .icon")].map((b) => b.getBoundingClientRect().height)]`);
  if (Math.abs(mf - 21.25) > .5) return `the line is set at ${mf} px`;
  const off = icons.filter((h) => Math.abs(h - 40) > .5);
  return off.length ? `icon buttons at ${off.join(", ")} px` : "";
});

// The line: the hint is a popover that takes no space in the flow
await criterion("hint popover takes no space", async () => {
  const before = await js(`document.querySelector(".cell > .ask").getBoundingClientRect().height`);
  await js(`(() => { const mf = document.querySelector("math-field"); mf.focus(); mf.value = "\\\\sin"; mf.dispatchEvent(new Event("input", { bubbles: true })); mf.dispatchEvent(new Event("selection-change")); })()`);
  await sleep(400);
  const [shown, after, pos] = await js(`[Boolean(document.querySelector(".ask .hint")?.textContent.trim()), document.querySelector(".cell > .ask").getBoundingClientRect().height, getComputedStyle(document.querySelector(".ask .hint")).position]`);
  await viewport(1280);
  await shot("typing-hint-1280-light");
  if (!shown) return "no hint for sin";
  return pos === "absolute" && Math.abs(after - before) < 1 ? "" : `hint is ${pos}, ask row ${before} -> ${after}`;
});

// The quadratic answered: steps open, as the reference's chain
await write(QUAD);
await criterion("quadratic: answered with a result and the chain", async () => {
  const [res, steps] = await js(`[document.querySelector(".card .result .math")?.textContent ?? "", document.querySelectorAll(".card ol.chain li").length]`);
  return res && steps >= 3 ? "" : `result ${JSON.stringify(res)}, ${steps} steps`;
});

// Five rules, 2: the steps are rendered open, and cannot be folded
await criterion("steps rendered open and never fold", async () => {
  const sel = `[...document.querySelectorAll(".card details.sec")].find((d) => d.querySelector("summary")?.textContent.trim().toLowerCase().startsWith("steps"))`;
  const open = await js(`${sel}?.open ?? null`);
  if (open !== true) return open === null ? "no Steps section" : "Steps is closed";
  await js(`${sel}.querySelector("summary").click()`);
  await sleep(100);
  return (await js(`${sel}?.open`)) ? "" : "Steps folded on a click";
});

// Complaint 4: no graph hangs after the answer; Plot is a section, closed by default
await criterion("plot is a closed section, no graph after the answer", async () => {
  const [plots, open, loose] = await js(`(() => { const d = [...document.querySelectorAll(".card details.sec")].filter((d) => d.querySelector("summary")?.textContent.trim().toLowerCase().startsWith("plot"));
    return [d.length, d.some((x) => x.open), document.querySelectorAll(".card > .also > *, .card > .graph").length]; })()`);
  if (loose) return `${loose} graph elements outside the sheet`;
  return plots === 1 && !open ? "" : `${plots} plot sections, open ${open}`;
});

// Answer sheet anatomy: the cell measured against answer-sheet.html, value for value
await criterion("sheet measures as the reference", async () => {
  const m = await js(`(() => { const cs = (s) => getComputedStyle(document.querySelector(s)); const out = [];
    const cell = cs(".card.cell"); if (cell.borderTopWidth !== "1px" || cell.borderTopLeftRadius !== "8px") out.push("cell border " + cell.borderTopWidth + " radius " + cell.borderTopLeftRadius);
    const r = cs(".card .result"); if (r.paddingLeft !== "24px" || r.rowGap !== "8px") out.push("result padding " + r.paddingLeft + " gap " + r.rowGap);
    const rm = cs(".card .result .math"); if (Math.abs(parseFloat(rm.fontSize) - 17 * 1.35) > .5) out.push("result at " + rm.fontSize);
    const ink = getComputedStyle(document.documentElement).getPropertyValue("--ink").trim();
    const st = cs(".card .verdict .strength"), ac = getComputedStyle(document.documentElement).getPropertyValue("--accent").trim();
    const hex = (c) => "#" + c.match(/\\d+/g).slice(0, 3).map((n) => (+n).toString(16).padStart(2, "0")).join("");
    if (hex(st.color) !== ac) out.push("strength in " + hex(st.color) + " not accent " + ac);
    if (hex(rm.color) !== ink) out.push("result in " + hex(rm.color) + " not ink");
    const su = cs(".card details.sec > summary"); if (su.textTransform !== "uppercase" || su.paddingLeft !== "24px" || Math.abs(parseFloat(su.letterSpacing) - .56) > .05) out.push("summary " + su.textTransform + " " + su.paddingLeft + " " + su.letterSpacing);
    const li = cs(".card ol.chain li"); if (li.gridTemplateColumns.split(" ")[0] !== "32px" || li.borderBottomStyle !== "dashed" || li.paddingTop !== "12px") out.push("chain row " + li.gridTemplateColumns + " " + li.borderBottomStyle + " " + li.paddingTop);
    const sm = cs(".card ol.chain li .math"); if (Math.abs(parseFloat(sm.fontSize) - 17) > .5) out.push("step formula at " + sm.fontSize);
    const op = cs(".card ol.chain li .op"); if (Math.abs(parseFloat(op.fontSize) - 14) > .5) out.push("step op at " + op.fontSize);
    const body = cs(".card details.sec .body"); if (body.paddingBottom !== "24px" || body.rowGap !== "16px") out.push("body " + body.paddingBottom + " " + body.rowGap);
    const dl = document.querySelector(".card dl.facts"); if (dl) { const f = getComputedStyle(dl); if (f.columnGap !== "24px" || f.rowGap !== "8px") out.push("facts gap " + f.rowGap + " " + f.columnGap); }
    return out; })()`);
  return m.join("; ");
});

// Numbers: no run of trailing zeros, no 16-digit tail, no -0 in a result or a step (the core's rule)
async function numbers() {
  return js(`(() => { const t = [...document.querySelectorAll(".card .result .math .katex-mathml annotation, .card ol.chain .math .katex-mathml annotation")].map((a) => a.textContent);
    return t.filter((s) => /\\d\\.\\d*?0{2,}(?![\\d])/.test(s) || /\\d{16,}/.test(s) || /(?<![\\d.])-0(?![\\d.])/.test(s)).slice(0, 3); })()`);
}
await criterion("numbers: no trailing zeros, 16-digit tails or -0 (quadratic)", async () => (await numbers()).join(" | "));

// Measured checks: no text under 14 px, at 1280
await criterion("no text under 14 px at 1280", async () => (await js(SMALLEST)).map(([e, px]) => `${e} ${px}px`).join(", "));

// Measured checks: every text token pair at or above 4.5:1 in both themes; dark is warm, not blue
await criterion("text token contrast at or above 4.5:1, both themes", async () => {
  const pairs = [["--ink", "--surface"], ["--ink-2", "--surface"], ["--ink-3", "--surface"], ["--accent", "--surface"], ["--ok", "--surface"], ["--warn", "--surface"],
    ["--ink-3", "--bg"], ["--accent", "--bg"], ["--ok", "--bg"], ["--ink-3", "--raised"], ["--accent", "--accent-soft"], ["--warn", "--warn-soft"], ["--ink", "--raised"], ["--ink-2", "--raised"]];
  const low = [];
  for (const which of ["light", "dark"]) {
    await theme(which);
    await sleep(100);
    for (const [fg, bg] of pairs) {
      const c = contrast(rgb(await token(fg)), rgb(await token(bg)));
      if (c < 4.5) low.push(`${which} ${fg} on ${bg} ${c.toFixed(2)}`);
    }
    if (which === "dark") {
      const [r, g, b] = rgb(await token("--bg"));
      if (!(r >= g && g >= b)) low.push(`dark bg ${await token("--bg")} is cool`);
    }
  }
  await theme("light");
  return low.join("; ");
});

// Complaint 14: Copy link confirms in place with a check icon, no toast, no tooltip
await criterion("copy link confirms in place, no toast", async () => {
  await js(`navigator.clipboard.writeText = () => Promise.resolve()`);
  await clickOn("#copy-link");
  await sleep(200);
  const [said, toast, check] = await js(`[document.querySelector("#said")?.textContent, !document.querySelector("#toast")?.hidden, /M5 12/.test(document.querySelector("#copy-link svg")?.innerHTML ?? "")]`);
  return said === "Link copied" && !toast && check ? "" : `said ${JSON.stringify(said)}, toast ${toast}, check ${check}`;
});

// Calculator, the line: Enter commits to the history, shown open on top, and the line empties
await criterion("Enter commits to the history, open on top, line empty", async () => {
  await js(`document.querySelector("math-field").focus()`);
  await sleep(150);
  await key("Enter");
  await sleep(600);
  const [v, n, open] = await js(`[document.querySelector("math-field")?.value, document.querySelectorAll(".history .hcell").length, document.querySelectorAll(".history .hcell .sheet-body > *").length > 0]`);
  return v === "" && n === 1 && open ? "" : `line ${JSON.stringify(v)}, ${n} entries, open ${open}`;
});

// History: up puts the question back in the line; x removes it with an Undo strip under the bar
await criterion("history: up reloads, x removes with Undo at the top", async () => {
  await clickOn('.history .hcell [aria-label="Put in the line"]');
  await sleep(300);
  const v = await line();
  if (!v || !v.replace(/\s/g, "").includes("2x^2")) return `up gave ${JSON.stringify(v)}`;
  await clickOn('.history .hcell [aria-label="Remove from history"]');
  await sleep(300);
  const [n, undo, top] = await js(`[document.querySelectorAll(".history .hcell").length, !document.querySelector("#undo").hidden, document.querySelector("#undo").getBoundingClientRect().top]`);
  return undo && top < 200 ? "" : `${n} entries, undo shown ${undo} at ${top}`;
});

// The line: Esc Esc clears the line and the sheet whole (complaints 3 and 5)
await criterion("Esc Esc clears the line and the sheet", async () => {
  await write(QUAD);
  await js(`document.querySelector("math-field").focus()`);
  await key("Escape");
  await key("Escape");
  await sleep(400);
  const [v, a] = await js(`[document.querySelector("math-field")?.value, document.querySelector(".card .answer")?.childElementCount]`);
  return v === "" && a === 0 ? "" : `line ${JSON.stringify(v)}, answer has ${a} children`;
});

// On-screen keyboard: fixed overlay, 3 layer tabs, 8 x 4 keys of 48, Esc closes it
await criterion("keyboard: fixed overlay of 3 layers, 8 x 4 keys 48 tall", async () => {
  await clickOn("#open-keys");
  await sleep(400);
  const m = await js(`(() => { const k = document.querySelector("#keys"); if (!k || k.hidden) return ["keyboard not shown"]; const out = [];
    const cs = getComputedStyle(k); if (cs.position !== "fixed" || cs.zIndex !== "30") out.push("position " + cs.position + " z " + cs.zIndex);
    const tabs = k.querySelectorAll(".keys-head [role=tab], .keys-head .layer").length; if (tabs !== 3) out.push(tabs + " layer tabs");
    const keys = [...k.querySelectorAll("button.key")]; if (keys.length !== 32) out.push(keys.length + " keys");
    const h = keys.map((b) => Math.round(b.getBoundingClientRect().height)); if (h.some((x) => x !== 48)) out.push("key heights " + [...new Set(h)].join("/"));
    return out; })()`);
  return m.join("; ");
});

// Complaint 16: a key types its symbol and the focus never leaves the line
await criterion("keyboard key types its symbol, focus stays in the line", async () => {
  await js(`document.querySelector("math-field").focus()`);
  const seven = await js(`(() => { const b = [...document.querySelectorAll("#keys button.key")].find((b) => b.textContent.trim() === "7"); if (!b) return null; const r = b.getBoundingClientRect(); return [r.x + r.width / 2, r.y + r.height / 2]; })()`);
  if (!seven) return "no key 7";
  await click(seven[0], seven[1]);
  await sleep(300);
  const [v, active] = await js(`[document.querySelector("math-field").value, document.activeElement?.tagName]`);
  return v === "7" && active === "MATH-FIELD" ? "" : `line ${JSON.stringify(v)}, active ${active}`;
});

// Five rules, 4: the keyboard adds no padding to the body and does not change the scroll height
await criterion("keyboard open: no body padding, scroll height unchanged", async () => {
  await js(`document.querySelector("#close-keys")?.click()`);
  await sleep(300);
  const before = await js(`[document.documentElement.scrollHeight, getComputedStyle(document.body).paddingBottom]`);
  await clickOn("#open-keys");
  await sleep(400);
  const after = await js(`[document.documentElement.scrollHeight, getComputedStyle(document.body).paddingBottom]`);
  return before[0] === after[0] && before[1] === after[1] ? "" : `scroll ${before[0]} -> ${after[0]}, padding ${before[1]} -> ${after[1]}`;
});
await criterion("no text under 14 px with the keyboard open", async () => (await js(SMALLEST)).map(([e, px]) => `${e} ${px}px`).join(", "));

// On-screen keyboard: Esc closes it, and the keyboard does not fold on a press
await criterion("Esc closes the keyboard", async () => {
  await key("Escape");
  await sleep(300);
  return (await js(`document.querySelector("#keys").hidden`)) ? "" : "keyboard still open after Esc";
});

// Complaint 17: sin, ln, lg and log keys are answered, each by its example
await criterion("sin, ln, lg, log keys answered", async () => {
  await clickOn("#open-keys");
  await sleep(300);
  await js(`[...document.querySelectorAll("#keys .layer")].find((t) => t.textContent === "f(x)")?.click()`);
  const keys = await js(`[...document.querySelectorAll("#keys button.key[data-example]")].filter((b) => /^(sin|ln|lg|log)/.test(b.textContent.trim())).map((b) => [b.textContent.trim(), b.dataset.example])`);
  await key("Escape");
  if (keys.length < 4) return `${keys.length} of sin, ln, lg, log on the keyboard`;
  const refused = [];
  for (const [label, example] of keys) {
    await write(example);
    const a = await answerText();
    if (!a || /refused|not answered|did not load/i.test(a) || !(await js(`Boolean(document.querySelector(".card .result .math"))`))) refused.push(label);
  }
  return refused.length ? `refused: ${refused.join(", ")}` : "";
});

// Five rules, 4: the reference drawer and the settings dialog float, nothing moves the page
await criterion("reference drawer and settings: no body padding, scroll height unchanged", async () => {
  await write(QUAD);
  const h = () => js(`[document.documentElement.scrollHeight, getComputedStyle(document.body).paddingBottom, getComputedStyle(document.body).paddingRight]`);
  const before = await h();
  await clickOn("#open-ref");
  await sleep(400);
  const drawer = await h();
  await key("Escape");
  await sleep(200);
  await clickOn("#open-settings");
  await sleep(300);
  const dialog = await h();
  await key("Escape");
  await sleep(200);
  const same = (a) => a[0] === before[0] && a[1] === before[1] && a[2] === before[2];
  return same(drawer) && same(dialog) ? "" : `before ${before}, drawer ${drawer}, settings ${dialog}`;
});

// Settings: the rows of the dialog, the keyboard with Auto, and no digits setting anywhere
await criterion("settings rows: theme, text size, keyboard auto/on/off, formulas, saved work; no digits", async () => {
  const m = await js(`(() => { const legends = [...document.querySelectorAll("#settings legend")].map((l) => l.textContent.trim().toLowerCase()); const out = [];
    for (const want of ["theme", "text size", "keyboard", "formulas", "saved work"]) if (!legends.some((l) => l.includes(want))) out.push("no " + want);
    if (legends.some((l) => /digit|precision/.test(l))) out.push("a digits setting");
    const kb = [...document.querySelectorAll('#settings [data-setting="keyboard"] button')].map((b) => b.value); if (kb.join() !== "auto,on,off") out.push("keyboard values " + kb.join("/"));
    return out; })()`);
  return m.join("; ");
});

// Complaint 6: a slow answer shows that the page is working, and the result still comes
await criterion("slow answer: busy state, then Still working with Stop", async () => {
  await js(`window.__hold = [12000]`);
  await write("x^2-2=0", { settle: 600 });
  const busy = await js(`document.querySelector(".card.pending") !== null`);
  await sleep(10200);
  const still = await js(`/still working/i.test(document.querySelector(".card .status")?.textContent ?? "") && Boolean(document.querySelector(".card .status button"))`);
  await sleep(2500);
  const done = await js(`!document.querySelector(".card.pending") && Boolean(document.querySelector(".card .result .math"))`);
  return busy && still && done ? "" : `busy ${busy}, still working ${still}, answered ${done}`;
});

// Complaint 6: of two questions in flight, the newest wins
await criterion("newest request wins", async () => {
  await js(`window.__hold = [3000]`);
  await write("3*3", { settle: 300 });
  await write("2*2", { settle: 5000 });
  const r = await js(`document.querySelector(".card .result .math")?.textContent ?? ""`);
  return r.includes("4") && !r.includes("9") ? "" : `result reads ${JSON.stringify(r)}`;
});

// Acceptance shots: the quadratic, a refusal, the inequality with its axis, the keyboard, at both
// widths and in both themes; the phone width measured for horizontal scroll and small text
const scenes = [["quadratic", QUAD], ["refusal", "\\sqrt{"], ["inequality", "x^2-4>0"]];
if (thrown.length === 0) for (const width of [1280, 400]) {
  await viewport(width);
  for (const which of ["light", "dark"]) {
    await theme(which);
    await open();
    if (width === 400 && which === "light") {
      await criterion("no horizontal page scroll at 400", async () => {
        const [sw, cw] = await js(`[document.documentElement.scrollWidth, document.documentElement.clientWidth]`);
        return sw <= cw ? "" : `scrollWidth ${sw} > ${cw}`;
      });
    }
    await shot(`empty-${width}-${which}`);
    for (const [name, latex] of scenes) {
      await write(latex);
      if (name === "inequality") await js(`[...document.querySelectorAll(".card details.sec")].find((d) => /plot/i.test(d.querySelector("summary")?.textContent ?? ""))?.setAttribute("open", "")`);
      await sleep(400);
      await shot(`${name}-${width}-${which}`);
      if (width === 400 && which === "light" && name === "quadratic") {
        await criterion("no horizontal page scroll at 400, answered", async () => {
          const [sw, cw] = await js(`[document.documentElement.scrollWidth, document.documentElement.clientWidth]`);
          return sw <= cw ? "" : `scrollWidth ${sw} > ${cw}`;
        });
        await criterion("no text under 14 px at 400", async () => (await js(SMALLEST)).map(([e, px]) => `${e} ${px}px`).join(", "));
      }
    }
    await write(QUAD);
    await clickOn("#open-keys");
    await sleep(400);
    await shot(`keyboard-${width}-${which}`);
    await key("Escape");
  }
}

// Complaint 18: the gallery is a page of its own, linked from the bar, loading no core
await criterion("gallery is its own page", async () => {
  const href = await js(`document.querySelector('.bar a[href="gallery/"]')?.getAttribute("href") ?? null`);
  if (!href) return "no link to gallery/ in the bar";
  await send("Page.navigate", { url: origin + "gallery/" });
  await sleep(2500);
  const [tiles, worker] = await js(`[document.querySelectorAll(".tile").length, document.querySelector('script[src*="core"]') !== null]`);
  return tiles > 0 && !worker ? "" : `${tiles} tiles, core script ${worker}`;
});

const passed = results.filter(([, ok]) => ok).length;
console.log(`criteria: ${passed} of ${results.length} pass`);
ws.close();
const gone = new Promise((r) => chrome.on("exit", r));
chrome.kill();
await gone;
server.close();
rmSync(profile, { recursive: true, force: true, maxRetries: 10, retryDelay: 300 });
process.exit(passed === results.length ? 0 : 1);
