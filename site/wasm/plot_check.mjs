// The graph's WebGL2 shader against the core, in a headless browser: serves site/wasm, opens
// plot_check.html over the DevTools protocol, and reads its report. WebGL2 runs on the browser's
// software renderer (SwiftShader), so no GPU is needed.
//
//   node site/wasm/plot_check.mjs [browser]    (default: the first of chromium, google-chrome found)
//
// It fails when the page cannot run, when any pixel's F differs from the core's fp32 value by more
// than the tolerance below, in pixels of height -- how far the drawn curve would move -- or when a
// case was compiled more than once: moving a slot must not compile a shader. It must also tell a
// wrong frame apart: one drawn with a slot a hundredth off its twin's has to exceed the tolerance.
// Measured on SwiftShader: 0.006 px for sin, 0.004 px for cos, sqrt and log at arguments up to 27,
// 8e-6 px for exp.
import { spawn, execFileSync } from 'node:child_process';
import { createServer } from 'node:http';
import { readFileSync, writeFileSync, mkdtempSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

const tolerance = 1 / 16;   // of a pixel
const here = path.dirname(fileURLToPath(import.meta.url));
const types = { '.html': 'text/html', '.js': 'text/javascript', '.wasm': 'application/wasm' };
const server = createServer((req, res) => {
  const file = path.join(here, path.basename(new URL(req.url, 'http://x').pathname));
  let body;
  try { body = readFileSync(file); } catch { res.writeHead(404); res.end(); return; }
  res.writeHead(200, { 'content-type': types[path.extname(file)] ?? 'application/octet-stream' });
  res.end(body);
});
await new Promise((ok) => server.listen(0, '127.0.0.1', ok));
const url = `http://127.0.0.1:${server.address().port}/plot_check.html`;

const which = (name) => { try { return execFileSync('which', [name], { stdio: ['ignore', 'pipe', 'ignore'] }).toString().trim(); } catch { return null; } };
const browser = process.argv[2] ?? ['chromium', 'chromium-browser', 'google-chrome', 'google-chrome-stable'].map(which).find(Boolean);
if (!browser) { console.log('plot check: no browser found'); process.exit(1); }
const profile = mkdtempSync(path.join(tmpdir(), 'plot-check-'));
const port = 9300 + Math.floor(Math.random() * 600);
const chrome = spawn(browser, ['--headless=new', '--no-first-run', '--no-sandbox', `--user-data-dir=${profile}`,
  '--use-angle=swiftshader', '--enable-unsafe-swiftshader', `--remote-debugging-port=${port}`, 'about:blank'], { stdio: 'ignore' });
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const done = (code) => {
  chrome.once('exit', () => { server.close(); rmSync(profile, { recursive: true, force: true }); process.exit(code); });
  chrome.kill();
};

let ws;
for (let i = 0; i < 100 && !ws; i++) {
  try {
    const page = (await (await fetch(`http://127.0.0.1:${port}/json`)).json()).find((t) => t.type === 'page');
    const socket = new WebSocket(page.webSocketDebuggerUrl);
    await new Promise((ok, bad) => { socket.onopen = ok; socket.onerror = bad; });
    ws = socket;
  } catch { await sleep(200); }
}
if (!ws) { console.log('plot check: the browser did not start'); done(1); }
let id = 0;
const waiting = new Map();
ws.onmessage = (m) => { const d = JSON.parse(m.data); if (d.id && waiting.has(d.id)) { waiting.get(d.id)(d); waiting.delete(d.id); } };
const send = (method, params = {}) => new Promise((ok) => { const n = ++id; waiting.set(n, ok); ws.send(JSON.stringify({ id: n, method, params })); });

await send('Page.enable');
await send('Page.navigate', { url });
let report;
for (let i = 0; i < 150 && !report; i++) {
  await sleep(200);
  const r = await send('Runtime.evaluate', { expression: 'window.result', awaitPromise: true, returnByValue: true });
  report = r.result?.result?.value;
}
if (!report) { console.log('plot check: the page gave no report'); done(1); }
// PLOT_PNG=<path> keeps the drawn frame, for a person to look at
if (report.png && process.env.PLOT_PNG) writeFileSync(process.env.PLOT_PNG, Buffer.from(report.png.split(',')[1], 'base64'));
delete report.png;
console.log(JSON.stringify(report));
const ok = report.ok && report.worst <= tolerance && report.compiled === report.cases.length && report.wrong > tolerance;
console.log(ok ? `plot check: ${report.runs} frames of ${report.compiled} shaders agree, worst ${report.worst.toExponential(2)} px`
               : `plot check: failed${report.why ? ': ' + report.why : ''}`);
done(ok ? 0 : 1);
