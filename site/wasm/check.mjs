// The core in WebAssembly answers as the host does: each line of questions.txt is asked of both,
// and the answers must agree -- the same JSON, numbers to 1e-9 of their size (the WebAssembly
// build has its own libm, musl's, whose last digits differ from glibc's), numbers inside a string
// likewise and the rest of the string exactly, a row of numbers at the scale of its largest, and
// `computed` and a route's `cost` the only fields allowed to differ: a root found to the last bit
// takes as many evaluations as the libm's last digits lead it to.
//
//   node site/wasm/check.mjs <host binary>
import { readFileSync } from 'node:fs';
import { execFileSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

const here = path.dirname(fileURLToPath(import.meta.url));
const lines = readFileSync(path.join(here, 'questions.txt'), 'utf8').split('\n').filter(l => l.length > 0);
const host = execFileSync(process.argv[2], { input: lines.join('\n') + '\n' }).toString().trim().split('\n');
const SpatiumCore = (await import(path.join(here, 'spatium_core.js'))).default;
const core = await SpatiumCore();
const respond = core.cwrap('spatium_respond', 'string', ['string', 'number']);
const preview = core.cwrap('spatium_preview', 'string', ['string']);
const evaluate = core.cwrap('spatium_evaluate', 'string', ['string']);

const number = /-?\d+(\.\d+)?(e[+-]?\d+)?/gi;
const close = (a, b) => a === b || Math.abs(a - b) <= 1e-9 * Math.max(Math.abs(a), Math.abs(b));
function agree(a, b, where) {
  if (typeof a === 'number' && typeof b === 'number') return close(a, b) || `${where}: ${a} vs ${b}`;
  if (typeof a === 'string' && typeof b === 'string') {
    if (a === b) return true;
    const na = a.match(number) ?? [], nb = b.match(number) ?? [];
    if (a.replace(number, '#') !== b.replace(number, '#') || na.length !== nb.length) return `${where}: "${a}" vs "${b}"`;
    for (let i = 0; i < na.length; ++i) if (!close(+na[i], +nb[i])) return `${where}: "${a}" vs "${b}"`;
    return true;
  }
  if (Array.isArray(a) && Array.isArray(b)) {
    if (a.length !== b.length) return `${where}: ${a.length} items vs ${b.length}`;
    // a row of numbers -- the coefficients of one polynomial -- is compared at the scale of its
    // largest, since the small ones are the libm's noise in the same absolute units
    if (a.length > 0 && a.every(x => typeof x === 'number') && b.every(x => typeof x === 'number')) {
      const scale = Math.max(...a.map(Math.abs), ...b.map(Math.abs));
      for (let i = 0; i < a.length; ++i)
        if (Math.abs(a[i] - b[i]) > 1e-9 * scale) return `${where}[${i}]: ${a[i]} vs ${b[i]} (scale ${scale})`;
      return true;
    }
    for (let i = 0; i < a.length; ++i) { const r = agree(a[i], b[i], `${where}[${i}]`); if (r !== true) return r; }
    return true;
  }
  if (a && b && typeof a === 'object' && typeof b === 'object') {
    const keys = Object.keys(a);
    if (keys.join() !== Object.keys(b).join()) return `${where}: keys ${keys} vs ${Object.keys(b)}`;
    for (const k of keys) {
      if (k === 'computed' || k === 'cost') continue;
      const r = agree(a[k], b[k], `${where}.${k}`);
      if (r !== true) return r;
    }
    return true;
  }
  return a === b || `${where}: ${JSON.stringify(a)} vs ${JSON.stringify(b)}`;
}

let failed = 0;
lines.forEach((line, i) => {
  const wasm = line.startsWith('preview ') ? preview(line.slice(8))
             : line.startsWith('evaluate ') ? evaluate(line.slice(9)) : respond(line, 0);
  const r = agree(JSON.parse(host[i]), JSON.parse(wasm), line);
  if (r !== true) { ++failed; console.log(`differs: ${r}`); }
  if ((line.includes(' for ') || line.startsWith('evaluate ')) && !line.startsWith('preview ') && !wasm.includes('"computed":"wasm"') && wasm.includes('"certificate"'))
    { ++failed; console.log(`differs: ${line}: the certificate does not say wasm`); }
});
console.log(`${lines.length - failed} of ${lines.length} agree`);
process.exit(failed ? 1 : 0);
