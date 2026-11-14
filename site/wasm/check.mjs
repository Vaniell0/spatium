// The core in WebAssembly answers as the host does: each line of questions.txt is asked of both,
// and the answers must agree -- the same JSON, numbers to 1e-9 of their size (the WebAssembly
// build has its own libm, musl's, whose last digits differ from glibc's), numbers inside a string
// likewise and the rest of the string exactly, a row of numbers at the scale of its largest, and
// `computed` and a route's `cost` (a candidate's, and the cost column of a sheet's table) the only
// fields allowed to differ: a root found to the last bit
// takes as many evaluations as the libm's last digits lead it to. A graph's numbers are its fp32
// twin, the device's, and are compared to four float steps of their scale (2^-22): musl's tanhf
// and glibc's differ in a float's last bit.
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
const rows = core.cwrap('spatium_rows', 'string', ['string']);
const evaluate = core.cwrap('spatium_evaluate', 'string', ['string']);
const plot = core.cwrap('spatium_plot', 'string', ['string']);
const zeta = core.cwrap('spatium_zeta', 'string', ['string']);

const number = /-?\d+(\.\d+)?(e[+-]?\d+)?/gi;
let tolerance = 1e-9;   // relative; a graph's is a float's, set per line below
const close = (a, b) => a === b || Math.abs(a - b) <= tolerance * Math.max(Math.abs(a), Math.abs(b));
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
        if (Math.abs(a[i] - b[i]) > tolerance * scale) return `${where}[${i}]: ${a[i]} vs ${b[i]} (scale ${scale})`;
      return true;
    }
    for (let i = 0; i < a.length; ++i) { const r = agree(a[i], b[i], `${where}[${i}]`); if (r !== true) return r; }
    return true;
  }
  if (a && b && typeof a === 'object' && typeof b === 'object') {
    // a sheet's table of routes: the cost column counts evaluations, as a candidate's cost does
    if (Array.isArray(a.columns) && Array.isArray(a.rows) && a.columns.includes('cost')) {
      const skip = a.columns.indexOf('cost');
      const drop = t => ({ columns: t.columns.filter((_, i) => i !== skip), rows: (t.rows ?? []).map(r => r.filter((_, i) => i !== skip)) });
      return agree(drop(a), drop(b), where);
    }
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
  const wasm = line.startsWith('rows ') ? rows(line.slice(5))
             : line.startsWith('evaluate ') ? evaluate(line.slice(9))
             : line.startsWith('plot ') ? plot(line.slice(5))
             : line.startsWith('zeta ') ? zeta(line.slice(5)) : respond(line, 0);
  tolerance = line.startsWith('plot ') ? 2 ** -22 : 1e-9;
  const r = agree(JSON.parse(host[i]), JSON.parse(wasm), line);
  if (r !== true) { ++failed; console.log(`differs: ${r}`); }
  if ((line.startsWith('evaluate ') || line.startsWith('rows ')) && !wasm.includes('"computed":"wasm"') && wasm.includes('"certificate"'))
    { ++failed; console.log(`differs: ${line}: the certificate does not say wasm`); }
});
console.log(`${lines.length - failed} of ${lines.length} agree`);

// A card's question is the rows the library writes (`question_latex`); read back as the page reads
// them, by the formula editor's parser, they ask the card's question again: the same rows, the same
// route and the same answer.
const ls = await import(path.join(here, '../vendor/latex-syntax/latex-syntax.min.js'));
const gallery = JSON.parse(readFileSync(path.join(here, '../gallery/gallery.json'), 'utf8'));
const reference = JSON.parse(readFileSync(path.join(here, '../data/reference.json'), 'utf8'));
// every example of the reference the library answered with rows: a question, a function, an operator
const examples = [...reference.questions.map(q => [q.form, q]), ...reference.functions.map(f => [f.name, f]),
                  ...reference.operators.map(o => [o.written, o])];
const cards = [...(gallery.cards ?? gallery).filter(c => c.question_latex),
               ...examples.filter(([, e]) => e.answer?.question_latex)
                 .map(([id, e]) => ({ id, question_latex: e.answer.question_latex, route: e.answer.route, refused: e.answer.refused, example: true }))];
let unread = 0;
for (const c of cards) {
  const a = JSON.parse(rows(JSON.stringify({ rows: c.question_latex.map(r => ls.parse(r)) })));
  const pieces = x => JSON.stringify((x.answer ?? []).map(p => p.latex));
  // a reference example keeps its route and refusal; a card the whole answer
  const keys = c.example ? ['route', 'question_latex', 'refused'] : ['route', 'question_latex', 'latex', 'values', 'value', 'refused'];
  const bad = keys.filter(k => JSON.stringify(a[k]) !== JSON.stringify(c[k]));
  if (!c.example && pieces(a) !== pieces(c)) bad.push('answer');
  if (bad.length) { ++unread; console.log(`differs: card ${c.id} read back from its rows: ${bad.join(', ')}`); }
}
console.log(`${cards.length - unread} of ${cards.length} questions of the gallery and the reference read back from their rows`);

// The reference's own notation: what is typed in the editor (Writing, the rows the library wrote)
// and a function as the editor writes it (`\lg`, `\log_3`, `\sqrt[3]`, rows written by hand), read
// back by the editor's parser, answer as the line of the language beside them does.
const notes = [...reference.writing.map(w => ({ id: w.typed, rows: w.answer.question_latex, answer: w.answer })),
               ...reference.functions.filter(f => f.rows).map(f => ({ id: f.head, rows: f.rows, answer: f.answer }))];
let misread = 0;
for (const n of notes) {
  const a = JSON.parse(rows(JSON.stringify({ rows: n.rows.map(r => ls.parse(r)) })));
  const bad = ['route', 'refused'].filter(k => JSON.stringify(a[k]) !== JSON.stringify(n.answer[k]));
  if (n.answer.value !== undefined && !(Math.abs(a.value - Number(n.answer.value)) <= 1e-9 * Math.abs(a.value))) bad.push(`value ${a.value}`);
  if (n.answer.label !== undefined && a.sheet?.label !== n.answer.label) bad.push(`label ${a.sheet?.label}`);
  if (bad.length) { ++misread; console.log(`differs: reference ${n.id} read from ${JSON.stringify(n.rows)}: ${bad.join(', ')}`); }
}
console.log(`${notes.length - misread} of ${notes.length} rows of the reference's notation read as the language they stand for`);
// What a sheet offers next (`Next`) is put in the line as the row the library wrote for it: read
// back by the editor's parser, that row asks what the offer's text asks -- the same route and the
// same answer.
let offers = 0, wrong = 0;
for (const line of lines.filter(l => !/^(rows|evaluate|plot) /.test(l))) {
  const sheet = JSON.parse(respond(line, 0)).sheet;
  for (const n of (sheet?.sections ?? []).flatMap(s => s.next ?? [])) {
    ++offers;
    const byText = JSON.parse(respond(n.text, 0)), byRow = JSON.parse(rows(JSON.stringify({ rows: [ls.parse(n.latex)] })));
    const bad = ['route', 'refused'].filter(k => JSON.stringify(byRow[k]) !== JSON.stringify(byText[k]));
    if (byRow.sheet?.result !== byText.sheet?.result) bad.push(`result ${byRow.sheet?.result} vs ${byText.sheet?.result}`);
    if (bad.length) { ++wrong; console.log(`differs: next of "${line}": ${n.text} read from ${n.latex}: ${bad.join(', ')}`); }
  }
}
console.log(`${offers - wrong} of ${offers} questions the sheets offer next read back from their rows`);
process.exit(failed || unread || misread || wrong ? 1 : 0);
