// The canvas: an infinite field of nodes, each a cell of the document written in the language,
// wired by the names they read, placed in layers by dependency, with Layers on the left (the
// Layout tree), a toolbar in the bar and the Inspector on the right. A formula node shows its line
// rendered as the core read it and is edited as text; a graph node draws the core's shader. What
// a node says comes from its answer; the page composes nothing.

import type { Answer } from "../core.js";
import { layout, nodesOf, wirePath, wires, type Node, type Placed } from "../canvas.js";
import { type Cell, type Doc, graph, math } from "../doc.js";
import { answerView, type Sheet } from "./answer.js";
import { $, el, icon, math as tex } from "./dom.js";
import { graphView, type GraphOptions } from "./graph.js";

export interface CanvasHost {
  source: () => boolean;
  changed(): void;
  curve: (cell: Cell) => GraphOptions["curve"];
  ask: (text: string, rows?: string[]) => void;
  undo(message: string, restore: () => void): void;
}

interface View { x: number; y: number; z: number }
type Tool = "select" | "hand";
const VIEW = "spatium.canvas.view";
const GLYPH = { formula: "fx", graph: "∿" } as const;
// every kind the design names; the ones not yet on the page are listed, so the menu is the plan
const KINDS: [string, string, boolean][] = [["Formula", "a line of the language", true], ["Graph", "the plot of a name", true],
  ["Solve", "routes and the result of an equation", false], ["Parameters", "names on sliders", false], ["Scene", "shapes in the language", false],
  ["View", "the scene rendered", false], ["Sound", "nodes chained into a wave, out as .wav", false]];

export class Canvas {
  readonly root = el("div", { class: "canvas", hidden: true });
  readonly tools = el("div", { class: "tools", role: "toolbar", "aria-label": "Canvas tools", hidden: true });
  private readonly layers = el("aside", { class: "panel layers", "aria-label": "Layers" });
  private readonly field = el("div", { class: "field", tabindex: "0", "aria-label": "Field" });
  private readonly sheet = el("div", { class: "sheet" });   // the transformed plane
  private readonly wiresSvg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  private readonly inspector = el("aside", { class: "panel inspector", "aria-label": "Inspector" });
  private readonly zoom = el("button", { type: "button", class: "zoom mono", "aria-label": "Zoom" });
  private nodes: Node[] = [];
  private placed: Placed[] = [];
  private boxes = new Map<string, HTMLElement>();
  private answers = new Map<Cell, Answer | null>();
  private above = new Map<Cell, number>();   // rows asked before the cell's own, which its answer's latex includes
  private drawn = new Map<Cell, Answer>();   // a graph cell's plot, as drawn
  private selected: string | null = null;
  private tool: Tool = "select";
  private view: View = { x: 0, y: 0, z: 1 };
  private editing: string | null = null;
  private slow: ReturnType<typeof setTimeout> | 0 = 0;

  constructor(private doc: Doc, private readonly host: CanvasHost) {
    try { this.view = { ...this.view, ...JSON.parse(localStorage.getItem(VIEW) ?? "{}") }; } catch { /* the default view */ }
    this.wiresSvg.setAttribute("class", "wires");
    this.sheet.append(this.wiresSvg);
    this.field.append(this.sheet);
    this.root.append(this.layers, this.field, this.inspector);
    this.toolbar();
    this.gestures();
    this.build();
  }

  get document() { return this.doc; }
  setDoc(doc: Doc) { this.doc = doc; this.selected = null; this.build(); requestAnimationFrame(() => this.fit()); }
  answerOf(cell: Cell): Answer | null { return this.answers.get(cell) ?? null; }

  // every node drawn again from the document, then placed by the layout once its height is known
  build(select?: string) {
    this.nodes = nodesOf(this.doc);
    for (const cell of [...this.answers.keys()]) if (!this.doc.cells.includes(cell)) this.answers.delete(cell);
    this.boxes.clear();
    for (const n of this.nodes) this.boxes.set(n.id, this.node(n));
    this.sheet.replaceChildren(this.wiresSvg, ...this.boxes.values());
    if (select !== undefined) this.selected = select;
    if (this.selected && !this.boxes.has(this.selected)) this.selected = null;
    this.place();
    this.tree();
    this.select(this.selected);
    this.apply();
  }

  // an answer for a cell: a formula node shows its line as read, its result and its strength; a
  // graph node its plot; a graph the page derives from a formula is not a node of its own
  show(cell: Cell, answer: Answer | null, above = 0, from?: Cell) {
    if (from) return;
    this.answers.set(cell, answer);
    this.above.set(cell, above);
    const n = this.nodes.find((x) => x.cell === cell);
    if (!n) return;
    const box = this.boxes.get(n.id);
    if (!box) return;
    box.classList.remove("pending");
    if (cell.kind === "graph" && answer) this.drawn.set(cell, answer);
    this.fill(n, box);
    this.place();
    if (this.selected === n.id) this.inspect();
  }

  pending(stop: () => void) {
    for (const box of this.boxes.values()) box.classList.add("pending");
    clearTimeout(this.slow);
    this.slow = setTimeout(() => {
      if (![...this.boxes.values()].some((b) => b.classList.contains("pending"))) return;
      this.field.append(el("p", { class: "status float", role: "status" }, "Still working",
        el("button", { type: "button", class: "quiet small", text: "Stop", onclick: () => { stop(); this.field.querySelector(".status")?.remove(); } })));
    }, 10000);
  }
  clearDerived() { this.field.querySelector(".status")?.remove(); }

  // cells added, each placed by the layout; the first new one selected
  append(cells: Cell[]) {
    const last = this.doc.cells.at(-1);
    if (last?.kind === "math" && !("rows" in last.written ? last.written.rows.join("") : last.written.text).trim()) this.doc.cells.pop();
    this.doc.cells.push(...cells);
    this.build(nodesOf(this.doc)[this.doc.cells.indexOf(cells[0])]?.id);
    this.host.changed();
    this.reveal(this.selected);
    if (cells[0].kind === "math") this.edit(this.selected!);
  }

  // ── nodes ───────────────────────────────────────────────────────────

  private node(n: Node): HTMLElement {
    const box = el("article", { class: `node ${n.kind}`, "data-id": n.id, tabindex: "0", style: `width: ${n.w}px` });
    const head = el("header", { class: "head" }, el("span", { class: "id", text: n.id }), el("span", { class: "kind-note" }));
    box.append(head, el("div", { class: "nbody" }));
    if (n.reads.length) box.append(el("span", { class: "pin in", "aria-hidden": "true" }));
    if (n.defines.length || n.kind === "formula") box.append(el("span", { class: "pin out", "aria-hidden": "true" }));
    box.addEventListener("pointerdown", (e) => {
      if (this.tool === "hand" || this.editing === n.id) return;
      this.select(n.id, e.shiftKey);
      if ((e.target as Element).closest(".head")) this.drag(n, box, e);
    });
    box.addEventListener("dblclick", () => { if (n.kind === "formula") this.edit(n.id); });
    box.addEventListener("keydown", (e) => { if (e.key === "Enter" && this.editing !== n.id && n.kind === "formula") { e.preventDefault(); this.edit(n.id); } });
    this.fill(n, box);
    return box;
  }

  private fill(n: Node, box: HTMLElement) {
    const body = box.querySelector(".nbody")!;
    const note = box.querySelector(".kind-note")!;
    const a = this.answers.get(n.cell);
    if (n.kind === "graph") {
      note.textContent = `along ${n.cell.kind === "graph" ? n.cell.of : ""}`;
      const plot = this.drawn.get(n.cell);
      if (!plot) { body.replaceChildren(el("p", { class: "muted", text: a?.refused ? String(a.refused) : "waits for its name" })); return; }
      if (!body.querySelector(".graph")) body.replaceChildren(graphView(plot, { height: 152, curve: this.host.curve(n.cell) }));
      return;
    }
    note.textContent = "DSL";
    if (this.editing === n.id) return;
    const s = a?.sheet as Sheet | undefined;
    const read = (a?.question_latex as string[] | undefined)?.slice(this.above.get(n.cell) ?? 0);
    const line = !n.text.trim() ? el("p", { class: "muted", text: "name = expr" })
      : read?.length ? el("div", { class: "line" }, ...read.map((r) => tex(r, { cls: "math", source: this.host.source() })))
      : el("code", { class: "src mono", text: n.text });
    const parts: (HTMLElement | null)[] = [line];
    if (a?.refused && !s) parts.push(el("p", { class: "told", text: String(a.refused) }));
    else if (s?.outcome === "refused") parts.push(el("p", { class: "told", text: String(a?.refused ?? s.label) }));
    else if (s?.result) parts.push(tex(s.result, { cls: "math res", display: true, source: this.host.source() }),
      el("p", { class: "gist", text: [s.verdict?.strength ?? s.label, ...(s.verdict?.notes ?? []).slice(0, 1)].join(" · ") }));
    body.replaceChildren(...parts.filter((p): p is HTMLElement => p !== null));
  }

  // the node's text as a field in its body: Enter commits, Shift+Enter a new line, Esc cancels
  private edit(id: string) {
    const n = this.nodes.find((x) => x.id === id);
    const box = this.boxes.get(id);
    if (!n || !box || n.cell.kind !== "math") return;
    this.editing = id;
    const cell = n.cell;
    const t = el("textarea", { class: "text mono", rows: String(Math.max(1, n.text.split("\n").length)), spellcheck: "false", "aria-label": `${id}: source` });
    t.value = n.text;
    let done = false;
    const finish = (commit: boolean) => {
      if (done) return;
      done = true;
      this.editing = null;
      if (commit && t.value !== n.text) {
        cell.written = { text: t.value };
        this.build(id);
        this.host.changed();
      } else this.fill(n, box);
      this.place();
      box.focus();
    };
    t.addEventListener("keydown", (e) => {
      if (e.key === "Enter" && !e.shiftKey) { e.preventDefault(); e.stopPropagation(); finish(true); }
      else if (e.key === "Escape") { e.preventDefault(); e.stopPropagation(); finish(false); }
    });
    t.addEventListener("blur", () => finish(true));
    box.querySelector(".nbody")!.replaceChildren(t);
    this.place();
    t.focus();
    t.setSelectionRange(t.value.length, t.value.length);
  }

  private select(id: string | null, add = false) {
    void add;
    this.selected = id;
    for (const [k, b] of this.boxes) b.classList.toggle("selected", k === id);
    for (const r of this.layers.querySelectorAll<HTMLElement>(".row")) r.setAttribute("aria-selected", String(r.dataset.id === id));
    this.inspect();
    this.wire();
  }

  private remove(id: string) {
    const n = this.nodes.find((x) => x.id === id);
    if (!n) return;
    const at = this.doc.cells.indexOf(n.cell);
    const was = [...this.doc.cells];
    this.doc.cells.splice(at, 1);
    this.build();
    this.host.changed();
    this.host.undo(`${id} removed`, () => { this.doc.cells.splice(0, this.doc.cells.length, ...was); this.build(id); this.host.changed(); });
  }

  // a header dragged moves the node, snapped to 8, and pins it
  private drag(n: Node, box: HTMLElement, e: PointerEvent) {
    const p = this.placed.find((x) => x.id === n.id);
    if (!p) return;
    const start = { x: e.clientX, y: e.clientY, px: p.x, py: p.y };
    let moved = false;
    const move = (ev: PointerEvent) => {
      const dx = (ev.clientX - start.x) / this.view.z, dy = (ev.clientY - start.y) / this.view.z;
      if (!moved && Math.hypot(dx, dy) < 4) return;
      moved = true;
      p.x = Math.round((start.px + dx) / 8) * 8;
      p.y = Math.round((start.py + dy) / 8) * 8;
      box.style.transform = `translate(${p.x}px, ${p.y}px)`;
      this.wire();
    };
    const up = () => {
      box.removeEventListener("pointermove", move);
      box.removeEventListener("pointerup", up);
      box.removeEventListener("pointercancel", up);
      box.releasePointerCapture(e.pointerId);
      if (!moved) return;
      (n.cell.settings as { at?: [number, number] }).at = [p.x, p.y];
      n.at = [p.x, p.y];
      this.place();
      this.inspect();
      this.host.changed();
    };
    box.setPointerCapture(e.pointerId);
    box.addEventListener("pointermove", move);
    box.addEventListener("pointerup", up);
    box.addEventListener("pointercancel", up);
    e.preventDefault();
  }

  // ── placement and wires ─────────────────────────────────────────────

  private place() {
    this.placed = layout(this.nodes, (id) => this.boxes.get(id)?.offsetHeight || 96);
    for (const p of this.placed) {
      const box = this.boxes.get(p.id);
      if (box) box.style.transform = `translate(${p.x}px, ${p.y}px)`;
    }
    this.wire();
  }

  private wire() {
    const at = new Map(this.placed.map((p) => [p.id, p]));
    this.wiresSvg.replaceChildren(...wires(this.nodes).map(([from, to]) => {
      const a = at.get(from), b = at.get(to);
      if (!a || !b) return null;
      const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
      path.setAttribute("d", wirePath(a.x + a.w, a.y + 20, b.x, b.y + 20));
      path.setAttribute("class", from === this.selected || to === this.selected ? "wire lit" : "wire");
      return path;
    }).filter((p): p is SVGPathElement => p !== null));
  }

  // Tidy: positions pinned by hand are released, and the layout lays everything out again
  private tidy() {
    for (const n of this.nodes) { delete (n.cell.settings as { at?: unknown }).at; delete n.at; }
    this.place();
    this.inspect();
    this.host.changed();
  }

  // ── the field: pan, zoom, fit ───────────────────────────────────────

  private apply() {
    this.sheet.style.transform = `translate(${this.view.x}px, ${this.view.y}px) scale(${this.view.z})`;
    this.field.style.backgroundPosition = `${this.view.x}px ${this.view.y}px`;
    this.field.style.backgroundSize = `${32 * this.view.z}px ${32 * this.view.z}px`;
    this.zoom.textContent = `${Math.round(this.view.z * 100)}%`;
    try { localStorage.setItem(VIEW, JSON.stringify(this.view)); } catch { /* not kept */ }
  }

  private zoomAt(z: number, cx: number, cy: number) {
    const nz = Math.min(2, Math.max(.25, z));
    const r = this.field.getBoundingClientRect(), px = cx - r.left, py = cy - r.top;
    this.view.x = px - ((px - this.view.x) * nz) / this.view.z;
    this.view.y = py - ((py - this.view.y) * nz) / this.view.z;
    this.view.z = nz;
    this.apply();
  }

  fit() {
    if (!this.placed.length) { this.view = { x: 0, y: 0, z: 1 }; this.apply(); return; }
    const x0 = Math.min(...this.placed.map((p) => p.x)), y0 = Math.min(...this.placed.map((p) => p.y));
    const x1 = Math.max(...this.placed.map((p) => p.x + p.w)), y1 = Math.max(...this.placed.map((p) => p.y + p.h));
    const r = this.field.getBoundingClientRect();
    const z = Math.min(1, Math.max(.25, Math.min((r.width - 64) / (x1 - x0), (r.height - 64) / (y1 - y0))));
    this.view = { z, x: (r.width - (x1 - x0) * z) / 2 - x0 * z, y: (r.height - (y1 - y0) * z) / 2 - y0 * z };
    this.apply();
  }

  // the node brought into view, the rest of the view unmoved when it already is
  private reveal(id: string | null) {
    const p = this.placed.find((x) => x.id === id);
    if (!p) return;
    const r = this.field.getBoundingClientRect(), z = this.view.z;
    const left = p.x * z + this.view.x, top = p.y * z + this.view.y;
    if (left < 0 || top < 0 || left + p.w * z > r.width || top + p.h * z > r.height) {
      this.view.x = r.width / 2 - (p.x + p.w / 2) * z;
      this.view.y = r.height / 2 - (p.y + p.h / 2) * z;
      this.apply();
    }
  }

  private gestures() {
    this.field.addEventListener("wheel", (e) => {
      e.preventDefault();
      if (e.ctrlKey || e.metaKey) this.zoomAt(this.view.z * Math.exp(-e.deltaY / 400), e.clientX, e.clientY);
      else { this.view.x -= e.deltaX; this.view.y -= e.deltaY; this.apply(); }
    }, { passive: false });
    this.field.addEventListener("pointerdown", (e) => {
      const onNode = (e.target as Element).closest(".node");
      if (onNode && this.tool !== "hand") return;
      if (!onNode && this.tool === "select" && e.button === 0) this.select(null);
      const start = { x: e.clientX - this.view.x, y: e.clientY - this.view.y };
      const move = (ev: PointerEvent) => { this.view.x = ev.clientX - start.x; this.view.y = ev.clientY - start.y; this.apply(); };
      const up = () => { this.field.removeEventListener("pointermove", move); this.field.removeEventListener("pointerup", up); this.field.classList.remove("panning"); };
      this.field.classList.add("panning");
      this.field.setPointerCapture(e.pointerId);
      this.field.addEventListener("pointermove", move);
      this.field.addEventListener("pointerup", up);
      e.preventDefault();
    });
    // shortcuts, whenever no text field has the focus and the canvas is shown
    document.addEventListener("keydown", (e) => {
      if (this.root.hidden || e.ctrlKey || e.metaKey || e.altKey) return;
      const a = document.activeElement;
      if (a instanceof HTMLInputElement || a instanceof HTMLTextAreaElement || a?.tagName === "MATH-FIELD" || (a as HTMLElement | null)?.isContentEditable) return;
      const k = e.key;
      if (k === "v" || k === "V") this.useTool("select");
      else if (k === "h" || k === "H") this.useTool("hand");
      else if (k === "n" || k === "N") this.tools.querySelector<HTMLButtonElement>(".add-node")?.click();
      else if (k === "0") this.fit();
      else if (k === "+" || k === "=") this.zoomStep(1.25);
      else if (k === "-") this.zoomStep(.8);
      else if (k === "T" && e.shiftKey) this.tidy();
      else if ((k === "Delete" || k === "Backspace") && this.selected) this.remove(this.selected);
      else return;
      e.preventDefault();
    });
  }

  private zoomStep(f: number) {
    const r = this.field.getBoundingClientRect();
    this.zoomAt(this.view.z * f, r.left + r.width / 2, r.top + r.height / 2);
  }

  private useTool(t: Tool) {
    this.tool = t;
    for (const b of this.tools.querySelectorAll<HTMLButtonElement>("[data-tool]")) b.setAttribute("aria-pressed", String(b.dataset.tool === t));
    this.field.classList.toggle("hand", t === "hand");
  }

  // ── toolbar, Layers, Inspector ──────────────────────────────────────

  private toolbar() {
    const button = (label: string, glyph: Parameters<typeof icon>[0] | string, on: () => void, props: Record<string, unknown> = {}) =>
      el("button", { type: "button", class: "icon", "aria-label": label, onclick: on, ...props },
        glyph.length <= 2 ? el("span", { class: "glyph mono", text: glyph }) : icon(glyph as Parameters<typeof icon>[0]));
    const menu = el("div", { class: "add-menu", role: "menu", hidden: true });
    const add = button("Add a node (N)", "plus", () => {
      menu.hidden = !menu.hidden;
      if (menu.hidden) return;
      const names = this.nodes.flatMap((n) => n.defines);
      menu.replaceChildren(...KINDS.map(([kind, what, live]) => {
        if (kind === "Graph") return el("button", { type: "button", role: "menuitem", disabled: !names.length || undefined,
          onclick: () => { menu.hidden = true; this.append([graph(names.at(-1)!)]); } },
          el("span", { text: `Graph of ${names.at(-1) ?? "…"}` }), el("span", { class: "what", text: names.length ? what : "needs a name defined" }));
        return el("button", { type: "button", role: "menuitem", disabled: !live || undefined,
          onclick: () => { menu.hidden = true; if (kind === "Formula") this.append([math({ text: "" })]); } },
          el("span", { text: kind }), el("span", { class: "what", text: live ? what : `${what} · later` }));
      }));
    }, { class: "icon add-node", "aria-haspopup": "menu" });
    document.addEventListener("pointerdown", (e) => { if (!menu.hidden && !(e.target as Element).closest(".add-menu, .add-node")) menu.hidden = true; });
    this.zoom.addEventListener("click", () => this.fit());
    this.tools.append(
      button("Select (V)", "V", () => this.useTool("select"), { "data-tool": "select", "aria-pressed": "true" }),
      button("Hand (H)", "H", () => this.useTool("hand"), { "data-tool": "hand", "aria-pressed": "false" }),
      el("span", { class: "add-wrap" }, add, menu),
      button("Tidy (Shift+T)", "≡", () => this.tidy()),
      button("Fit (0)", "▣", () => this.fit()),
      this.zoom);
  }

  private tree() {
    const rows = this.nodes.map((n) => el("button", { type: "button", class: "row", role: "treeitem", "data-id": n.id, "aria-selected": String(n.id === this.selected),
      onclick: () => { this.select(n.id); this.reveal(n.id); } },
      el("span", { class: "glyph mono", text: GLYPH[n.kind] }), el("span", { class: "name", text: n.id }),
      n.defines.length ? el("span", { class: "defs mono", text: n.defines.join(", ") }) : null));
    this.layers.replaceChildren(el("h2", { text: "Layers" }),
      rows.length ? el("div", { role: "tree", "aria-label": "Layout" }, ...rows) : el("p", { class: "muted", text: "No nodes. Add one with N, or an example from the reference." }));
  }

  private inspect() {
    const n = this.nodes.find((x) => x.id === this.selected);
    if (!n) { this.inspector.replaceChildren(el("h2", { text: "Inspector" }), el("p", { class: "muted", text: "Select a node to see its properties and its answer." })); return; }
    const p = this.placed.find((x) => x.id === n.id);
    const a = this.answers.get(n.cell) ?? null;
    const link = (id: string) => el("button", { type: "button", class: "link mono", text: id, onclick: () => { this.select(id); this.reveal(id); } });
    const props = el("dl", { class: "props" });
    const row = (label: string, ...value: (HTMLElement | string | null)[]) => props.append(el("dt", { text: label }), el("dd", {}, ...value));
    if (n.kind === "formula") {
      const t = el("textarea", { class: "text mono", rows: "3", spellcheck: "false", "aria-label": "Source" });
      t.value = n.text;
      t.addEventListener("change", () => { if (n.cell.kind === "math") { n.cell.written = { text: t.value }; this.build(n.id); this.host.changed(); } });
      row("Source", t);
      row("Defines", n.defines.length ? el("span", { class: "mono", text: n.defines.join(", ") }) : el("span", { class: "gap", text: "–" }));
    } else row("Of", el("span", { class: "mono", text: n.cell.kind === "graph" ? n.cell.of : "" }));
    row("Reads", ...(n.reads.length ? n.reads.flatMap((r, k) => [k ? " " : null, link(r)]) : [el("span", { class: "gap", text: "–" })]));
    const readers = this.nodes.filter((x) => x.reads.includes(n.id));
    row("Read by", ...(readers.length ? readers.flatMap((r, k) => [k ? " " : null, link(r.id)]) : [el("span", { class: "gap", text: "–" })]));
    row("X Y", el("span", { class: "mono", text: p ? `${p.x} ${p.y}` : "" }), n.at ? el("button", { type: "button", class: "quiet small", text: "Unpin",
      onclick: () => { delete (n.cell.settings as { at?: unknown }).at; delete n.at; this.place(); this.inspect(); this.host.changed(); } }) : el("span", { class: "muted", text: " placed" }));
    row("W", el("span", { class: "mono", text: String(n.w) }));
    const tabs = el("div", { class: "route-tabs", role: "tablist" });
    const panel = el("div", { role: "tabpanel" });
    const pages: [string, () => HTMLElement | DocumentFragment][] = [
      ["Properties", () => props],
      ["Answer", () => a ? answerView(a, n.text, { steps: true, source: this.host.source(), ask: this.host.ask }) : el("p", { class: "muted", text: "Not answered yet." })],
      ["Trace", () => el("p", { class: "muted", text: "The trace of this node's asks and answers: later." })]];
    const pick = (k: number) => { tabs.querySelectorAll("button").forEach((b, j) => b.setAttribute("aria-selected", String(j === k))); panel.replaceChildren(pages[k][1]()); this.inspector.dataset.tab = String(k); };
    pages.forEach(([t], k) => tabs.append(el("button", { type: "button", role: "tab", text: t, onclick: () => pick(k) })));
    this.inspector.replaceChildren(el("h2", {}, el("span", { class: "id", text: n.id }), " ", el("span", { class: "kind", text: n.kind === "formula" ? "Formula" : "Graph" })), tabs, panel);
    pick(Number(this.inspector.dataset.tab ?? 0) < pages.length ? Number(this.inspector.dataset.tab ?? 0) : 0);
  }
}

// the toolbar goes into the bar after the page tabs, as the design places it
export function mountTools(tools: HTMLElement) { $(".preset")?.after(tools); }
