// The core the page asks: the library built to WebAssembly, run in a worker (site/core-worker.js)
// so a long synthesis never holds the page. One question at a time goes to the worker, and each
// asker keeps only its newest; a question made stale while the worker is on it -- a synthesis
// typed past -- is given a moment, then the worker is ended and a new one started. That is how a
// synthesis is cancelled. Nothing here reads or changes an answer.

export type Answer = Record<string, unknown>;
export type Kind = "rows" | "text" | "plot";
export type State = "loading" | "ready" | "busy" | "failed";

export interface CoreOptions {
  glue: string;            // the Emscripten glue, wasm/spatium_core.js
  wasm: string;            // the module, wasm/spatium_core.wasm
  worker: string;          // core-worker.js
  cache: boolean;          // keep the core in the browser's cache
  onState?: (state: State, detail: string) => void;
}

const CACHE = "spatium-core";
const STALE_MS = 400;

interface Waiting { kind: Kind; payload: string; resolve: (a: Answer | null) => void }
interface Inflight { key: object; id: number; resolve: (a: Answer | null) => void; timer: ReturnType<typeof setTimeout> | 0 }

export class Core {
  ready = false;
  failed: string | null = null;
  version = "";
  size = 0;
  from: "cache" | "network" | "" = "";

  private glueUrl = "";
  private wasm: ArrayBuffer | null = null;
  private worker: Worker | null = null;
  private waiting = new Map<object, Waiting>();
  private inflight: Inflight | null = null;
  private asked = 0;

  constructor(private readonly options: CoreOptions) {}

  private state(s: State, detail: string) {
    this.options.onState?.(s, detail);
  }

  async load(): Promise<void> {
    this.state("loading", "Core loading");
    try {
      const [glue, wasm] = await Promise.all([this.fetchBytes(this.options.glue), this.fetchBytes(this.options.wasm)]);
      await this.use(glue.bytes, wasm.bytes);
      this.from = wasm.from;
      if (this.failed) throw new Error(this.failed);
      if (wasm.from === "network") {
        await this.store(this.options.glue, glue.bytes, "text/javascript");
        await this.store(this.options.wasm, wasm.bytes, "application/wasm");
      }
    } catch (e) {
      this.failed = e instanceof Error ? e.message : String(e);
      this.state("failed", `Core did not load: ${this.failed}`);
    }
  }

  // the core on the server, when it is not the one running: true when it was replaced
  async update(): Promise<boolean> {
    const [glue, wasm] = await Promise.all([this.fetchBytes(this.options.glue, true), this.fetchBytes(this.options.wasm, true)]);
    if (this.ready && this.version && (await digest(wasm.bytes)) === this.version) return false;
    await this.store(this.options.glue, glue.bytes, "text/javascript");
    await this.store(this.options.wasm, wasm.bytes, "application/wasm");
    await this.use(glue.bytes, wasm.bytes);
    this.from = "network";
    return true;
  }

  // The newest question of `key` (a cell, a graph): resolves to its answer, or to null when a
  // newer question of the same key took its place.
  ask(key: object, kind: Kind, payload: string): Promise<Answer | null> {
    return new Promise((resolve) => {
      this.waiting.get(key)?.resolve(null);
      this.waiting.set(key, { kind, payload, resolve });
      const f = this.inflight;
      if (f && f.key === key && !f.timer)
        f.timer = setTimeout(() => {
          if (this.inflight !== f) return;
          f.resolve(null);
          this.inflight = null;
          void this.spawn();
        }, STALE_MS);
      this.pump();
    });
  }

  // a key no longer asked for: its question is dropped
  forget(key: object) {
    this.waiting.get(key)?.resolve(null);
    this.waiting.delete(key);
  }

  private pump() {
    if (!this.ready || this.inflight || this.waiting.size === 0 || !this.worker) return;
    const [key, q] = this.waiting.entries().next().value as [object, Waiting];
    this.waiting.delete(key);
    this.inflight = { key, id: ++this.asked, resolve: q.resolve, timer: 0 };
    this.state("busy", "Core working");
    this.worker.postMessage({ type: "ask", id: this.inflight.id, kind: q.kind, payload: q.payload });
  }

  private answered(id: number, answer: Answer) {
    const f = this.inflight;
    if (!f || f.id !== id) return;
    if (f.timer) clearTimeout(f.timer);
    this.inflight = null;
    this.state("ready", "Core ready");
    f.resolve(answer);
    this.pump();
  }

  // the glue and the module are kept and replaced together, so a cached module never meets a newer glue
  private async use(glue: ArrayBuffer, wasm: ArrayBuffer) {
    if (this.glueUrl) URL.revokeObjectURL(this.glueUrl);
    this.glueUrl = URL.createObjectURL(new Blob([glue], { type: "text/javascript" }));
    this.wasm = wasm;
    this.version = await digest(wasm);
    this.size = wasm.byteLength;
    await this.spawn();
  }

  // a worker with the core in it; a new one replaces one ended mid-question
  private spawn(): Promise<void> {
    this.worker?.terminate();
    this.ready = false;
    const w = new Worker(this.options.worker);
    this.worker = w;
    return new Promise((done) => {
      w.onmessage = (e: MessageEvent) => {
        const m = e.data as { type: string; why?: string; id?: number; answer?: Answer };
        if (m.type === "ready") {
          this.ready = true;
          this.failed = null;
          this.state("ready", "Core ready");
          done();
          this.pump();
        } else if (m.type === "failed") {
          this.failed = m.why ?? "the core did not start";
          this.state("failed", `Core failed: ${this.failed}`);
          done();
        } else if (m.type === "answer") this.answered(m.id ?? -1, m.answer ?? {});
      };
      w.onerror = (e) => {
        this.failed = e.message || "the worker did not start";
        this.state("failed", this.failed);
        done();
      };
      w.postMessage({ type: "start", glue: this.glueUrl, wasm: this.wasm!.slice(0) });
    });
  }

  private async fetchBytes(path: string, fresh = false): Promise<{ bytes: ArrayBuffer; from: "cache" | "network" }> {
    if (this.options.cache && !fresh && "caches" in globalThis) {
      try {
        const hit = await (await caches.open(CACHE)).match(path);
        if (hit) return { bytes: await hit.arrayBuffer(), from: "cache" };
      } catch {}
    }
    const res = await fetch(path, { cache: fresh ? "no-store" : "default" });
    if (!res.ok) throw new Error(`${path} answered ${res.status}`);
    return { bytes: await res.arrayBuffer(), from: "network" };
  }

  private async store(path: string, bytes: ArrayBuffer, type: string) {
    if (!this.options.cache || !("caches" in globalThis)) return;
    try {
      await (await caches.open(CACHE)).put(path, new Response(bytes.slice(0), { headers: { "content-type": type } }));
    } catch {}
  }
}

async function digest(bytes: ArrayBuffer): Promise<string> {
  if (!crypto.subtle) return "";
  const h = new Uint8Array(await crypto.subtle.digest("SHA-256", bytes));
  return [...h.slice(0, 4)].map((b) => b.toString(16).padStart(2, "0")).join("");
}
