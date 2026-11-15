// The core in a worker, so a long synthesis never holds the page: the page sends a question, the
// worker answers it with the library's JSON. A question made stale by a newer one is not waited
// for -- the page ends this worker and starts another (web/src/core.ts), which is how a synthesis is
// cancelled. Nothing here reads or changes an answer.

let respond = null, rows = null, plot = null, zeta = null;

self.onmessage = async (e) => {
  const m = e.data;
  if (m.type === "start") {
    try {
      importScripts(m.glue);
      const core = await self.SpatiumCore({ wasmBinary: new Uint8Array(m.wasm) });
      respond = core.cwrap("spatium_respond", "string", ["string", "number"]);
      rows = core.cwrap("spatium_rows", "string", ["string"]);
      plot = core.cwrap("spatium_plot", "string", ["string"]);
      zeta = core.cwrap("spatium_zeta", "string", ["string"]);
      self.postMessage({ type: "ready" });
    } catch (err) {
      self.postMessage({ type: "failed", why: String(err && err.message || err) });
    }
    return;
  }
  if (m.type === "ask") {
    let answer;
    try {
      answer = JSON.parse(m.kind === "rows" ? rows(m.payload) : m.kind === "plot" ? plot(m.payload) : m.kind === "zeta" ? zeta(m.payload) : respond(m.payload, 0));
    } catch (err) {
      answer = { refused: `the core failed: ${err && err.message || err}` };
    }
    self.postMessage({ type: "answer", id: m.id, answer });
  }
};
