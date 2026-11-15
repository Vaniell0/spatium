// The on-screen keyboard: a panel of the page's own over the line, with the layers of keys.ts
// as tabs. A key never takes the focus -- its press is prevented, so the line stays focused and
// the key types into it -- and the panel takes no room from the page. Formulas on the keys are
// rendered once, when the panel is built.

import { type Key, layers } from "../keys.js";
import { el, icon, math } from "./dom.js";

type MathField = HTMLElement & { executeCommand(command: unknown): void };

export interface KeysHost {
  field(): MathField | null;   // the row being written
  commit(): void;              // the line is done
  close(): void;
}

export class Keys {
  readonly root: HTMLElement;
  private readonly body = el("div", { class: "keys-body" });
  private readonly tabs: HTMLButtonElement[];

  constructor(private readonly host: KeysHost) {
    this.tabs = layers.map((l, k) => el("button", { type: "button", class: "layer", role: "tab", text: l.name, onclick: () => this.show(k) }));
    const close = el("button", { type: "button", class: "icon small", id: "close-keys", "aria-label": "Close keyboard", onclick: () => host.close() }, icon("x"));
    this.root = el("section", { class: "keys", id: "keys", "aria-label": "Keyboard", hidden: true },
      el("div", { class: "keys-head", role: "tablist" }, ...this.tabs, close), this.body);
    // a press anywhere on the panel leaves the focus where it is
    this.root.addEventListener("pointerdown", (e) => e.preventDefault());
    this.show(0);
  }

  private show(which: number) {
    this.tabs.forEach((t, k) => t.setAttribute("aria-selected", String(k === which)));
    this.body.replaceChildren(...layers[which].keys.map((key) => this.key(key)));
  }

  private key(key: Key): HTMLElement {
    const b = el("button", { type: "button", class: key.command === "commit" ? "key enter" : "key", "data-example": key.example,
      "aria-label": key.tex ? key.insert : undefined, onclick: () => this.press(key) }, key.tex ? math(key.label) : key.label);
    return b;
  }

  private press(key: Key) {
    if (key.command === "commit") { this.host.commit(); return; }
    const mf = this.host.field();
    if (!mf) return;
    if (key.command) mf.executeCommand({ left: "moveToPreviousChar", right: "moveToNextChar", backspace: "deleteBackward" }[key.command]);
    else mf.executeCommand(["insert", key.insert, { format: "latex", selectionMode: "placeholder" }]);
  }
}
