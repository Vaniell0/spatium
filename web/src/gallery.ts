// The gallery, a page of its own: solved questions grouped by how strongly each answer holds, read
// from gallery/gallery.json, which the library writes. No core is loaded here; each tile is a link
// that opens its question in the calculator.

import { math, toLink } from "./doc.js";
import { $, el, math as tex } from "./view/dom.js";

interface Card { id: string; question?: string; question_latex?: string[]; typed?: boolean;
  sheet?: { guarantee: string; label: string; ok: boolean; result: string } }
interface Gallery { groups: { guarantee: string; title: string }[]; cards: Card[] }

async function tile(c: Card): Promise<HTMLElement> {
  const rows = c.question_latex?.length ? c.question_latex : [c.question!];
  const href = `../#w=${await toLink({ preset: "calculator", cells: [math(c.question_latex?.length ? { rows } : { text: c.question! })] })}`;
  const s = c.sheet!;
  return el("a", { class: "tile", href },
    el("span", { class: "q" }, ...rows.map((r) => (c.question_latex?.length ? tex(r) : el("span", { class: "mono", text: r })))),
    s.result ? tex(s.result, { cls: "math res" }) : null,
    el("span", { class: s.ok ? "label ok" : "label", text: s.label }));
}

async function main() {
  const box = $("#gallery")!;
  let g: Gallery;
  try {
    g = await (await fetch("gallery.json")).json();
  } catch {
    box.append(el("p", { class: "state", text: "The gallery did not load. Reload the page to try again." }));
    return;
  }
  // by how strongly each answer holds, in the gallery's own order of groups
  for (const group of g.groups) {
    const picks = g.cards.filter((c) => c.question && c.typed !== false && c.sheet?.guarantee === group.guarantee);
    if (!picks.length) continue;
    box.append(el("section", { class: "group" }, el("h2", { text: group.title }),
      el("div", { class: "tiles" }, ...(await Promise.all(picks.map(tile))))));
  }
}

void main();
