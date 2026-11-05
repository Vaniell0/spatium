// The on-screen keyboard's one layout, given to MathLive: the functions the language has as keys
// of their own (no shift layer that swaps them in and out), digits, operators and the arrows.
// Each function key inserts its call with a placeholder for the argument, as the reader takes it.

export interface Key { label: string; insert: string; example: string }   // example: the key's call filled in

export const functionKeys: Key[] = [
  { label: "sin", insert: "\\sin\\left(#0\\right)", example: "\\sin\\left(1\\right)" },
  { label: "cos", insert: "\\cos\\left(#0\\right)", example: "\\cos\\left(1\\right)" },
  { label: "ln", insert: "\\ln\\left(#0\\right)", example: "\\ln\\left(2\\right)" },
  { label: "lg", insert: "\\lg\\left(#0\\right)", example: "\\lg\\left(100\\right)" },
  { label: "log_□", insert: "\\log_{#0}\\left(#0\\right)", example: "\\log_{2}\\left(8\\right)" },
  { label: "√", insert: "\\sqrt{#0}", example: "\\sqrt{2}" },
  { label: "□²", insert: "#@^{2}", example: "3^{2}" },
  { label: "□^□", insert: "#@^{#0}", example: "2^{10}" },
];

// MathLive's layout object: rows of keycaps, a keycap a string of LaTeX or {label, insert}
export function layout() {
  const f = functionKeys.map((k) => ({ label: k.label, insert: k.insert, class: "small" }));
  return {
    label: "123",
    rows: [
      [...f.slice(0, 5), "\\pi", "e", "[backspace]"],
      ["7", "8", "9", "\\div", "(", ")", f[5], "x"],
      ["4", "5", "6", "\\times", f[6], f[7], "\\le", "y"],
      ["1", "2", "3", "-", "=", "<", ">", "\\ge"],
      ["0", ".", ",", "+", "[left]", "[right]", { label: "[return]", width: 2 }],
    ],
  };
}
