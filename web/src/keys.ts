// The on-screen keyboard's layers: three of them, eight keys by four each, with no shift layer
// that swaps keys in and out. A key inserts LaTeX with placeholders for what is still to be
// written, as the reader takes it, or runs one of the line's own commands. A key whose insert is
// a call carries an example of that call filled in, which a test reads through the parser: no key
// inserts what the page would then refuse.

export interface Key {
  label: string;        // what the key shows; LaTeX where `tex` is set, so it is rendered once
  tex?: boolean;
  insert?: string;      // LaTeX put at the caret, `#0` a placeholder, `#@` the selection
  command?: "left" | "right" | "backspace" | "commit";
  example?: string;     // the insert filled in, for the test
}

const k = (label: string, insert = label, example?: string): Key => ({ label, insert, example });
const tex = (label: string, insert: string, example?: string): Key => ({ label, tex: true, insert, example });
const fn = (name: string, example: string): Key => k(name, `\\${name}\\left(#0\\right)`, example);
const LEFT: Key = { label: "←", command: "left" }, RIGHT: Key = { label: "→", command: "right" };
const BACK: Key = { label: "⌫", command: "backspace" }, ENTER: Key = { label: "⏎", command: "commit" };

export const layers: { name: string; keys: Key[] }[] = [
  { name: "123", keys: [
    k("7"), k("8"), k("9"), k("÷", "\\div"), k("("), k(")"), LEFT, RIGHT,
    k("4"), k("5"), k("6"), k("×", "\\times"), k("["), k("]"), tex("\\square^{\\square}", "#@^{#0}", "2^{10}"), tex("\\square^{2}", "#@^{2}", "3^{2}"),
    k("1"), k("2"), k("3"), k("−", "-"), tex("\\frac{\\square}{\\square}", "\\frac{#0}{#0}", "\\frac{1}{2}"), tex("\\sqrt{\\square}", "\\sqrt{#0}", "\\sqrt{2}"), k("π", "\\pi"), k("y"),
    k("0"), k("."), k(","), k("+"), k("="), k("x"), BACK, ENTER,
  ] },
  { name: "f(x)", keys: [
    fn("sin", "\\sin\\left(1\\right)"), fn("cos", "\\cos\\left(1\\right)"), fn("ln", "\\ln\\left(2\\right)"), fn("lg", "\\lg\\left(100\\right)"),
    tex("\\log_{\\square}", "\\log_{#0}\\left(#0\\right)", "\\log_{2}\\left(8\\right)"), tex("e^{\\square}", "e^{#0}", "e^{1}"),
    tex("\\left|\\square\\right|", "\\left|#0\\right|", "\\left|-2\\right|"), tex("\\sqrt{\\square}", "\\sqrt{#0}", "\\sqrt{2}"),
    k("min", "\\min\\left(#0, #0\\right)", "\\min\\left(1, 2\\right)"), k("max", "\\max\\left(#0, #0\\right)", "\\max\\left(1, 2\\right)"),
    tex("\\int_{\\square}^{\\square}", "\\int_{#0}^{#0} #0\\,dx", "\\int_{0}^{1} x^{2}\\,dx"), k("∈", "\\in"),
    k("<"), k(">"), k("≤", "\\le"), k("≥", "\\ge"),
    k("π", "\\pi"), k("e"), k("∞", "\\infty"), k("x"), k("t"), k("p"), k("("), k(")"),
    k("["), k("]"), k(","), k("="), LEFT, RIGHT, BACK, ENTER,
  ] },
  { name: "abc", keys: [
    ..."abcdefghijklmnopqrstuvwxyz".split("").map((c) => k(c)),
    k(","), k("="), LEFT, RIGHT, BACK, ENTER,
  ] },
];

// every key that inserts a call, with the example the test reads
export const functionKeys = layers.flatMap((l) => l.keys).filter((key): key is Key & { example: string } => Boolean(key.example));
