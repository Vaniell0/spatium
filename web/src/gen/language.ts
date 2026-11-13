// Written by examples/site_reference.cpp from the library's own tables; do not edit.
// Regenerate: site_reference --write site/data/reference.json web/src/gen/language.ts docs/language.md

export const functions = ["sin", "cos", "log", "exp", "sqrt", "abs", "min", "max", "less"] as const;
export type FunctionName = (typeof functions)[number];

export const operators = ["a + b, a - b", "a * b, a / b", "2x, 2(a + b), (a)(b)", "-a", "a^b", "(a)", "a = b"] as const;
export type Operator = (typeof operators)[number];

export const constants = ["pi", "inf", "nan"] as const;
export type Constant = (typeof constants)[number];

export const statements = ["value", "parameter", "root", "integral", "shape"] as const;
export type StatementKind = (typeof statements)[number];

export const relations = ["=", "<", "<=", ">", ">="] as const;
export type Relation = (typeof relations)[number];

export const routes = ["at infinity", "evaluate", "expression", "closed form", "linear system", "interval method", "every root, proved", "root on a bracket", "integral", "synthesis"] as const;
export type Route = (typeof routes)[number];

export const arithmetics = ["double", "float", "rational", "decimal", "real50"] as const;
export type Arithmetic = (typeof arithmetics)[number];

export const shapes = ["torus", "cylinder", "sphere", "flake", "offset", "scatter", "union"] as const;
export type ShapeOperation = (typeof shapes)[number];
