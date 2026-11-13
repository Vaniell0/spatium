// The reference of the question language, as the page reads it: `site/data/reference.json`, and
// the same words as TypeScript types, `web/src/gen/language.ts`, so the page's code cannot name
// a route, a relation or a shape the library does not have, and `web/src/registry.ts` must give
// each one the library has a place on the page.
//
// Written from the library's own tables (`functions`, `operators`, `constants` in
// io/field_parse.hpp; `statement_kinds`, `relations`, `shape_operations` in io/document.hpp;
// `route_table`, `arithmetics` in spatium/answer/evaluate.hpp), so the reference cannot name a word the
// library does not know or miss one it does; every example is answered through `respond`, the
// page's door, so each one shown is one the library reads, and its value and LaTeX are the
// library's. `--check` is a ctest.
//
//   site_reference --write site/data/reference.json web/src/gen/language.ts docs/language.md
//   site_reference --check site/data/reference.json web/src/gen/language.ts docs/language.md

#include <spatium/answer/evaluate.hpp>
#include <spatium/answer/respond.hpp>
#include <spatium/io/document.hpp>
#include <spatium/io/field_mathjson.hpp>
#include <spatium/io/field_parse.hpp>
#include <spatium/io/field_text.hpp>
#include <spatium/io/json.hpp>

#include <algorithm>
#include <cstdio>
#include <format>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

namespace bd = spatium::io::build;
using spatium::io::JsonValue;
using Object = JsonValue::ObjectType;
using Array = JsonValue::ArrayType;

// what the library makes of an example: its LaTeX, the rows a page opens it as, and its value to
// twelve digits (the last ones are the libm's, and the page shows as many as its setting says), or
// its route
JsonValue answered(std::string_view example) {
    const JsonValue r = spatium::answer::respond(example);
    Object out;
    for (const auto& [k, v] : r.as_object()) {
        if (k == "route" || k == "refused" || k == "question_latex") out.emplace_back(k, v);
        if (k == "value") out.emplace_back("value", std::format("{:.12g}", v.as_number()));
        if (k == "sheet")
            for (const auto& [sk, sv] : v.as_object())
                if (sk == "guarantee" || sk == "label") out.emplace_back(sk, sv);
    }
    if (const auto f = bd::parse(example); f && std::string_view(example).find('=') == std::string_view::npos)
        out.emplace_back("latex", *bd::to_latex(*f));
    return out;
}

// The editor's notation of a function, where it is not the language's own: `\ln`, `\lg`,
// `\log_b`, `\tan`, `\sqrt[n]`. `rows` is the example as the editor writes it, `example` the line
// of the language it reads as; site/wasm/check.mjs reads the rows back with the editor's parser and
// asks that they answer as the line does. Every MathJSON head the reader takes is here, or is a
// form the Writing and Operators sections show (`forms` below), or main fails.
struct Notation {
    const char* head;      // MathJSON's name, as field_mathjson.hpp has it
    const char* form;      // LaTeX, as a person writes it
    const char* rows;      // the example, in the editor
    const char* example;   // the same, in the language
    const char* caption;   // where the language has no function of this name
    const char* domain;
};
constexpr Notation notation[] = {
    {"Sin", "\\sin x", "\\sin\\left(\\frac{\\pi}{6}\\right)", "sin(pi / 6)", "", ""},
    {"Cos", "\\cos x", "\\cos\\left(\\frac{\\pi}{3}\\right)", "cos(pi / 3)", "", ""},
    {"Tan", "\\tan x", "\\tan\\left(\\frac{\\pi}{4}\\right)", "sin(pi / 4) / cos(pi / 4)", "tangent, sin x / cos x",
     "every x where cos x is not 0"},
    {"Exp", "\\exp x", "\\exp\\left(1\\right)", "exp(1)", "", ""},
    {"Ln", "\\ln x", "\\ln\\left(2\\right)", "log(2)", "", ""},
    {"Lg", "\\lg x", "\\lg\\left(100\\right)", "log(100) / log(10)", "logarithm in base 10; \\log x alone is the same",
     "x > 0"},
    {"Log", "\\log_{b} x", "\\log_{3}\\left(81\\right)", "log(81) / log(3)", "logarithm in base b", "x > 0, b > 0, b not 1"},
    {"Lb", "\\log_{2} x", "\\log_{2}\\left(8\\right)", "log(8) / log(2)", "logarithm in base 2", "x > 0"},
    {"Sqrt", "\\sqrt{x}", "\\sqrt{2}", "sqrt(2)", "", ""},
    {"Root", "\\sqrt[n]{x}", "\\sqrt[3]{8}", "8^(1/3)", "the n-th root", "x > 0"},
    {"Abs", "\\left|x\\right|", "\\left|-3\\right|", "abs(-3)", "", ""},
    {"Min", "\\min\\left(x, y\\right)", "\\min\\left(2, 3\\right)", "min(2, 3)", "", ""},
    {"Max", "\\max\\left(x, y\\right)", "\\max\\left(2, 3\\right)", "max(2, 3)", "", ""},
    {"Boole", "\\left[x < y\\right]", "\\left[1 < 2\\right]", "less(1, 2)", "", ""},
};
// the heads that are forms, not functions: the Writing and Operators sections show them
constexpr std::string_view forms[] = {"Add", "Subtract", "Negate", "Multiply", "InvisibleOperator", "Divide", "Rational",
                                      "Half", "Power", "Square", "Delimiter", "Equal", "Less", "LessEqual", "Greater",
                                      "GreaterEqual", "Error"};

// a function's derivative in each argument, as its row of the table writes it; none where the
// notation is no function of the language
Array derivative(const bd::Function* f) {
    Array out;
    if (f)
        for (int k = 0; k < f->arity; ++k) out.push_back(std::string(f->derivative[static_cast<std::size_t>(k)]));
    return out;
}

// the language's function a notation's example calls first, if it is one
const bd::Function* function_of(std::string_view example) {
    for (const auto& f : bd::functions)
        if (example.starts_with(std::string(f.name) + "(")) return &f;
    return nullptr;
}

std::string reference() {
    // every head the reader takes -- a form, a composition, or a function of the language by its row
    std::vector<std::string_view> heads;
    for (const auto& h : bd::mathjson_heads) heads.push_back(h.name);
    for (const auto& f : bd::functions)
        if (!f.mathjson.empty()) heads.push_back(f.mathjson);
    for (const auto h : heads) {
        bool placed = std::ranges::find(forms, h) != std::end(forms);
        for (const auto& n : notation) placed = placed || h == n.head;
        if (!placed) throw std::runtime_error(std::format("the editor's {} has no row in the reference", h));
    }
    Array writing, functions, operators, constants, questions;
    // how the formula editor takes keystrokes: what is typed, and the line of the language the rows
    // it shows read as; the rows are the library's own (`question_latex`)
    struct Typed { const char* typed; const char* reads; const char* caption; };
    for (const Typed& t : {Typed{"sqrt", "sqrt(2)", "a square root"}, Typed{"/", "1 / 3", "a fraction"},
                           Typed{"^", "2^10", "a power; the arrow keys leave it"}, Typed{"pi", "pi", "the number pi"},
                           Typed{"<=", "x: x^2 <= 4", "an inequality; >= gives the other"},
                           Typed{"in", "x: x^2 = 2 on [0, 2]", "where to look for the roots"},
                           Typed{"Shift+Enter", "x, y: x + y = 3, x - y = 1", "a new row; rows of equations are a system"}})
        writing.push_back(Object{{"typed", t.typed}, {"reads", t.reads}, {"caption", t.caption}, {"answer", answered(t.reads)}});
    for (const Notation& n : notation) {
        const bd::Function* f = *n.caption ? nullptr : function_of(n.example);
        if (!*n.caption && !f) throw std::runtime_error(std::format("{}: no caption, and no function of the language", n.head));
        std::string reads;   // what the reader makes of it, as field_mathjson.hpp says
        for (const auto& h : bd::mathjson_heads)
            if (h.name == n.head) reads = std::string(h.caption);
        for (const auto& g : bd::functions)
            if (g.mathjson == n.head) reads = std::string(g.name) + (g.arity == 1 ? "(a)" : "(a, b)");
        functions.push_back(Object{{"head", n.head}, {"form", n.form}, {"reads", reads}, {"name", f ? std::string(f->name) : std::string()}, {"arity", f ? f->arity : 1},
                                   {"caption", f ? std::string(f->caption) : std::string(n.caption)},
                                   {"domain", f ? std::string(f->domain) : std::string(n.domain)}, {"derivative", derivative(f)},
                                   {"rows", Array{n.rows}}, {"example", n.example}, {"answer", answered(n.example)}});
    }
    // a function of the language no notation calls: the language's text only
    for (const auto& f : bd::functions) {
        bool shown = false;
        for (const auto& n : notation) shown = shown || std::string_view(n.example).starts_with(std::string(f.name) + "(");
        if (!shown)
            functions.push_back(Object{{"name", std::string(f.name)}, {"arity", f.arity}, {"caption", std::string(f.caption)},
                                       {"domain", std::string(f.domain)}, {"derivative", derivative(&f)}, {"example", std::string(f.example)},
                                       {"answer", answered(f.example)}});
    }
    for (const auto& o : bd::operators)
        operators.push_back(Object{{"written", std::string(o.written)}, {"caption", std::string(o.caption)},
                                   {"example", std::string(o.example)}, {"answer", answered(o.example)}});
    for (const auto& c : bd::constants)
        constants.push_back(Object{{"name", std::string(c.name)}, {"caption", std::string(c.caption)}});
    // the kinds of question, each written as the rows of a formula editor (`form`, LaTeX), with an
    // example the library answers; the page opens an example as the rows the answer gives back
    struct Kind { const char* form; const char* caption; const char* example; };
    for (const Kind& k : {Kind{"f", "an expression, evaluated, its digits proved", "2(3 + 4)^2 - 1/4"},
                          Kind{"f = g", "an equation in one unknown: in closed form with its steps where it is a quadratic, "
                                        "and every real root proved", "2x^2 - 3x - 5 = 0"},
                          Kind{"f = g,\\ x \\in [a, b]", "every root on the interval, each proved alone in its own",
                               "x: sin(x) = 0.5 on [0, 10]"},
                          Kind{"f > g", "an inequality: where it holds, by the interval method", "x: x^2 - 4 > 0"},
                          Kind{"f = g \\\\ h = k", "a system: one row per equation, an equation per unknown",
                               "x, y: x + y = 3, x - y = 1"},
                          Kind{"p \\in [c, d]", "a parameter: a row below that reads it is answered as a function of p, "
                                              "synthesised and certified", "p in [0.5, 4]\nx: x^2 - p = 0 on [0, 10]"},
                          Kind{"I = \\int_a^b f\\,dt", "an integral; a limit may be a parameter above",
                               "x in [0.5, 6]\nSi = integral sin(t)/t dt on [0, x]"},
                          Kind{"\\frac{d}{dx} f", "a derivative, by the rule of each operation, checked against a dual number",
                               "d/dx sin(x)^2"},
                          Kind{"\\frac{d}{dx} f,\\ x = a", "a derivative at a point, its digits proved", "d/dx x^3 / 3 at 2"}})
        questions.push_back(Object{{"form", k.form}, {"caption", k.caption}, {"example", k.example},
                                   {"answer", answered(k.example)}});
    // the canvas's text: a line of the language a definition, each example read by the library and
    // written back as it reads it, with its row of LaTeX; every kind of statement has a line here
    Array language, told;
    std::vector<bd::StatementKind> seen;
    struct Line { const char* form; const char* caption; const char* example; };
    for (const Line& l : {Line{"name = f", "a value, named; a line below reads it by its name", "r = 0.5\nA = pi * r^2"},
                          Line{"x: f = g", "the roots of an equation in x", "x: 2x^2 - 3x - 5 = 0"},
                          Line{"x: f = g on [a, b]", "every root on the interval", "x: sin(x) = 0.5 on [0, 10]"},
                          Line{"x: f > g", "where an inequality holds; also <, <= and >=", "x: x^2 - 4 > 0"},
                          Line{"x, y: f = g, h = k", "a system: an equation for each unknown", "x, y: x + y = 3, x - y = 1"},
                          Line{"p in [a, b]", "a parameter: a line that reads it is answered as a function of it",
                               "p in [0.5, 4]\nx: x^2 - p = 0 on [0, 10]"},
                          Line{"I = integral f dt on [a, b]", "an integral in the variable after d", "I = integral t^2 dt on [0, 1]"},
                          Line{"g = d/dx f", "a derivative in the variable after d/d: a function of it, or its value at a point after at",
                               "g = d/dx sin(x)^2\nh = d/dx x^3 at 2"},
                          Line{"S = torus(R, r)", "a shape, from the operations below; another reads it by name",
                               "S0 = torus(2, 0.5)\nS1 = offset(S0, 0.1)"},
                          Line{"# words", "a comment, to the end of the line", "r = 2  # the radius"}}) {
        const auto doc = bd::read_document(l.example);
        if (!doc) throw std::runtime_error(std::format("{}: {}", l.example, doc.error().message));
        Array latex;
        for (const auto& st : doc->statements) {
            seen.push_back(st.kind);
            if (const auto rows = bd::to_latex_rows(st))
                for (const auto& row : *rows) latex.push_back(row);
        }
        language.push_back(Object{{"form", l.form}, {"caption", l.caption}, {"example", l.example}, {"latex", latex}});
    }
    for (const auto k : bd::statement_kinds)
        if (std::ranges::find(seen, k) == seen.end())
            throw std::runtime_error(std::format("no line of the language shows a {}", bd::statement_kind_name(k)));
    // what the library says of a document it does not take, in its own words
    for (const auto& [caption, example] : {std::pair{"a name read before a line defines it", "A = pi * r^2"},
                                           std::pair{"a name defined twice", "r = 1\nr = 2"},
                                           std::pair{"a cycle: a line reads only the lines above it", "a = b + 1\nb = a - 1"}}) {
        const auto doc = bd::read_document(example);
        if (doc) throw std::runtime_error(std::format("{}: read, and should not be", example));
        told.push_back(Object{{"caption", caption}, {"example", example}, {"told", doc.error().message}});
    }
    Array statements, relations, routes, arithmetics, shapes, guarantees;
    for (const auto k : bd::statement_kinds) statements.push_back(Object{{"name", bd::statement_kind_name(k)}});
    for (const auto r : bd::relations)
        relations.push_back(Object{{"written", std::string(r)}, {"latex", bd::detail::relation_latex(r)}});
    for (const auto& r : spatium::answer::route_table) routes.push_back(Object{{"name", r.name}, {"caption", r.caption}});
    for (const auto a : spatium::answer::arithmetics) arithmetics.push_back(Object{{"name", spatium::answer::arithmetic_name(a)}});
    for (const auto& o : bd::shape_operations)
        shapes.push_back(Object{{"name", std::string(o.name)}, {"arguments", std::string(o.arguments)},
                                {"required", o.required}, {"variadic", o.variadic}, {"caption", std::string(o.caption)}});
    for (const auto& g : spatium::answer::guarantees)
        guarantees.push_back(Object{{"name", std::string(g.name)}, {"title", std::string(g.title)}, {"promise", std::string(g.promise)}});
    const std::pair<const char*, const Array*> sections[] = {
        {"writing", &writing},       {"functions", &functions},   {"operators", &operators},
        {"constants", &constants},   {"questions", &questions},   {"language", &language},
        {"told", &told},             {"guarantees", &guarantees}, {"statements", &statements},
        {"relations", &relations},   {"routes", &routes},         {"arithmetics", &arithmetics},
        {"shapes", &shapes}};
    std::string out = "{\n";
    for (std::size_t s = 0; s < std::size(sections); ++s) {
        const auto& [name, rows] = sections[s];
        out += std::string("  \"") + name + "\": [\n";
        for (std::size_t i = 0; i < rows->size(); ++i) out += "    " + (*rows)[i].dump() + (i + 1 < rows->size() ? ",\n" : "\n");
        out += s + 1 < std::size(sections) ? "  ],\n" : "  ]\n";
    }
    return out + "}\n";
}

// The words of the language as TypeScript: a constant list each, and its union type. The page's
// registry is a record over each union, so a word the library gains is a type error there until
// the page gives it a place.
std::string language_ts() {
    std::string out = "// Written by examples/site_reference.cpp from the library's own tables; do not edit.\n"
                      "// Regenerate: site_reference --write site/data/reference.json web/src/gen/language.ts docs/language.md\n";
    const auto list = [&](const char* name, const char* type, const std::vector<std::string>& words) {
        out += std::format("\nexport const {} = [", name);
        for (std::size_t i = 0; i < words.size(); ++i) out += (i ? ", " : "") + JsonValue(words[i]).dump();
        out += std::format("] as const;\nexport type {} = (typeof {})[number];\n", type, name);
    };
    std::vector<std::string> w;
    const auto names = [&](auto&& range, auto&& name) {
        w.clear();
        for (const auto& x : range) w.push_back(std::string(name(x)));
        return w;
    };
    list("functions", "FunctionName", names(bd::functions, [](const auto& f) { return f.name; }));
    list("operators", "Operator", names(bd::operators, [](const auto& o) { return o.written; }));
    list("constants", "Constant", names(bd::constants, [](const auto& c) { return c.name; }));
    list("statements", "StatementKind", names(bd::statement_kinds, bd::statement_kind_name));
    list("relations", "Relation", names(bd::relations, [](std::string_view r) { return r; }));
    list("routes", "Route", names(spatium::answer::route_table, [](const auto& r) { return r.name; }));
    list("arithmetics", "Arithmetic", names(spatium::answer::arithmetics, spatium::answer::arithmetic_name));
    list("shapes", "ShapeOperation", names(bd::shape_operations, [](const auto& o) { return o.name; }));
    return out;
}

// The same reference as a page of the docs, `docs/language.md`, written from the JSON above so the
// two cannot differ: LaTeX in `$...$`, the language's text in code.
std::string language_md(const std::string& json) {
    const JsonValue r = *spatium::io::parse(json);
    const auto rows = [&](const char* name) -> const Array& {
        for (const auto& [k, v] : r.as_object())
            if (k == name) return v.as_array();
        throw std::runtime_error(name);
    };
    // a table cell: a pipe would end it, a line break too
    const auto cell = [](std::string s) {
        std::string out;
        for (const char c : s) out += c == '|' ? std::string("\\|") : c == '\n' ? std::string("<br>") : std::string(1, c);
        return out;
    };
    const auto tex = [&](const std::string& s) { return s.empty() ? std::string() : cell("$" + s + "$"); };
    const auto code = [&](const std::string& s) { return s.empty() ? std::string() : "`" + cell(s) + "`"; };
    const auto shown = [&](const JsonValue& row) {
        std::string out;
        if (const JsonValue* a = row.find("answer"); a && a->find("question_latex"))
            for (const auto& l : a->find("question_latex")->as_array()) out += (out.empty() ? "" : "<br>") + tex(l.as_string());
        return out;
    };
    std::string out = "# The language\n\n"
                      "Written by `examples/site_reference.cpp` from the library's own tables, as the site's reference "
                      "is; do not edit. Every example here is answered by the library, and its formula is the one "
                      "it writes.\n\n## Writing a question\n\nIn the formula editor, what is typed becomes the formula "
                      "shown, which reads as the line of the language beside it.\n\n| type | shows | reads as | |\n|---|---|---|---|\n";
    for (const auto& w : rows("writing"))
        out += std::format("| {} | {} | {} | {} |\n", code(w.string_or("typed", "")), shown(w), code(w.string_or("reads", "")),
                           cell(w.string_or("caption", "")));
    out += "\n## Questions\n\n| form | what it asks | example | answer |\n|---|---|---|---|\n";
    for (const auto& q : rows("questions"))
        out += std::format("| {} | {} | {} | {} |\n", tex(q.string_or("form", "")), cell(q.string_or("caption", "")), shown(q),
                           cell(q.find("answer")->string_or("label", "")));
    out += "\n## Functions\n\nThe editor's notation, the language's text, where each has a value, its derivative in each "
           "argument a and b, and an example.\n\n"
           "| written | in the language | | defined for | derivative | example |\n|---|---|---|---|---|---|\n";
    for (const auto& f : rows("functions")) {
        std::string d;
        if (const JsonValue* ds = f.find("derivative"))
            for (const auto& x : ds->as_array()) d += (d.empty() ? "" : ", ") + code(x.as_string());
        out += std::format("| {} | {} | {} | {} | {} | {} |\n", tex(f.string_or("form", "")),
                           code(f.string_or("name", "").empty() ? f.string_or("example", "") : f.string_or("name", "")),
                           cell(f.string_or("caption", "")), cell(f.string_or("domain", "")), d, code(f.string_or("example", "")));
    }
    out += "\n## Operators\n\nTightest last.\n\n| written | | example |\n|---|---|---|\n";
    for (const auto& o : rows("operators"))
        out += std::format("| {} | {} | {} |\n", code(o.string_or("written", "")), cell(o.string_or("caption", "")),
                           code(o.string_or("example", "")));
    out += "\n## Constants\n\n| name | |\n|---|---|\n";
    for (const auto& c : rows("constants"))
        out += std::format("| {} | {} |\n", code(c.string_or("name", "")), cell(c.string_or("caption", "")));
    out += "\n## The language of the canvas\n\nA document is lines, a definition each; a line reads the names the "
           "lines above it define.\n\n| line | | example | as a formula |\n|---|---|---|---|\n";
    for (const auto& l : rows("language")) {
        std::string latex;
        for (const auto& x : l.find("latex")->as_array()) latex += (latex.empty() ? "" : "<br>") + tex(x.as_string());
        out += std::format("| {} | {} | {} | {} |\n", code(l.string_or("form", "")), cell(l.string_or("caption", "")),
                           code(l.string_or("example", "")), latex);
    }
    out += "\nShapes:\n\n| operation | |\n|---|---|\n";
    for (const auto& sh : rows("shapes"))
        out += std::format("| {} | {} |\n", code(sh.string_or("name", "") + "(" + sh.string_or("arguments", "") + ")"),
                           cell(sh.string_or("caption", "")));
    out += "\n`n` is a number, `s` a shape.\n\nWhat the library says of a document it does not take:\n\n"
           "| | document | told |\n|---|---|---|\n";
    for (const auto& t : rows("told"))
        out += std::format("| {} | {} | {} |\n", cell(t.string_or("caption", "")), code(t.string_or("example", "")),
                           cell(t.string_or("told", "")));
    out += "\n## Answers\n\nHow strongly an answer holds, strongest first.\n\n| | promises |\n|---|---|\n";
    for (const auto& g : rows("guarantees"))
        out += std::format("| {} | {} |\n", cell(g.string_or("title", "")), cell(g.string_or("promise", "")));
    return out;
}

// the file at `path` as it is now, or empty
std::string contents(const char* path) {
    std::ifstream in(path);
    std::stringstream was;
    was << in.rdbuf();
    return was.str();
}

} // namespace

int main(int argc, char** argv) {
    const std::string_view mode = argc > 4 ? argv[1] : "";
    if (mode != "--write" && mode != "--check") {
        std::fprintf(stderr, "usage: site_reference --write|--check <reference.json> <language.ts> <language.md>\n");
        return 2;
    }
    const std::string json = reference();
    const std::string files[] = {json, language_ts(), language_md(json)};
    int stale = 0;
    for (int k = 0; k < 3; ++k) {
        if (mode == "--write") std::ofstream(argv[2 + k]) << files[k];
        else if (contents(argv[2 + k]) != files[k]) {
            std::fprintf(stderr, "%s is out of date; regenerate: site_reference --write %s %s %s\n", argv[2 + k], argv[2], argv[3], argv[4]);
            stale = 1;
        }
    }
    return stale;
}
