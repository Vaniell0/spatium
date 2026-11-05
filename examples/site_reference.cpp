// The reference of the question language, as the page reads it: `site/data/reference.json`, and
// the same words as TypeScript types, `web/src/gen/language.ts`, so the page's code cannot name
// a route, a relation or a shape the library does not have, and `web/src/registry.ts` must give
// each one the library has a place on the page.
//
// Written from the library's own tables (`functions`, `operators`, `constants` in
// io/field_parse.hpp; `statement_kinds`, `relations`, `shape_operations` in io/document.hpp;
// `routes`, `arithmetics` in rsc/spatium/evaluate.hpp), so the reference cannot name a word the
// library does not know or miss one it does; every example is answered through `respond`, the
// page's door, so each one shown is one the library reads, and its value and LaTeX are the
// library's. `--check` is a ctest.
//
//   site_reference --write site/data/reference.json web/src/gen/language.ts
//   site_reference --check site/data/reference.json web/src/gen/language.ts

#include <rsc/spatium/evaluate.hpp>
#include <rsc/spatium/respond.hpp>
#include <spatium/io/document.hpp>
#include <spatium/io/field_parse.hpp>
#include <spatium/io/field_text.hpp>
#include <spatium/io/json.hpp>

#include <cstdio>
#include <format>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
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
    const JsonValue r = rsc::respond(example);
    Object out;
    for (const auto& [k, v] : r.as_object()) {
        if (k == "route" || k == "refused" || k == "question_latex") out.emplace_back(k, v);
        if (k == "value") out.emplace_back("value", std::format("{:.12g}", v.as_number()));
    }
    if (const auto f = bd::parse(example); f && std::string_view(example).find('=') == std::string_view::npos)
        out.emplace_back("latex", *bd::to_latex(*f));
    return out;
}

std::string reference() {
    Array functions, operators, constants, questions;
    for (const auto& f : bd::functions)
        functions.push_back(Object{{"name", std::string(f.name)}, {"arity", f.arity}, {"caption", std::string(f.caption)},
                                   {"domain", std::string(f.domain)}, {"example", std::string(f.example)},
                                   {"answer", answered(f.example)}});
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
                          Kind{"p \\in [c, d]", "a parameter: a row below that reads it is answered as a function of p, "
                                              "synthesised and certified", "p in [0.5, 4]\nx: x^2 - p = 0 on [0, 10]"},
                          Kind{"I = \\int_a^b f\\,dt", "an integral; a limit may be a parameter above",
                               "x in [0.5, 6]\nSi = integral sin(t)/t dt on [0, x]"}})
        questions.push_back(Object{{"form", k.form}, {"caption", k.caption}, {"example", k.example},
                                   {"answer", answered(k.example)}});
    Array statements, relations, routes, arithmetics, shapes;
    for (const auto k : bd::statement_kinds) statements.push_back(Object{{"name", bd::statement_kind_name(k)}});
    for (const auto r : bd::relations)
        relations.push_back(Object{{"written", std::string(r)}, {"latex", bd::detail::relation_latex(r)}});
    for (const auto r : rsc::routes) routes.push_back(Object{{"name", rsc::route_name(r)}, {"caption", rsc::route_caption(r)}});
    for (const auto a : rsc::arithmetics) arithmetics.push_back(Object{{"name", rsc::arithmetic_name(a)}});
    for (const auto& o : bd::shape_operations)
        shapes.push_back(Object{{"name", std::string(o.name)}, {"arguments", std::string(o.arguments)},
                                {"required", o.required}, {"variadic", o.variadic}, {"caption", std::string(o.caption)}});
    const std::pair<const char*, const Array*> sections[] = {
        {"functions", &functions},   {"operators", &operators}, {"constants", &constants},
        {"questions", &questions},   {"statements", &statements}, {"relations", &relations},
        {"routes", &routes},         {"arithmetics", &arithmetics}, {"shapes", &shapes}};
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
                      "// Regenerate: site_reference --write site/data/reference.json web/src/gen/language.ts\n";
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
    list("routes", "Route", names(rsc::routes, rsc::route_name));
    list("arithmetics", "Arithmetic", names(rsc::arithmetics, rsc::arithmetic_name));
    list("shapes", "ShapeOperation", names(bd::shape_operations, [](const auto& o) { return o.name; }));
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
    const std::string_view mode = argc > 3 ? argv[1] : "";
    if (mode != "--write" && mode != "--check") {
        std::fprintf(stderr, "usage: site_reference --write|--check <reference.json> <language.ts>\n");
        return 2;
    }
    const std::string files[] = {reference(), language_ts()};
    int stale = 0;
    for (int k = 0; k < 2; ++k) {
        if (mode == "--write") std::ofstream(argv[2 + k]) << files[k];
        else if (contents(argv[2 + k]) != files[k]) {
            std::fprintf(stderr, "%s is out of date; regenerate: site_reference --write %s %s\n", argv[2 + k], argv[2], argv[3]);
            stale = 1;
        }
    }
    return stale;
}
