// The reference of the question language, as the page reads it: `site/data/reference.json`.
//
// Written from the reader's own tables (`functions`, `operators`, `constants` in
// io/field_parse.hpp), so the reference cannot name a function the reader does not know or miss
// one it does; every example is answered through `respond`, the page's door, so each one shown
// is one the library reads, and its value and LaTeX are the library's. `--check` is a ctest.
//
//   site_reference --write site/data/reference.json
//   site_reference --check site/data/reference.json

#include <rsc/spatium/respond.hpp>
#include <spatium/io/field_parse.hpp>
#include <spatium/io/field_text.hpp>
#include <spatium/io/json.hpp>

#include <cstdio>
#include <format>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>

namespace {

namespace bd = spatium::io::build;
using spatium::io::JsonValue;
using Object = JsonValue::ObjectType;
using Array = JsonValue::ArrayType;

// what the library makes of an example: its LaTeX and its value to twelve digits (the last ones
// are the libm's, and the page shows as many as its setting says), or its route
JsonValue answered(std::string_view example) {
    const JsonValue r = rsc::respond(example);
    Object out;
    for (const auto& [k, v] : r.as_object()) {
        if (k == "route" || k == "refused") out.emplace_back(k, v);
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
    // the four kinds of question, each with an example the library answers
    struct Kind { const char* form; const char* caption; const char* example; };
    for (const Kind& k : {Kind{"an expression", "evaluated, as a calculator does", "2(3 + 4)^2 - 1/4"},
                          Kind{"f = g", "an equation in one unknown of degree at most two, solved in closed form with "
                                        "its steps and checked by Vieta's formulas", "2x^2 - 3x - 5 = 0"},
                          Kind{"solve f = g for x in [a, b], p in [c, d]",
                               "the root on [a, b] as a function of the parameter p, synthesised and certified",
                               "solve x^2 - p = 0 for x in [0, 10], p in [0.5, 4]"},
                          Kind{"integrate f for t in [a, b], p in [c, d]",
                               "the integral as a function of p; a limit may be p itself",
                               "integrate sin(t)/t for t in [0, x], x in [0.5, 6]"}})
        questions.push_back(Object{{"form", k.form}, {"caption", k.caption}, {"example", k.example},
                                   {"answer", answered(k.example)}});
    std::string out = "{\n";
    const auto section = [&](const char* name, const Array& rows, bool last) {
        out += std::string("  \"") + name + "\": [\n";
        for (std::size_t i = 0; i < rows.size(); ++i) out += "    " + rows[i].dump() + (i + 1 < rows.size() ? ",\n" : "\n");
        out += last ? "  ]\n" : "  ],\n";
    };
    section("functions", functions, false);
    section("operators", operators, false);
    section("constants", constants, false);
    section("questions", questions, true);
    return out + "}\n";
}

} // namespace

int main(int argc, char** argv) {
    const std::string_view mode = argc > 2 ? argv[1] : "";
    if (mode != "--write" && mode != "--check") {
        std::fprintf(stderr, "usage: site_reference --write|--check <path to reference.json>\n");
        return 2;
    }
    const std::string now = reference();
    if (mode == "--write") {
        std::ofstream(argv[2]) << now;
        return 0;
    }
    std::ifstream in(argv[2]);
    std::stringstream was;
    was << in.rdbuf();
    if (was.str() == now) return 0;
    std::fprintf(stderr, "reference.json is out of date; regenerate: site_reference --write %s\n", argv[2]);
    return 1;
}
