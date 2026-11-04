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
