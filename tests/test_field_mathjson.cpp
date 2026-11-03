// MathJSON, what a LaTeX editor on a page writes, read as the text of the language.
//
// Each claim is held against a second route: a document whose field is MathJSON must be the same
// document as the one written in the language by hand (`same_document`, field by field), and the
// text it is read as must be that hand-written text, with parentheses only where they are needed.
// The cases it must fail: a function the language does not have, named with where it stood; the
// editor's own error; a complex number; the wrong number of arguments; a name the language cannot
// read; a number written in a form it does not write.

#include <catch2/catch_test_macros.hpp>
#include <spatium/io/document.hpp>
#include <spatium/io/field_mathjson.hpp>
#include <spatium/io/json.hpp>

#include <algorithm>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace bd = spatium::io::build;
using spatium::io::JsonValue;

namespace {

JsonValue json(const std::string& text) {
    auto v = spatium::io::parse(text);
    REQUIRE(v);
    return *v;
}

// a document of the parameter p and a value a, the value written as MathJSON
bd::Document from_mathjson(const std::string& value) {
    auto doc = bd::document_from_json(json(
        R"({"statements": [{"kind": "parameter", "name": "p", "lo": 1, "hi": 2}, {"kind": "value", "name": "a", "value": )" +
        value + "}]}"));
    INFO(value);
    REQUIRE(doc);
    return *doc;
}

std::string refusal(const std::string& value) {
    auto doc = bd::document_from_json(json(
        R"({"statements": [{"kind": "parameter", "name": "p", "lo": 1, "hi": 2}, {"kind": "value", "name": "a", "value": )" +
        value + "}]}"));
    INFO(value);
    REQUIRE_FALSE(doc);
    return doc.error().message;
}

} // namespace

TEST_CASE("every function of the table reads as the text a person writes", "[io][mathjson]") {
    const std::vector<std::pair<std::string, std::string>> cases{
        {R"(["Add", "p", 1])", "p + 1"},
        {R"(["Add", "p", -1])", "p - 1"},
        {R"(["Add", ["Multiply", 2, ["Square", "p"]], ["Multiply", -3, "p"], -5])", "2 * p^2 - 3 * p - 5"},
        {R"(["Add", "p", ["Multiply", -1, "p", "p"]])", "p - p * p"},
        {R"(["Add", "p", ["Negate", ["Add", "p", 1]]])", "p - (p + 1)"},
        {R"(["Subtract", "p", ["Add", "p", 1]])", "p - (p + 1)"},
        {R"(["Negate", ["Power", "p", 2]])", "-p^2"},
        {R"(["Negate", ["Add", "p", 1]])", "-(p + 1)"},
        {R"(["Multiply", ["Add", "p", 1], "p"])", "(p + 1) * p"},
        {R"(["Divide", 1, ["Multiply", 2, "p"]])", "1 / (2 * p)"},
        {R"(["Rational", 1, 3])", "1 / 3"},
        {R"(["Half", "p"])", "p / 2"},
        {R"(["Power", -2, "p"])", "(-2)^p"},
        {R"(["Power", "p", ["Rational", 1, 2]])", "p^(1 / 2)"},
        {R"(["Power", "p", -1])", "p^(-1)"},
        {R"(["Root", "p", 3])", "p^(1 / 3)"},
        {R"(["Sqrt", ["Add", "p", 1]])", "sqrt(p + 1)"},
        {R"(["Exp", "p"])", "exp(p)"},
        {R"(["Ln", "p"])", "log(p)"},
        {R"(["Log", "p"])", "log(p) / log(10)"},
        {R"(["Log", "p", 2])", "log(p) / log(2)"},
        {R"(["Sin", "p"])", "sin(p)"},
        {R"(["Cos", "p"])", "cos(p)"},
        {R"(["Tan", "p"])", "sin(p) / cos(p)"},
        {R"(["Abs", ["Negate", "p"]])", "abs(-p)"},
        {R"(["Max", "p", 1])", "max(p, 1)"},
        {R"(["Min", "p", 1])", "min(p, 1)"},
        {R"(["Delimiter", ["Add", "p", 1]])", "(p + 1)"},
        {R"({"sym": "Pi"})", "pi"},
        {R"(["Multiply", 2, {"sym": "Pi"}])", "2 * pi"},
        {R"({"sym": "ExponentialE"})", "exp(1)"},
        {R"(["Add", "p", {"sym": "PositiveInfinity"}])", "p + inf"},
        {R"(["Multiply", "p", {"sym": "NaN"}])", "p * nan"},
        {R"({"num": "0.1"})", "0.1"},
        {R"({"num": "-2.5e-3"})", "-2.5e-3"},
        {R"(2.5)", "2.5"},
        {R"(-3)", "-3"},
        {R"(["Multiply", 0.00001, "p"])", "1e-05 * p"},
        // as the editor's parser writes them before any rewriting (latex-syntax 0.151.0, measured
        // in a browser): 2p^2 - 3p - 5, a product with no sign, log_2, log_10 and e
        {R"(["Subtract", ["Subtract", ["InvisibleOperator", 2, ["Power", "p", 2]], ["InvisibleOperator", 3, "p"]], 5])", "2 * p^2 - 3 * p - 5"},
        {R"(["InvisibleOperator", 2, "Pi", "p"])", "2 * pi * p"},
        {R"(["InvisibleOperator", ["Delimiter", ["Add", "p", 1]], ["Delimiter", ["Subtract", "p", 2]]])", "(p + 1) * (p - 2)"},
        {R"(["Add", "p", ["InvisibleOperator", -3, "p"]])", "p - 3 * p"},
        {R"(["Lb", "p"])", "log(p) / log(2)"},
        {R"(["Lg", "p"])", "log(p) / log(10)"},
        {R"(["Power", "e", ["InvisibleOperator", ["Negate", "p"], ["Power", "p", 2]]])", "exp(-p * p^2)"},
        {R"(["Power", {"sym": "ExponentialE"}, 2])", "exp(2)"},
        {R"(["Multiply", 2, "e"])", "2 * exp(1)"},
    };
    for (const auto& [mathjson, text] : cases) {
        INFO(mathjson);
        auto t = bd::mathjson_text(json(mathjson));
        REQUIRE(t);
        CHECK(*t == text);
        auto by_hand = bd::read_document("p in [1, 2]\na = " + text);
        REQUIRE(by_hand);
        CHECK(bd::same_document(from_mathjson(mathjson), *by_hand));
    }

    // the table and the cases agree: every function and constant the table names has a case above
    // (an equation and the editor's error have theirs below)
    const auto covered = [&](std::string_view quoted) {
        return std::ranges::any_of(cases, [&](const auto& c) { return c.first.find(quoted) != std::string::npos; });
    };
    for (const bd::MathJsonHead& h : bd::mathjson_heads)
        if (h.name != "Equal" && h.name != "Error") CHECK(covered("\"" + std::string(h.name) + "\""));
    for (const bd::MathJsonConstant& c : bd::mathjson_constants)
        CHECK(covered("\"" + std::string(c.name) + "\""));
}

TEST_CASE("an equation as MathLive writes it is the root a person writes", "[io][mathjson]") {
    // the canonical MathJSON of the LaTeX 2x^2-3x-5=0
    auto doc = bd::document_from_json(json(R"({"statements": [{"kind": "root", "name": "x", "lhs":
        ["Equal", ["Add", ["Multiply", 2, ["Square", "x"]], ["Multiply", -3, "x"], -5], 0]}]})"));
    REQUIRE(doc);
    auto by_hand = bd::read_document("x: 2x^2 - 3x - 5 = 0");
    REQUIRE(by_hand);
    CHECK(bd::same_document(*doc, *by_hand));
    CHECK(bd::to_text(doc->statements[0]) == "x: 2 * x^2 - 3 * x - 5 = 0");

    // a subscripted name is a name, and the ends of an interval may be numbers
    auto sub = bd::document_from_json(json(R"({"statements": [{"kind": "parameter", "name": "x_1", "lo": 0, "hi": 0.5},
        {"kind": "root", "name": "y", "lhs": ["Equal", ["Square", "y"], {"sym": "x_1"}], "lo": 0, "hi": 1}]})"));
    REQUIRE(sub);
    auto sub_by_hand = bd::read_document("x_1 in [0, 0.5]\ny: y^2 = x_1 on [0, 1]");
    REQUIRE(sub_by_hand);
    CHECK(bd::same_document(*sub, *sub_by_hand));
}

TEST_CASE("a row of the editor is a line of a document", "[io][mathjson]") {
    // each row as the editor's parser writes it (latex-syntax 0.151.0, measured in a browser),
    // beside the LaTeX typed and the line it reads as
    const std::vector<std::pair<std::string, std::string>> rows{
        {R"(["Element", "p", ["Interval", {"num": "0.5"}, 4]])", "p in [0.5, 4]"},                       // p\in[0.5,4]
        {R"(["Delimiter", ["Sequence", ["Equal", ["Subtract", ["Power", "x", 2], "p"], {"num": "0"}],
             ["Element", "x", ["Interval", {"num": "0"}, 10]]], {"str": ","}])", "x: x^2 - p = 0 on [0, 10]"},  // x^2-p=0,\ x\in[0,10]
        {R"(["Equal", "I", ["Integrate", ["Power", "e", ["InvisibleOperator", ["Negate", "p"], ["Power", "t", 2]]],
             ["Tuple", "t", 0, 1]]])", "I = integral exp(-p * t^2) dt on [0, 1]"},                      // I=\int_0^1e^{-pt^2}\,dt
        {R"(["Equal", "a", ["Add", 2, 2]])", "a = 2 + 2"},
        {R"(["Integrate", "x", ["Tuple", "x", 0, 1]])", "answer = integral x dx on [0, 1]"},             // \int_0^1 x\,dx
        {R"(["Add", 2, 2])", "2 + 2"},
    };
    for (const auto& [row, line] : rows) {
        INFO(row);
        auto t = bd::mathjson_row_text(json(row));
        REQUIRE(t);
        CHECK(*t == line);
    }
    // the second route: the first three rows as a document are the document written by hand, field by field
    std::string text;
    for (std::size_t k = 0; k < 3; ++k) text += *bd::mathjson_row_text(json(rows[k].first)) + "\n";
    auto doc = bd::read_document(text);
    REQUIRE(doc);
    auto by_hand = bd::read_document("p in [0.5, 4]\nx: x^2 - p = 0 on [0, 10]\nI = integral exp(-p * t^2) dt on [0, 1]");
    REQUIRE(by_hand);
    CHECK(bd::same_document(*doc, *by_hand));

    // the cases it must fail, each with where it stood
    const auto refused = [](const std::string& row) {
        auto t = bd::mathjson_row_text(json(row));
        INFO(row);
        REQUIRE_FALSE(t);
        return t.error().message;
    };
    CHECK(refused(R"(["Element", "Pi", ["Interval", 0, 1]])") == "at [1]: 'Pi' is a constant, not a name to define");
    CHECK(refused(R"(["Element", "p", ["Interval", ["Open", 0], 1]])").starts_with("at [2][1]: an open end"));
    CHECK(refused(R"(["Element", "p", ["List", 0, 1]])") == "at [2]: an interval is written [a, b]");
    CHECK(refused(R"(["Delimiter", ["Sequence", ["Element", "x", ["Interval", 0, 1]], ["Element", "y", ["Interval", 0, 1]]], {"str": ","}])")
              .starts_with("at [1]: a row holds one statement"));
    CHECK(refused(R"(["Integrate", "x", ["Tuple", "x"]])") == "at [2]: an integral here has both ends, \\int_a^b");
    CHECK(refused(R"(["Equal", ["Add", "p", 1], ["Integrate", "x", ["Tuple", "x", 0, 1]]])") == "at [1]: a name is a symbol, such as x");
}

TEST_CASE("what the language cannot say is refused by name and place", "[io][mathjson]") {
    const auto has = [](const std::string& message, const std::string& part) {
        INFO(message);
        CHECK(message.find(part) != std::string::npos);
    };
    has(refusal(R"(["Integrate", "p"])"), "'Integrate' is not in the language");
    has(refusal(R"(["Add", "p", ["Sin", ["D", "p"]]])"), "at [2][1]: 'D' is not in the language");
    has(refusal(R"(["Add", "p", ["Sin", ["D", "p"]]])"), "document: a, value");
    has(refusal(R"(["Error", "'missing'"])"), "the editor could not read the formula ('missing')");
    has(refusal(R"(["Multiply", "p", {"sym": "ImaginaryUnit"}])"), "at [2]: 'ImaginaryUnit' is a complex number");
    has(refusal(R"(["Sqrt", "p", "p"])"), "'Sqrt' takes 1 argument, not 2");
    has(refusal(R"(["Add", "p"])"), "'Add' takes 2 or more arguments, not 1");
    has(refusal(R"j({"num": "1.(3)"})j"), "'1.(3)' is not a number the language writes");
    has(refusal(R"({"sym": "2x"})"), "the symbol '2x' is not a name the language reads");
    has(refusal(R"(["Equal", "p", 1])"), "an equation stands only as a whole formula");
    has(refusal(R"({"what": 1})"), "not {\"num\"}, {\"sym\"} or {\"fn\"}");

    // read as text, the same refusals the language has: a name not defined above
    has(refusal(R"(["Add", "p", "q"])"), "q");
}
