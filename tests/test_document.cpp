// A document of named definitions, read from text and from JSON, printed back as it was written,
// and the question a calculator reads as one: the same Query whether the question is read word by
// word or through its document. Its shapes are a scene: the same Trace, node for node and line for
// line of its LaTeX, whether built from the document or by calls.

#include <rsc/spatium/query.hpp>
#include <spatium/io/document.hpp>
#include <spatium/io/document_trace.hpp>
#include <spatium/io/json.hpp>
#include <spatium/io/trace_text.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <format>
#include <string>
#include <vector>

namespace bd = spatium::io::build;
using Catch::Matchers::WithinAbs;

namespace {

const char* const sample = "r = 0.5\n"
                           "A = pi * r^2\n"
                           "p in [0.5, 4]\n"
                           "x: x^2 - p = 0 on [0, 10]\n"
                           "I = integral exp(-p * t^2) dt on [0, 1]\n"
                           "S = integral sin(t)/t dt on [0, p]\n";

} // namespace

TEST_CASE("document: a definition a line, each kind read as written", "[document]") {
    const auto doc = bd::read_document(sample);
    REQUIRE(doc.has_value());
    REQUIRE(doc->statements.size() == 6);
    using K = bd::StatementKind;
    CHECK(doc->statements[0].kind == K::Value);
    CHECK(doc->statements[1].value.text == "pi * r^2");
    CHECK(doc->statements[2].kind == K::Parameter);
    CHECK(doc->statements[3].kind == K::Root);
    CHECK(doc->statements[3].name == "x");
    CHECK(doc->statements[4].kind == K::Integral);
    CHECK(doc->statements[4].variable == "t");
    CHECK(doc->statements[4].value.text == "exp(-p * t^2)");
    CHECK(doc->statements[5].hi.text == "p");
    // a value reads the names above it, and the field is what the text says
    CHECK(bd::same_value(doc->statements[1].value.field, *bd::parse("pi * r^2")));
}

TEST_CASE("document: printed back as written, and read back the same through text and JSON", "[document]") {
    const auto doc = bd::read_document("# a comment line\n\n  r   =  0.5   # the radius\n" + std::string(sample).substr(8));
    REQUIRE(doc.has_value());
    CHECK(bd::to_text(*doc) == sample);   // spaces around the definition and comments are not kept
    const auto again = bd::read_document(bd::to_text(*doc));
    REQUIRE(again.has_value());
    CHECK(bd::same_document(*doc, *again));
    // the JSON carries the same texts, read by the same reader; and it survives being written out
    const auto from_json = bd::document_from_json(bd::to_json(*doc));
    REQUIRE(from_json.has_value());
    CHECK(bd::same_document(*doc, *from_json));
    const auto reread = spatium::io::parse(bd::to_json(*doc).dump());
    REQUIRE(reread.has_value());
    const auto through_dump = bd::document_from_json(*reread);
    REQUIRE(through_dump.has_value());
    CHECK(bd::same_document(*doc, *through_dump));
    // the case it must fail: a different interval is a different document
    const auto other = bd::read_document(std::string(sample).replace(std::string(sample).find("[0, 10]"), 7, "[0, 11]"));
    REQUIRE(other.has_value());
    CHECK_FALSE(bd::same_document(*doc, *other));
}

TEST_CASE("document: a system names its unknowns, then an equation each", "[document]") {
    const char* text = "a = 2\nx, y: x + a*y = 5, 3*x - y = 1\n";
    const auto doc = bd::read_document(text);
    REQUIRE(doc.has_value());
    CHECK(bd::to_text(*doc) == text);
    // each unknown is found, and is the system's; a line below may read it
    REQUIRE(doc->find("y") != nullptr);
    CHECK(doc->find("y") == doc->find("x"));
    CHECK(doc->find("y")->also.front().value.text == "3*x - y");
    const auto from_json = bd::document_from_json(bd::to_json(*doc));
    REQUIRE(from_json.has_value());
    CHECK(bd::same_document(*doc, *from_json));
    // the cases it must fail: an equation short, a name twice, an unknown no equation reads, a
    // different equation a different document
    const auto message = [](const char* t) { const auto d = bd::read_document(t); return d ? std::string("(read)") : d.error().message; };
    CHECK(message("x, y: x + y = 1").starts_with("line 1: parse:"));
    CHECK(message("x, x: x = 1, x = 2") == "line 1: 'x' is defined twice");
    CHECK(message("x, y: x = 1, x = 2") == "line 1: no equation of the system reads y");
    const auto other = bd::read_document("a = 2\nx, y: x + a*y = 5, 3*x - y = 2\n");
    REQUIRE(other.has_value());
    CHECK_FALSE(bd::same_document(*doc, *other));
}

TEST_CASE("document: refused with the line, and the column for a parse error", "[document]") {
    const auto message = [](const char* text) {
        const auto d = bd::read_document(text);
        return d ? std::string("(read)") : d.error().message;
    };
    CHECK(message("p in [0, 1]\nx: x^2 - q = 0 on [0, 1]") == "line 2: x reads 'q', which no line above defines");
    CHECK(message("x: x^2 - p = 0 on [0, 1]\np in [0, 1]") == "line 1: x reads 'p', which no line above defines");
    CHECK(message("a = 1\na = 2") == "line 2: 'a' is defined twice");
    CHECK(message("sin = 2") == "line 1: 'sin' is a word of the language, not a name");
    CHECK(message("a = 1\np in [0, a]") == "line 2: an end of the parameter p reads 'a'; it is a number");
    CHECK(message("p in [0, 1]\nx: p^2 = 1 on [0, 1]") == "line 2: the equation of x does not read x");
    CHECK(message("I = integral exp(t) on [0, 1]").starts_with("line 1: parse:"));
    CHECK(message("p in [0, 1]\n\nx: x^2 - = 0 on [0, 1]") == "line 3: parse: at column 10, expected a number, a name, a function or '(', found '='");
    CHECK(message("x = 1 2") == "line 1: parse: at column 7, expected the end of the line, found '2'");
}

TEST_CASE("query_of: every question the gallery asks, written as a document", "[document][rsc]") {
    struct Asked { const char* document; const char* wanted; rsc::Ask ask; const char* unknown; const char* parameter;
                   double x_lo, x_hi, p_lo, p_hi; bool bound_reads_parameter; };
    for (const Asked& a : std::initializer_list<Asked>{
             {"p in [0.5, 4]\nx: x^2 - p = 0 on [0, 10]", "x", rsc::Ask::Solve, "x", "p", 0, 10, 0.5, 4, false},
             {"p in [0.5, 4]\nx: x^2 + p = 0 on [0, 10]", "x", rsc::Ask::Solve, "x", "p", 0, 10, 0.5, 4, false},
             {"p in [1, 8]\nx: x^3 = p on [0, 3]", "x", rsc::Ask::Solve, "x", "p", 0, 3, 1, 8, false},
             {"p in [0.5, 1]\nx: x - cos(p * x) = 0 on [0, 2]", "x", rsc::Ask::Solve, "x", "p", 0, 2, 0.5, 1, false},
             {"p in [0, 3]\nI = integral exp(p * log(x)) dx on [0, 1]", "I", rsc::Ask::Integrate, "x", "p", 0, 1, 0, 3, false},
             {"p in [0.1, 4]\nI = integral exp(-p * x^2) dx on [0, 1]", "I", rsc::Ask::Integrate, "x", "p", 0, 1, 0.1, 4, false},
             {"x in [0.5, 6]\nSi = integral sin(t)/t dt on [0, x]", "Si", rsc::Ask::Integrate, "t", "x", 0, 0, 0.5, 6, true},
             {"p in [1, 2]\nI = integral t dt on [p, 2 p]", "I", rsc::Ask::Integrate, "t", "p", 0, 0, 1, 2, true},
         }) {
        INFO(a.document);
        const auto doc = bd::read_document(a.document);
        REQUIRE(doc.has_value());
        const auto q = rsc::query_of(*doc, a.wanted);
        REQUIRE(q.has_value());
        CHECK(q->ask == a.ask);
        CHECK(q->unknown == a.unknown);
        CHECK(q->parameter == a.parameter);
        CHECK(q->p_lo == a.p_lo);
        CHECK(q->p_hi == a.p_hi);
        CHECK((q->x_lo_at.has_value() || q->x_hi_at.has_value()) == a.bound_reads_parameter);
        if (!a.bound_reads_parameter) {
            CHECK(q->x_lo == a.x_lo);
            CHECK(q->x_hi == a.x_hi);
        }
    }
    // what a document cannot ask is refused when it is read: a name nothing defines, an interval backwards
    const auto open = bd::read_document("p in [0.5, 4]\nx: x^2 - q = 0 on [0, 10]");
    REQUIRE_FALSE(open.has_value());
    CHECK(open.error().message.find("reads 'q', which no line above defines") != std::string::npos);
    const auto backwards = bd::read_document("p in [4, 1]\nx: x^2 - p = 0 on [0, 10]");
    REQUIRE(backwards.has_value());
    CHECK(rsc::query_of(*backwards, "x").error().message == "p in [4, 1] is empty");
    // and a question with no parameter left free is a number, not a function: no Query
    const auto fixed = bd::read_document("x: x^2 - 2 = 0 on [0, 10]");
    REQUIRE(fixed.has_value());
    CHECK_FALSE(rsc::query_of(*fixed, "x").has_value());
}

TEST_CASE("to_latex: a statement as the row of a formula editor", "[document]") {
    const auto doc = bd::read_document(
        "p in [0.5, 4]\nx: x^2 - p = 0 on [0, 10]\ny: 2y^2 - 3y - 5 = 0\na = 2 * p + 1\nI = integral sin(t)/t dt on [0, p]");
    REQUIRE(doc.has_value());
    std::vector<std::string> rows;
    for (const auto& d : doc->statements) {
        const auto row = bd::to_latex(d);
        REQUIRE(row.has_value());
        rows.push_back(*row);
    }
    for (const auto& r : rows) INFO(r);
    CHECK(rows[0] == "p \\in \\left[0.5, 4\\right]");
    CHECK(rows[1].ends_with(", x \\in \\left[0, 10\\right]"));
    CHECK(rows[1].starts_with("x^{2} - p = 0"));
    CHECK(rows[2].find("\\in") == std::string::npos);
    CHECK(rows[3].starts_with("a = "));
    CHECK(rows[4].starts_with("I = \\int_{0}^{p} "));
    CHECK(rows[4].ends_with("\\,dt"));
    // a shape is a scene's and has no row
    const auto scene = bd::read_document("S = sphere(1)");
    REQUIRE(scene.has_value());
    CHECK_FALSE(bd::to_latex(scene->statements.front()).has_value());
}

TEST_CASE("query_of: a document's root with the values it reads put in place", "[document][rsc]") {
    const auto doc = bd::read_document("a = 2\nb = a + 1\np in [1, 3]\nx: x^2 - b * p = 0 on [0, 10]\ny: y - x = 0 on [0, 1]");
    REQUIRE(doc.has_value());
    const auto q = rsc::query_of(*doc, "x");
    REQUIRE(q.has_value());
    CHECK(q->parameter == "p");
    CHECK_THAT(rsc::reference(rsc::bound(*q), 3.0), WithinAbs(3.0, 1e-12));   // x^2 = 3 * 3
    // the cases it must fail: a value is not a question, and a root that reads another root is
    // not one yet
    CHECK_FALSE(rsc::query_of(*doc, "a").has_value());
    CHECK_FALSE(rsc::query_of(*doc, "y").has_value());
    CHECK_FALSE(rsc::query_of(*doc, "z").has_value());
}

namespace {

// the gallery's shell, then a dusting of spheres on it, and both together
const char* const shell_document = "R = 2\n"
                                   "r = 0.5\n"
                                   "h: (r + h)^2 - 2 * r^2 = 0 on [0, 1]\n"
                                   "S0 = torus(R, r)\n"
                                   "S1 = offset(S0, h)\n"
                                   "S2 = sphere(r / 10)\n"
                                   "S3 = scatter(S2, S1, 40, 7)\n"
                                   "S4 = union(S1, S3)\n";

std::string scene_latex(const char* text, const std::vector<std::pair<std::string, double>>& chosen = {}) {
    const auto doc = bd::read_document(text);
    if (!doc) return doc.error().message;
    const auto scene = bd::to_trace(*doc, chosen);
    if (!scene) return scene.error().message;
    const auto latex = bd::to_latex(scene->trace, scene->root());
    return latex ? *latex : latex.error().message;
}

} // namespace

TEST_CASE("document: shapes read as written, printed back and read back through JSON", "[document]") {
    const auto doc = bd::read_document(shell_document);
    REQUIRE(doc.has_value());
    REQUIRE(doc->statements.size() == 8);
    const auto& offset = doc->statements[4];
    CHECK(offset.kind == bd::StatementKind::Shape);
    CHECK(offset.operation == "offset");
    REQUIRE(offset.args.size() == 2);
    CHECK(offset.args[0].text == "S0");
    CHECK(doc->statements[5].args[0].text == "r / 10");
    CHECK(bd::to_text(*doc) == shell_document);
    const auto from_json = bd::document_from_json(bd::to_json(*doc));
    REQUIRE(from_json.has_value());
    CHECK(bd::same_document(*doc, *from_json));
    // the case it must fail: another argument is another document
    const auto other = bd::read_document(std::string(shell_document).replace(std::string(shell_document).find("40"), 2, "41"));
    REQUIRE(other.has_value());
    CHECK_FALSE(bd::same_document(*doc, *other));
}

TEST_CASE("document: a shape refused by what is wrong with it", "[document]") {
    const auto message = [](const char* text) {
        const auto d = bd::read_document(text);
        return d ? std::string("(read)") : d.error().message;
    };
    CHECK(message("S = torus(2)") == "line 1: torus takes 2 arguments, given 1");
    CHECK(message("S = flake(1, 1)") == "line 1: flake takes 3 to 4 arguments, given 2");
    CHECK(message("S = union()") == "line 1: union takes 1 or more arguments, given 0");
    CHECK(message("r = 1\nS = offset(r, 0.1)") == "line 2: argument 1 of offset is a shape, given 'r'");
    CHECK(message("S = sphere(1)\nT = sphere(S)") == "line 2: argument 1 of sphere reads the shape S as a number");
    CHECK(message("S = sphere(1)\na = S + 1") == "line 2: a reads the shape S as a number");
    CHECK(message("S = sphere(1)\nT = offset(S, q)") == "line 2: T reads 'q', which no line above defines");
    // every operation's name is a word of the language
    for (const bd::ShapeOperation& op : bd::shape_operations)
        CHECK(message(std::format("{} = 1", op.name).c_str()) == std::format("line 1: '{}' is a word of the language, not a name", op.name));
    CHECK(message("S = sphere 1") == "line 1: parse: at column 12, expected '(' after sphere, found '1'");
    CHECK(message("S = sphere(1") == "line 1: parse: at column 13, expected ',' or ')', found the end");
}

TEST_CASE("to_trace: a document's scene is the scene its calls build", "[document][build_dsl]") {
    const auto doc = bd::read_document(shell_document);
    REQUIRE(doc.has_value());
    const auto scene = bd::to_trace(*doc);
    REQUIRE(scene.has_value());

    // the same scene by calls, the equation as the document's reader reads it
    bd::Trace<double> t;
    const auto R = t.let("R", 2.0), r = t.let("r", 0.5);
    const auto h = t.let("h", *bd::parse("(r + h)^2 - 2 * r^2"), 0.0, 1.0);
    const auto shell = t.offset(t.torus(R, r), h.field());
    const auto dust = t.scatter(t.sphere(r / 10), shell, 40, 7);
    const auto both = t.compose({shell, dust});

    const auto& a = scene->trace;
    REQUIRE(a.size() == t.size());
    CHECK(scene->root() == both.index);
    for (std::size_t i = 0; i < t.size(); ++i) {
        INFO(i);
        const auto &x = a.node(i), &y = t.node(i);
        CHECK(x.kind == y.kind);
        CHECK(x.recipe == y.recipe);
        CHECK(x.base == y.base);
        CHECK(x.item == y.item);
        CHECK(x.target == y.target);
        CHECK(x.count == y.count);
        CHECK(x.seed == y.seed);
        CHECK(x.children == y.children);
    }
    REQUIRE(a.definitions().size() == t.definitions().size());
    for (std::size_t i = 0; i < t.definitions().size(); ++i)
        CHECK(a.definitions()[i].value == t.definitions()[i].value);   // to the bit: the same root found the same way
    const auto by_document = bd::to_latex(a, scene->root()), by_calls = bd::to_latex(t, both.index);
    REQUIRE(by_document.has_value());
    REQUIRE(by_calls.has_value());
    CHECK(*by_document == *by_calls);

    // the case it must fail: another tube is another recipe, and another line of LaTeX
    const auto wider = bd::read_document(std::string(shell_document).replace(std::string(shell_document).find("torus(R, r)"), 11, "torus(R, 2 * r)"));
    REQUIRE(wider.has_value());
    const auto other = bd::to_trace(*wider);
    REQUIRE(other.has_value());
    CHECK(other->trace.node(0).recipe != t.node(0).recipe);
    CHECK(*bd::to_latex(other->trace, other->root()) != *by_calls);
}

TEST_CASE("to_trace: a parameter is drawn at the value chosen for it, and what a scene cannot take is refused",
          "[document][build_dsl]") {
    const char* const text = "p in [0.5, 4]\nS0 = sphere(p)\n";
    const std::string latex = scene_latex(text, {{"p", 1.0}});
    CHECK(latex == "\\begin{aligned}\np &= 1,\\quad p \\in \\left[0.5, 4\\right] \\\\\nS_{0} &= \\operatorname{sphere}\\left(1\\right)\n\\end{aligned}");
    CHECK(scene_latex(text) == "p: a parameter on [0.5, 4], drawn at a value it is not given");
    CHECK(scene_latex(text, {{"p", 5.0}}) == "p: choose(p): the value is not on its interval");
    CHECK(scene_latex("I = integral t dt on [0, 1]\nS = sphere(I)") == "I: an integral, which is not a number a scene takes");
    CHECK(scene_latex("S = sphere(1)\nT = scatter(S, S, 2.5)") == "T: scatter: '2.5' is 2.5, not a whole number from 0 to 1000000000");
    CHECK(scene_latex("r = 1") == "the document has no shape");
    // a builder's own refusal, named by the line it came from: an open base has no offset
    CHECK(scene_latex("S = flake(1, 1, 0.1)\nT = offset(S, 0.1)").starts_with("T: offset: the base surface has an edge"));
}
