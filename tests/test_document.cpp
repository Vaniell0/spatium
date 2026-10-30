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

TEST_CASE("ask: a question read through its document is the question read word by word", "[document][rsc]") {
    // every question the gallery and the reference ask, and the forms around them
    const std::vector<std::string> questions{
        "solve x^2 - p = 0 for x in [0, 10], p in [0.5, 4]",
        "solve x^2 + p = 0 for x in [0, 10], p in [0.5, 4]",
        "solve x^3 = p for x in [0, 3], p in [1, 8]",
        "solve x - cos(p * x) for x in [0, 2], p in [0.5, 1]",
        "integrate exp(p * log(x)) for x in [0, 1], p in [0, 3]",
        "integrate exp(-p * x^2) for x in [0, 1], p in [0.1, 4]",
        "integrate sin(t)/t for t in [0, x], x in [0.5, 6]",
        "integrate 2 / sqrt(pi) * exp(-t^2) for t in [0, x], x in [0.25, 3]",
        "integrate cos(pi * t^2 / 2) for t in [0, x], x in [0.25, 2]",
        "integrate t for t in [p, 2 p], p in [1, 2]",
    };
    const auto same_bound = [](const std::optional<rsc::Expr>& a, const std::optional<rsc::Expr>& b) {
        return a.has_value() == b.has_value() && (!a || bd::same_value(*a, *b));
    };
    for (const auto& q : questions) {
        INFO(q);
        const auto by_words = rsc::detail::ask_by_words(q);
        const auto by_document = rsc::ask(q);
        REQUIRE(by_words.has_value());
        REQUIRE(by_document.has_value());
        CHECK(by_words->ask == by_document->ask);
        CHECK(bd::same_value(by_words->f, by_document->f));
        CHECK(by_words->unknown == by_document->unknown);
        CHECK(by_words->parameter == by_document->parameter);
        CHECK(by_words->x_lo == by_document->x_lo);
        CHECK(by_words->x_hi == by_document->x_hi);
        CHECK(same_bound(by_words->x_lo_at, by_document->x_lo_at));
        CHECK(same_bound(by_words->x_hi_at, by_document->x_hi_at));
        CHECK(by_words->p_lo == by_document->p_lo);
        CHECK(by_words->p_hi == by_document->p_hi);
    }
    // both refuse what neither can ask
    for (const char* bad : {"solve x^2 - p = 0 for x in [0, 10]", "solve x^2 - p = 0 for x in [0, 10], x in [0, 4]",
                            "solve x^2 - p = 0 for x in [0, 10], p in [4, 1]", "integrate t for t in [p, 0], p in [1, 2]",
                            "find x"}) {
        INFO(bad);
        CHECK_FALSE(rsc::detail::ask_by_words(bad).has_value());
        CHECK_FALSE(rsc::ask(bad).has_value());
    }
    // a name nothing defines is refused when the question is read, where the word reader let it
    // through to the answer
    const auto open = rsc::ask("solve x^2 - q = 0 for x in [0, 10], p in [0.5, 4]");
    REQUIRE_FALSE(open.has_value());
    CHECK(open.error().message == "x reads 'q', which no line above defines");
    // an error names the column of what was typed, not of the document it became
    const auto e = rsc::ask("solve x^2 - = 0 for x in [0, 10], p in [0.5, 4]");
    REQUIRE_FALSE(e.has_value());
    CHECK(e.error().message.starts_with("parse: at column 13,"));
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
