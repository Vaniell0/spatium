// A document of named definitions, read from text and from JSON, printed back as it was written,
// and the question a calculator reads as one: the same Query whether the question is read word by
// word or through its document.

#include <rsc/spatium/query.hpp>
#include <spatium/io/document.hpp>
#include <spatium/io/json.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
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
