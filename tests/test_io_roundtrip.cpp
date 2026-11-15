// One statement, every way it is written and read: each writer's output read back is the statement
// it started from. A statement reaches the library as a line of a document, as the document's JSON,
// as a shorthand a calculator takes (`respond`), and leaves it as text, as JSON and as the LaTeX row
// a page shows. Where two of these disagree, a question means one thing typed and another
// answered; each case below is one kind of statement in each form it can take.

#include <spatium/answer/respond.hpp>
#include <spatium/io/document.hpp>
#include <spatium/io/json.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <algorithm>
#include <iterator>
#include <optional>
#include <ranges>
#include <string>
#include <vector>

namespace bd = spatium::io::build;

namespace {

struct Case {
    const char* above;   // lines it reads, defined first
    const char* line;    // the statement, as a document line
    const char* typed;   // the same statement as a calculator's shorthand; empty when it has none
};

// every kind of statement, and every form a root takes
const Case cases[] = {
    {"", "r = 0.5", ""},
    {"", "a = 2(3 + 4)^2 - 1/4", ""},
    {"", "answer = 2 + 2", "2 + 2"},
    {"", "answer = 3 * 4", "3 * 4 ="},
    {"a = 2\n", "f = a^2 + 1", "f = a^2 + 1"},
    {"", "p in [0.5, 4]", ""},
    {"p in [0.5, 4]\n", "x: x^2 - p = 0 on [0, 10]", "x^2 - p = 0 on [0, 10]"},
    {"", "y: 2y^2 - 3y - 5 = 0", "2y^2 - 3y - 5 = 0"},
    {"", "x: x^2 - 4 > 0", "x^2 - 4 > 0"},
    {"", "x: x^2 - 4 >= 0 on [-5, 5]", "x^2 - 4 >= 0 on [-5, 5]"},
    {"", "x, y: x + 2y = 5, 3x - y = 1", "x + 2y = 5\n3x - y = 1"},
    {"", "x, y: x = 2y, y = 1", "x = 2y\ny = 1"},   // x reads a name below it: a system, not two definitions
    {"p in [0.5, 4]\n", "I = integral exp(-p * t^2) dt on [0, 1]", ""},
    {"r = 0.5\n", "S0 = torus(2, r)", ""},
    {"", "answer = d/dx sin(x)^2", "d/dx sin(x)^2"},
    {"", "answer = d/dx x^3 at 2", "d/dx x^3 at 2"},
    {"p in [0.5, 4]\n", "g = d/dp exp(-p^2)", ""},
};

bd::Document read(const Case& c) {
    auto doc = bd::read_document(std::string(c.above) + c.line);
    REQUIRE(doc.has_value());
    return std::move(*doc);
}

// the document an answer says it read
std::optional<bd::Document> answered_document(const spatium::io::JsonValue& answer) {
    const spatium::io::JsonValue* d = answer.find("document");
    if (!d) return std::nullopt;
    auto doc = bd::document_from_json(*d);
    if (!doc) return std::nullopt;
    return std::move(*doc);
}

std::size_t count(const std::string& s, std::string_view what) {
    std::size_t n = 0;
    for (auto at = s.find(what); at != std::string::npos; at = s.find(what, at + what.size())) ++n;
    return n;
}

} // namespace

TEST_CASE("round trip: a statement read back from its text is the statement", "[io][roundtrip]") {
    const Case& c = GENERATE(from_range(std::begin(cases), std::end(cases)));
    INFO(c.line);
    const auto doc = read(c);
    const auto again = bd::read_document(bd::to_text(doc));
    REQUIRE(again.has_value());
    CHECK(bd::same_document(doc, *again));
    CHECK(bd::to_text(*again) == bd::to_text(doc));
}

TEST_CASE("round trip: a statement read back from its JSON is the statement", "[io][roundtrip]") {
    const Case& c = GENERATE(from_range(std::begin(cases), std::end(cases)));
    INFO(c.line);
    const auto doc = read(c);
    const auto direct = bd::document_from_json(bd::to_json(doc));
    REQUIRE(direct.has_value());
    CHECK(bd::same_document(doc, *direct));
    const auto dumped = spatium::io::parse(bd::to_json(doc).dump());
    REQUIRE(dumped.has_value());
    const auto through = bd::document_from_json(*dumped);
    REQUIRE(through.has_value());
    CHECK(bd::same_document(doc, *through));
}

TEST_CASE("round trip: the LaTeX rows of a statement hold all of it", "[io][roundtrip]") {
    const Case& c = GENERATE(from_range(std::begin(cases), std::end(cases)));
    INFO(c.line);
    const auto doc = read(c);
    const bd::Statement& s = doc.statements.back();
    const auto rows = bd::to_latex_rows(s);
    REQUIRE(rows.has_value());
    // a shape is a scene's, and has no row; a system an equation a row; anything else one row
    CHECK(rows->size() == (s.kind == bd::StatementKind::Shape ? 0 : 1 + s.also.size()));
    for (const auto& row : *rows) {
        INFO(row);
        if (s.kind == bd::StatementKind::Root && s.relation == "=") CHECK(count(row, "=") == 1);
    }
    // one row alone never stands for a system: it is refused, not written as its first equation
    if (!s.also.empty()) CHECK_FALSE(bd::to_latex(s).has_value());
}

TEST_CASE("round trip: the document an answer says it read is the one typed", "[io][roundtrip]") {
    const Case& c = GENERATE(from_range(std::begin(cases), std::end(cases)));
    INFO(c.line);
    const auto doc = read(c);
    const auto answer = spatium::answer::respond(std::string(c.above) + c.line);
    INFO(answer.dump());
    const auto said = answered_document(answer);
    REQUIRE(said.has_value());
    CHECK(bd::same_document(doc, *said));
}

TEST_CASE("round trip: a shorthand is read as the statement it stands for", "[io][roundtrip]") {
    const Case& c = GENERATE(from_range(std::begin(cases), std::end(cases)));
    if (std::string_view(c.typed).empty()) return;
    INFO(c.typed);
    const auto doc = read(c);
    const auto answer = spatium::answer::respond(std::string(c.above) + c.typed);
    INFO(answer.dump());
    const auto said = answered_document(answer);
    REQUIRE(said.has_value());
    CHECK(bd::same_document(doc, *said));
}

TEST_CASE("round trip: every form a statement is written in reads one of the cases", "[io][roundtrip]") {
    // the cases above are complete when each form of the table reads one of their lines whole
    for (const bd::StatementForm& form : bd::statement_forms) {
        INFO(form.line);
        const bool read_one = std::ranges::any_of(cases, [&](const Case& c) {
            std::vector<std::string_view> lines{c.line};
            for (const auto l : std::views::split(std::string_view(c.typed), '\n')) lines.emplace_back(l.begin(), l.end());
            return std::ranges::any_of(lines, [&](std::string_view l) { return !l.empty() && bd::detail::read_form(form, l).statement.has_value(); });
        });
        CHECK(read_one);
    }
}

TEST_CASE("round trip: the rows an answer writes for its question read back as it", "[io][roundtrip]") {
    const Case& c = GENERATE(from_range(std::begin(cases), std::end(cases)));
    INFO(c.line);
    const auto doc = read(c);
    const auto answer = spatium::answer::respond(std::string(c.above) + c.line);
    const spatium::io::JsonValue* rows = answer.find("question_latex");
    if (doc.statements.back().kind == bd::StatementKind::Shape) return;
    REQUIRE(rows);
    REQUIRE(rows->is_array());
    // one row a statement, and a system one row more for each unknown past the first
    std::size_t expected = 0;
    for (const auto& s : doc.statements) expected += s.kind == bd::StatementKind::Shape ? 0 : 1 + s.also.size();
    CHECK(rows->as_array().size() == expected);
}
