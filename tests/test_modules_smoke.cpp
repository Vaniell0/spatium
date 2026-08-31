// Smoke-test consumer: a single `import` from spatium.core proves the
// modules helper + std.cppm chain is fully wired. The earlier
// stand-alone spatium.smoke module was a Phase-0 canary; once we
// have ten real modules importing each other it is redundant.

#include <catch2/catch_test_macros.hpp>

import spatium.core;

TEST_CASE("modules smoke: import spatium.core compiles", "[modules][smoke]") {
    using ::spatium::Error;
    using ::spatium::ErrorCode;
    Error e{ErrorCode::DegenerateInput, "smoke"};
    REQUIRE(e.code == ErrorCode::DegenerateInput);
}
