#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <algorithm>
#  include <format>
#  include <iostream>
#  include <string>
#  include <utility>
#  include <vector>
#endif

SPATIUM_EXPORT namespace spatium::io {

struct Table {
    std::string col1_header;
    std::string col2_header;
    std::vector<std::pair<std::string, std::string>> rows;

    inline Table(std::string h1, std::string h2)
        : col1_header(std::move(h1)), col2_header(std::move(h2)) {}

    template<typename T>
    Table& row(std::string key, const T& value) {
        rows.emplace_back(std::move(key), std::format("{}", value));
        return *this;
    }

    inline Table& row(std::string key, std::string value) {
        rows.emplace_back(std::move(key), std::move(value));
        return *this;
    }

    inline void print(std::ostream& os = std::cout) const {
        std::size_t w1 = col1_header.size();
        std::size_t w2 = col2_header.size();
        for (const auto& [k, v] : rows) {
            w1 = std::max(w1, k.size());
            w2 = std::max(w2, v.size());
        }
        w1 += 1; w2 += 1;

        // Content cell is "│ {:<w1} │" → w1 + 2 visible chars between separators.
        auto hr_top = std::format("\u250c{:\u2500>{}}\u252c{:\u2500>{}}\u2510",
                                  "", w1 + 2, "", w2 + 2);
        auto hr_mid = std::format("\u251c{:\u2500>{}}\u253c{:\u2500>{}}\u2524",
                                  "", w1 + 2, "", w2 + 2);
        auto hr_bot = std::format("\u2514{:\u2500>{}}\u2534{:\u2500>{}}\u2518",
                                  "", w1 + 2, "", w2 + 2);

        os << hr_top << '\n';
        os << std::format("\u2502 {:<{}} \u2502 {:<{}} \u2502",
                          col1_header, w1, col2_header, w2) << '\n';
        os << hr_mid << '\n';
        for (const auto& [k, v] : rows) {
            os << std::format("\u2502 {:<{}} \u2502 {:<{}} \u2502",
                              k, w1, v, w2) << '\n';
        }
        os << hr_bot << '\n';
    }
};

inline void section(const std::string& title, std::ostream& os = std::cout) {
    os << "\n\u2500\u2500 " << title << ' ';
    std::size_t pad = (title.size() < 60) ? 60 - title.size() : 0;
    for (std::size_t i = 0; i < pad; ++i) os << "\u2500";
    os << '\n';
}

} // namespace spatium::io
