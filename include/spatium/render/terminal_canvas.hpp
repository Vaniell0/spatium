#pragma once

// Presents a computed RGB frame (raytracer/software-rendered output) into
// a terminal, live/animated -- the terminal counterpart to write_image.hpp
// (which targets a file). Two problems a naive per-cell std::cout loop
// gets wrong, both handled here:
//
//   - A terminal character cell is not square (commonly ~1:2, width:height),
//     so one color per cell stretches the image vertically. Fixing this
//     doesn't need measuring the font: splitting each cell into two
//     independently-colored halves via the U+2580 upper-half-block glyph
//     turns each cell into two roughly-square sub-pixels stacked
//     vertically -- the standard technique behind every ANSI-truecolor
//     terminal image viewer.
//   - Writing escape codes straight to std::cout one cell at a time causes
//     visible tearing during animation. write_terminal_frame appends into
//     one std::string instead, so the caller flushes a whole frame as a
//     single write.

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/algebra/vector.hpp>
#  include <algorithm>
#  include <cstdio>
#  include <string>
#  include <string_view>
#  include <vector>
#  ifdef __unix__
#    include <sys/ioctl.h>
#    include <unistd.h>
#  endif
#endif

SPATIUM_EXPORT namespace spatium::render {

struct TerminalSize {
    int cols = 100;
    int rows = 40;
    // width/height of one character cell. ~0.5 is a typical monospace
    // guess; query_terminal_size() replaces it with an exact value
    // whenever the terminal reports pixel geometry.
    double cell_aspect = 0.5;
};

// Queries the real terminal size via ioctl(TIOCGWINSZ), including the
// exact character-cell aspect ratio when the terminal fills in the pixel
// geometry fields (xterm, kitty, alacritty, foot do; some multiplexers
// report zero there). Returns `fallback` unchanged on a non-POSIX build,
// a non-TTY stdout, or when the ioctl itself fails.
inline TerminalSize query_terminal_size(const TerminalSize& fallback = {}) {
#ifdef __unix__
    struct winsize w{};
    if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &w) == 0 && w.ws_col > 0 && w.ws_row > 0) {
        TerminalSize size;
        size.cols = w.ws_col;
        size.rows = w.ws_row;
        if (w.ws_xpixel > 0 && w.ws_ypixel > 0) {
            double cell_w = static_cast<double>(w.ws_xpixel) / static_cast<double>(w.ws_col);
            double cell_h = static_cast<double>(w.ws_ypixel) / static_cast<double>(w.ws_row);
            size.cell_aspect = cell_w / cell_h;
        } else {
            size.cell_aspect = fallback.cell_aspect;
        }
        return size;
    }
#endif
    return fallback;
}

// RAII cursor hide/restore. Plain "print the show-cursor code after the
// render loop" doesn't run on Ctrl+C or an early return out of an
// infinite loop -- this restores from the destructor instead, so the
// terminal cursor comes back regardless of how the caller exits.
class TerminalCursorGuard {
public:
    TerminalCursorGuard() {
        std::fputs("\033[?25l", stdout);
        std::fflush(stdout);
    }
    ~TerminalCursorGuard() {
        std::fputs("\033[?25h", stdout);
        std::fflush(stdout);
    }
    TerminalCursorGuard(const TerminalCursorGuard&) = delete;
    TerminalCursorGuard& operator=(const TerminalCursorGuard&) = delete;
};

namespace detail {

inline void append_color_code(std::string& out, const char* prefix,
                               const spatium::Vec<double, 3>& c) {
    int r = std::clamp(static_cast<int>(c[0]), 0, 255);
    int g = std::clamp(static_cast<int>(c[1]), 0, 255);
    int b = std::clamp(static_cast<int>(c[2]), 0, 255);
    out += prefix;
    out += std::to_string(r);
    out += ';';
    out += std::to_string(g);
    out += ';';
    out += std::to_string(b);
    out += 'm';
}

} // namespace detail

// Appends one rendered frame to `out` as half-block ANSI truecolor rows.
// `frame` is row-major, RGB in [0,255] per the render/color.hpp
// convention (frame[row][col]). Rows are consumed two at a time -- an
// odd row count renders its last row's bottom half as default background.
inline void write_terminal_frame(std::string& out,
                                  const std::vector<std::vector<spatium::Vec<double, 3>>>& frame) {
    out += "\033[H";
    std::size_t height = frame.size();
    for (std::size_t j = 0; j < height; j += 2) {
        const auto& top = frame[j];
        const std::vector<spatium::Vec<double, 3>>* bottom =
            (j + 1 < height) ? &frame[j + 1] : nullptr;
        std::size_t width = top.size();
        for (std::size_t i = 0; i < width; ++i) {
            detail::append_color_code(out, "\033[38;2;", top[i]); // fg = top sub-pixel
            if (bottom != nullptr && i < bottom->size())
                detail::append_color_code(out, "\033[48;2;", (*bottom)[i]); // bg = bottom sub-pixel
            else
                out += "\033[49m"; // no bottom half this row: default bg
            out += "\xE2\x96\x80"; // U+2580 UPPER HALF BLOCK
        }
        out += "\033[0m\n";
    }
}

// Classic dim -> bright luminance ramp (Bresenham/a1k0n ASCII-donut
// lineage). Any caller-supplied palette works just as well --
// write_terminal_frame_ascii doesn't assume this one.
inline constexpr std::string_view kDefaultAsciiPalette = ".,-~:;=!*#$@";

// Appends one rendered frame to `out` as plain ASCII, one character per
// cell -- no color, no half-block packing. `intensity` is row-major;
// each value is clamped into [0,1] and mapped to the nearest character
// in `palette` (index 0 = dimmest, last = brightest). A negative value
// means "nothing drawn here" and renders as a literal space, distinct
// from the darkest palette character. `palette` can be any ordered
// dim->bright character set -- swap it for a blockier or sparser look,
// there's nothing special about the default.
inline void write_terminal_frame_ascii(std::string& out,
                                        const std::vector<std::vector<double>>& intensity,
                                        std::string_view palette = kDefaultAsciiPalette) {
    out += "\033[H";
    for (const auto& row : intensity) {
        for (double v : row) {
            if (v < 0.0) {
                out += ' ';
                continue;
            }
            double clamped = std::clamp(v, 0.0, 1.0);
            std::size_t idx = static_cast<std::size_t>(
                clamped * static_cast<double>(palette.size() - 1) + 0.5);
            idx = std::min(idx, palette.size() - 1);
            out += palette[idx];
        }
        out += '\n';
    }
}

} // namespace spatium::render
