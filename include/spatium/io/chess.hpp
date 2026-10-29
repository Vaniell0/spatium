#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#  include <spatium/core/error.hpp>
#  include <spatium/discrete/position.hpp>
#  include <array>
#  include <charconv>
#  include <cstdint>
#  include <cstdlib>
#  include <optional>
#  include <string>
#  include <string_view>
#endif

// ── Chess notation ────────────────────────────────────────────────
//
// The text forms of discrete/position.hpp's positions and moves: FEN for a
// position, `e2e4` / `e7e8q` for a move as engines exchange it, SAN (`Nf3`,
// `exd5`, `O-O`, `e8=Q+`) as people write it. Reading text is this file's
// job; whether the position it describes may stand is decided by
// `Position::from_setup`, which every FEN read goes through.
//
// The en passant square is written after every double push, as the original
// FEN does, not only when a capture is possible.

SPATIUM_EXPORT namespace spatium::io {

constexpr char piece_letter(discrete::PieceType t) {
    return "pnbrqk?"[static_cast<int>(t)];
}

inline std::string square_name(discrete::Square s) {
    return std::string{static_cast<char>('a' + s.file()), static_cast<char>('1' + s.rank())};
}

inline std::optional<discrete::Square> parse_square(std::string_view text) {
    if (text.size() != 2 || text[0] < 'a' || text[0] > 'h' || text[1] < '1' || text[1] > '8')
        return std::nullopt;
    return discrete::Square::at(text[0] - 'a', text[1] - '1');
}

// `e2e4`, or `e7e8q` for a promotion.
inline std::string to_uci(const discrete::Move& m) {
    std::string s = square_name(m.from) + square_name(m.to);
    if (m.promotion != discrete::PieceType::None)
        s += piece_letter(m.promotion);
    return s;
}

// ── FEN ───────────────────────────────────────────────────────────

inline std::string to_fen(const discrete::Position& pos) {
    using namespace discrete;
    std::string out;
    for (int rank = 7; rank >= 0; --rank) {
        int empty = 0;
        for (int file = 0; file < 8; ++file) {
            const Square s = Square::at(file, rank);
            const PieceType t = pos.piece_on(s);
            if (t == PieceType::None) {
                ++empty;
                continue;
            }
            if (empty > 0)
                out += static_cast<char>('0' + empty);
            empty = 0;
            const char letter = piece_letter(t);
            out += pos.color_on(s) == Color::White ? static_cast<char>(letter - 'a' + 'A') : letter;
        }
        if (empty > 0)
            out += static_cast<char>('0' + empty);
        if (rank > 0)
            out += '/';
    }
    out += pos.side_to_move() == Color::White ? " w " : " b ";
    const std::uint8_t castling = pos.castling_rights();
    if (castling == 0)
        out += '-';
    if (castling & white_kingside)
        out += 'K';
    if (castling & white_queenside)
        out += 'Q';
    if (castling & black_kingside)
        out += 'k';
    if (castling & black_queenside)
        out += 'q';
    out += ' ';
    out += pos.en_passant() ? square_name(*pos.en_passant()) : std::string("-");
    out += ' ' + std::to_string(pos.halfmove_clock()) + ' ' + std::to_string(pos.fullmove_number());
    return out;
}

// Four fields at least (placement, side, castling, en passant); the two clocks default to
// 0 and 1. A text that parses but describes no playable position fails with from_setup's
// reason.
inline Result<discrete::Position> from_fen(std::string_view fen) {
    using namespace discrete;
    const auto fail = [](std::string message) -> Result<Position> {
        return std::unexpected(Error{ErrorCode::ParseError, "FEN: " + std::move(message)});
    };

    std::array<std::string_view, 6> field{};
    std::size_t count = 0;
    for (std::size_t i = 0; i < fen.size() && count < 6;) {
        while (i < fen.size() && fen[i] == ' ')
            ++i;
        const std::size_t start = i;
        while (i < fen.size() && fen[i] != ' ')
            ++i;
        if (i > start)
            field[count++] = fen.substr(start, i - start);
    }
    if (count < 4)
        return fail("expected at least four fields (placement, side, castling, en passant)");

    PositionSetup setup;

    // Placement, from rank 8 down to rank 1.
    int rank = 7, file = 0;
    for (char ch : field[0]) {
        if (ch == '/') {
            if (file != 8 || rank == 0)
                return fail("a rank does not have eight squares");
            --rank;
            file = 0;
        } else if (ch >= '1' && ch <= '8') {
            file += ch - '0';
            if (file > 8)
                return fail("a rank has more than eight squares");
        } else {
            const bool white = ch >= 'A' && ch <= 'Z';
            const char lower = white ? static_cast<char>(ch - 'A' + 'a') : ch;
            const std::size_t t = std::string_view("pnbrqk").find(lower);
            if (t == std::string_view::npos)
                return fail(std::string("unknown piece '") + ch + "'");
            if (file > 7)
                return fail("a rank has more than eight squares");
            setup.squares[static_cast<std::size_t>(Square::at(file, rank).index)] =
                PositionSetup::Occupant{white ? Color::White : Color::Black,
                                        static_cast<PieceType>(t)};
            ++file;
        }
    }
    if (rank != 0 || file != 8)
        return fail("the placement does not cover all eight ranks");

    if (field[1] == "w")
        setup.side = Color::White;
    else if (field[1] == "b")
        setup.side = Color::Black;
    else
        return fail("the side to move is not w or b");

    if (field[2] != "-") {
        for (char ch : field[2]) {
            const std::uint8_t bit = ch == 'K'   ? white_kingside
                                     : ch == 'Q' ? white_queenside
                                     : ch == 'k' ? black_kingside
                                     : ch == 'q' ? black_queenside
                                                 : std::uint8_t{0};
            if (bit == 0 || (setup.castling & bit))
                return fail("the castling field is not a set of K, Q, k, q");
            setup.castling |= bit;
        }
    }

    if (field[3] != "-") {
        const auto ep = parse_square(field[3]);
        if (!ep)
            return fail("the en passant square is not a square");
        setup.en_passant = ep;
    }

    const auto number = [](std::string_view text, int& value) {
        const auto [end, ec] = std::from_chars(text.data(), text.data() + text.size(), value);
        return ec == std::errc{} && end == text.data() + text.size();
    };
    if (count >= 5 && !number(field[4], setup.halfmove))
        return fail("the halfmove clock is not a number");
    if (count >= 6 && !number(field[5], setup.fullmove))
        return fail("the move number is not a number");

    return Position::from_setup(setup);
}

// ── Standard algebraic notation ───────────────────────────────────

// `Nf3`, `exd5`, `O-O`, `e8=Q+`, `Qd1#`: the piece, as much of its square as tells it from the
// other pieces of its kind that reach the same square, `x` for a capture, `+` or `#` for a check
// that is not or is mate.
inline std::string to_san(const discrete::Position& pos, const discrete::Move& m) {
    using namespace discrete;
    const PieceType piece = pos.piece_on(m.from);
    std::string out;
    if (piece == PieceType::King && std::abs(m.to.file() - m.from.file()) == 2) {
        out = m.to.file() > m.from.file() ? "O-O" : "O-O-O";
    } else {
        const bool capture =
            pos.piece_on(m.to) != PieceType::None ||
            (piece == PieceType::Pawn && m.from.file() != m.to.file());
        if (piece == PieceType::Pawn) {
            if (capture)
                out += static_cast<char>('a' + m.from.file());
        } else {
            out += static_cast<char>(piece_letter(piece) - 'a' + 'A');
            bool other = false, same_file = false, same_rank = false;
            for (const Move& o : legal_moves(pos)) {
                if (o.to != m.to || o.from == m.from || pos.piece_on(o.from) != piece)
                    continue;
                other = true;
                same_file |= o.from.file() == m.from.file();
                same_rank |= o.from.rank() == m.from.rank();
            }
            if (other) {
                if (!same_file)
                    out += static_cast<char>('a' + m.from.file());
                else if (!same_rank)
                    out += static_cast<char>('1' + m.from.rank());
                else
                    out += square_name(m.from);
            }
        }
        if (capture)
            out += 'x';
        out += square_name(m.to);
        if (m.promotion != PieceType::None)
            out += std::string("=") + static_cast<char>(piece_letter(m.promotion) - 'a' + 'A');
    }
    const Position next = pos.made(m);
    if (next.in_check())
        out += legal_moves(next).size ? '+' : '#';
    return out;
}

// The legal move a SAN text names there; `+`, `#`, `!` and `?` are ignored, and so are a
// promotion's `=` and the `0` of a castling written with zeros.
inline std::optional<discrete::Move> parse_san(const discrete::Position& pos, std::string_view text) {
    const auto clean = [](std::string_view t) {
        std::string out;
        for (char ch : t)
            if (ch != '+' && ch != '#' && ch != '!' && ch != '?' && ch != '=' && ch != ' ')
                out += ch == '0' ? 'O' : ch;
        return out;
    };
    const std::string want = clean(text);
    for (const discrete::Move& m : discrete::legal_moves(pos))
        if (clean(to_san(pos, m)) == want)
            return m;
    return std::nullopt;
}

} // namespace spatium::io
