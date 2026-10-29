#pragma once
// What the chess examples share: a move as a row of numbers. The first block is the geometry of the
// move from the piece's vectors (discrete/move_features.hpp); the other six are this program's idea
// of what a move is worth -- what it captures, whether it checks or mates, the piece it leaves
// where it can be taken, what it promotes to, the value of the piece that moves.

#include <spatium/discrete/move_features.hpp>
#include <spatium/discrete/position.hpp>

#include <array>
#include <vector>

namespace chess_rows {

using namespace spatium::discrete;

inline constexpr std::size_t geometric = geometric_feature_count<2>;
inline constexpr std::size_t columns = geometric + 6;
using Row = std::array<double, columns>;

inline double value_of(PieceType t) {
    switch (t) {
    case PieceType::Pawn:
        return 1;
    case PieceType::Knight:
    case PieceType::Bishop:
        return 3;
    case PieceType::Rook:
        return 5;
    case PieceType::Queen:
        return 9;
    default:
        return 0;
    }
}

// The piece as the vectors it moves by. A pawn has a direction, so it is built by hand.
inline Piece<2> piece_for(PieceType t, Color c) {
    switch (t) {
    case PieceType::Knight:
        return knight<2>();
    case PieceType::Bishop:
        return bishop<2>();
    case PieceType::Rook:
        return rook<2>();
    case PieceType::Queen:
        return queen<2>();
    case PieceType::King:
        return king<2>();
    default: {
        const int up = c == Color::White ? 1 : -1;
        return Piece<2>::of(
            {{{0, up}, Reach::Leap},
             {{0, 2 * up}, Reach::Leap},
             {{1, up}, Reach::Leap},
             {{-1, up}, Reach::Leap}}
        );
    }
    }
}

// What the move takes: the piece on the destination, or the pawn passed by en passant.
inline double captured(const Position& p, const Move& m) {
    const PieceType victim = p.piece_on(m.to);
    if (victim != PieceType::None)
        return value_of(victim);
    const bool en_passant = p.piece_on(m.from) == PieceType::Pawn && m.from.file() != m.to.file();
    return en_passant ? 1.0 : 0.0;
}

struct Judged {
    double captured, check, mate, hanging, promoted, moved;
};

inline Judged judge(const Position& p, const Move& m) {
    const Color mover = p.side_to_move();
    const PieceType piece = p.piece_on(m.from);
    const Position after = p.made(m);
    Judged j{};
    j.captured = captured(p, m);
    j.check = after.in_check() ? 1.0 : 0.0;
    j.mate = after.in_check() && legal_moves(after).size == 0 ? 1.0 : 0.0;
    const PieceType landed = m.promotion != PieceType::None ? m.promotion : piece;
    j.hanging = after.attacked(m.to, opposite(mover)) ? value_of(landed) : 0.0;
    j.promoted = m.promotion != PieceType::None ? value_of(m.promotion) : 0.0;
    j.moved = value_of(piece);
    return j;
}

inline Row row_of(const Position& p, const Move& m) {
    static const Shape<2> board{{8, 8}};
    const Color mover = p.side_to_move();
    const PieceType piece = p.piece_on(m.from);
    const auto geo = geometric_features(
        piece_for(piece, mover), board, {m.from.file(), m.from.rank()}, {m.to.file(), m.to.rank()}
    );
    const Judged j = judge(p, m);
    const std::array<double, 6> game{j.captured / 9, j.check,        j.mate,
                                     j.hanging / 9,  j.promoted / 9, j.moved / 9};
    return join(geo, game);
}

inline std::vector<Row> rows_of(const Position& p, const MoveList& moves) {
    std::vector<Row> rows;
    rows.reserve(static_cast<std::size_t>(moves.size));
    for (const Move& m : moves)
        rows.push_back(row_of(p, m));
    return rows;
}

inline double material(const Position& p, Color c) {
    double total = 0;
    for (PieceType t : all_piece_types)
        total += value_of(t) * p.pieces(c, t).count();
    return total;
}

} // namespace chess_rows
