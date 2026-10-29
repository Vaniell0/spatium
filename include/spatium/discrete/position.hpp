#pragma once

#include <spatium/_export_macro.hpp>
#ifndef SPATIUM_BUILDING_MODULE
#include <spatium/core/error.hpp>
#include <spatium/discrete/attacks.hpp>
#include <spatium/discrete/bitboard.hpp>

#include <array>
#include <cstdint>
#include <algorithm>
#include <optional>
#include <span>
#include <string>
#include <vector>
#endif

// ── A position of the game on the 8x8 board ───────────────────────
//
// `Position` holds who stands where as six sets of squares by piece type and two
// by colour, plus the four things a placement does not say: whose turn it is, the
// castling rights, the square an en passant capture would land on, and the two
// clocks. `made(move)` gives the position after a move by copying -- a position
// is about a hundred bytes, so a move is not undone, the earlier position is just
// kept.
//
// The rules are complete: castling, en passant, promotion, check. `legal_moves`
// generates every move of the side to move and keeps those that do not leave its
// own king attacked; `perft` counts the leaves of the tree of legal moves, the
// number every rules implementation is checked against.
//
// What a position may be is decided here (`from_setup`); how one is written down
// -- FEN, SAN, `e2e4` -- is io/chess.hpp, which reads text into a `PositionSetup`
// and hands it to `from_setup`.

SPATIUM_EXPORT namespace spatium::discrete {

enum class PieceType : std::uint8_t { Pawn, Knight, Bishop, Rook, Queen, King, None };

inline constexpr std::array<PieceType, 6> all_piece_types{PieceType::Pawn,   PieceType::Knight,
                                                          PieceType::Bishop, PieceType::Rook,
                                                          PieceType::Queen,  PieceType::King};

// Castling rights as bits of one byte.
inline constexpr std::uint8_t white_kingside = 1, white_queenside = 2, black_kingside = 4,
                              black_queenside = 8;

struct Move {
    Square from;
    Square to;
    PieceType promotion = PieceType::None;

    constexpr bool operator==(const Move&) const = default;
};

// Everything a position says, before the rules have looked at it.
struct PositionSetup {
    struct Occupant {
        Color color;
        PieceType type;
        constexpr bool operator==(const Occupant&) const = default;
    };
    std::array<std::optional<Occupant>, square_count> squares{};
    Color side = Color::White;
    std::uint8_t castling = 0;
    std::optional<Square> en_passant;
    int halfmove = 0;
    int fullmove = 1;
};

class Position {
  public:
    // A position the rules can be played from: one king a side, at most sixteen pieces a
    // side, pawns off the first and last ranks, the side that just moved not left in check,
    // an en passant square on the rank a double push of the other side passes, the clocks
    // not negative and the move number at least 1.
    static Result<Position> from_setup(const PositionSetup& setup);
    static Position initial();

    Color side_to_move() const { return side_; }
    std::uint8_t castling_rights() const { return castling_; }
    std::optional<Square> en_passant() const { return en_passant_; }
    int halfmove_clock() const { return halfmove_; }
    int fullmove_number() const { return fullmove_; }

    Bitboard pieces(Color c, PieceType t) const {
        return types_[static_cast<std::size_t>(t)] & colors_[static_cast<std::size_t>(c)];
    }
    Bitboard occupied(Color c) const { return colors_[static_cast<std::size_t>(c)]; }
    Bitboard occupied() const { return colors_[0] | colors_[1]; }

    // What stands on a square; None on an empty one.
    PieceType piece_on(Square s) const {
        for (PieceType t : all_piece_types)
            if (types_[static_cast<std::size_t>(t)].contains(s))
                return t;
        return PieceType::None;
    }
    // Precondition: the square is occupied.
    Color color_on(Square s) const {
        return colors_[0].contains(s) ? Color::White : Color::Black;
    }

    // Precondition: the side has a king (every position from_setup accepts does).
    Square king_square(Color c) const { return pieces(c, PieceType::King).lowest(); }

    // Whether any piece of colour `by` attacks the square.
    bool attacked(Square s, Color by) const {
        const Bitboard occ = occupied();
        const Bitboard diagonal = pieces(by, PieceType::Bishop) | pieces(by, PieceType::Queen);
        const Bitboard straight = pieces(by, PieceType::Rook) | pieces(by, PieceType::Queen);
        // A pawn of colour `by` attacks s exactly when s is in the pattern of the opposite
        // colour from there.
        return !((pawn_attacks(opposite(by), s) & pieces(by, PieceType::Pawn)) |
                 (knight_attacks(s) & pieces(by, PieceType::Knight)) |
                 (king_attacks(s) & pieces(by, PieceType::King)) |
                 (bishop_attacks(s, occ) & diagonal) | (rook_attacks(s, occ) & straight))
                    .empty();
    }

    bool in_check() const { return attacked(king_square(side_), opposite(side_)); }

    // The position after a move. Precondition: the move is pseudo-legal here (what
    // `legal_moves` and its generator produce); it is not re-checked.
    Position made(const Move& m) const;

    // The same game seen from the other side of the table: ranks flipped, colours
    // swapped, the turn passed. Every rule is symmetric under it, castling included,
    // so a position and its flip have the same number of moves at every depth.
    Position color_flipped() const;

    // The left-right mirror. A symmetry of the rules only without castling rights:
    // the king starts on e, which has no mirror image.
    Position mirrored_files() const;

    bool operator==(const Position&) const = default;

  private:
    void place(Color c, PieceType t, Square s) {
        types_[static_cast<std::size_t>(t)] = types_[static_cast<std::size_t>(t)].with(s);
        colors_[static_cast<std::size_t>(c)] = colors_[static_cast<std::size_t>(c)].with(s);
    }
    void remove(Color c, PieceType t, Square s) {
        types_[static_cast<std::size_t>(t)] = types_[static_cast<std::size_t>(t)].without(s);
        colors_[static_cast<std::size_t>(c)] = colors_[static_cast<std::size_t>(c)].without(s);
    }

    std::array<Bitboard, 6> types_{};
    std::array<Bitboard, 2> colors_{};
    Color side_ = Color::White;
    std::uint8_t castling_ = 0;
    std::optional<Square> en_passant_;
    int halfmove_ = 0;
    int fullmove_ = 1;
};

// ── Moves ─────────────────────────────────────────────────────────

// A fixed buffer: no allocation per position. 512 is above the most moves any position of at
// most sixteen pieces a side can offer (a queen reaches 27 squares, a pawn four promotions on
// each of three targets).
struct MoveList {
    std::array<Move, 512> moves;
    int size = 0;

    void push(const Move& m) { moves[static_cast<std::size_t>(size++)] = m; }
    const Move* begin() const { return moves.data(); }
    const Move* end() const { return moves.data() + size; }
};

namespace detail {

inline void push_pawn_move(MoveList& out, Square from, Square to, bool promotes) {
    if (!promotes) {
        out.push({from, to, PieceType::None});
        return;
    }
    for (PieceType t :
         {PieceType::Queen, PieceType::Rook, PieceType::Bishop, PieceType::Knight})
        out.push({from, to, t});
}

// Every move the rules allow before the question of leaving one's own king attacked.
inline void generate_pseudo_legal(const Position& pos, MoveList& out) {
    const Color us = pos.side_to_move(), them = opposite(us);
    const Bitboard own = pos.occupied(us), enemy = pos.occupied(them), all = pos.occupied();

    // Knights, bishops, rooks, queens and the king: the reach, minus our own pieces.
    const auto add_targets = [&](Square from, Bitboard targets) {
        for (Square to : targets & ~own)
            out.push({from, to, PieceType::None});
    };
    for (Square from : pos.pieces(us, PieceType::Knight))
        add_targets(from, knight_attacks(from));
    for (Square from : pos.pieces(us, PieceType::Bishop))
        add_targets(from, bishop_attacks(from, all));
    for (Square from : pos.pieces(us, PieceType::Rook))
        add_targets(from, rook_attacks(from, all));
    for (Square from : pos.pieces(us, PieceType::Queen))
        add_targets(from, queen_attacks(from, all));
    for (Square from : pos.pieces(us, PieceType::King))
        add_targets(from, king_attacks(from));

    // Pawns: a push, a double push from the home rank, captures, en passant, promotions.
    const int up = us == Color::White ? 8 : -8;
    const int home_rank = us == Color::White ? 1 : 6;
    const int last_rank = us == Color::White ? 7 : 0;
    for (Square from : pos.pieces(us, PieceType::Pawn)) {
        const Square one(from.index + up);
        if (!all.contains(one)) {
            push_pawn_move(out, from, one, one.rank() == last_rank);
            const Square two(from.index + 2 * up);
            if (from.rank() == home_rank && !all.contains(two))
                out.push({from, two, PieceType::None});
        }
        for (Square to : pawn_attacks(us, from) & enemy)
            push_pawn_move(out, from, to, to.rank() == last_rank);
        if (pos.en_passant() && pawn_attacks(us, from).contains(*pos.en_passant()))
            out.push({from, *pos.en_passant(), PieceType::None});
    }

    // Castling: the rights, the king and rook at home, the squares between empty, and the
    // king not in check, not crossing an attacked square, not landing on one.
    const int base = us == Color::White ? 0 : 56;
    const Square king_home(base + 4);
    if (pos.pieces(us, PieceType::King).contains(king_home)) {
        const std::uint8_t kingside = us == Color::White ? white_kingside : black_kingside;
        const std::uint8_t queenside = us == Color::White ? white_queenside : black_queenside;
        const auto safe = [&](int a, int b, int c) {
            return !pos.attacked(Square(base + a), them) &&
                   !pos.attacked(Square(base + b), them) &&
                   !pos.attacked(Square(base + c), them);
        };
        if ((pos.castling_rights() & kingside) &&
            pos.pieces(us, PieceType::Rook).contains(Square(base + 7)) &&
            !all.contains(Square(base + 5)) && !all.contains(Square(base + 6)) && safe(4, 5, 6))
            out.push({king_home, Square(base + 6), PieceType::None});
        if ((pos.castling_rights() & queenside) &&
            pos.pieces(us, PieceType::Rook).contains(Square(base)) &&
            !all.contains(Square(base + 1)) && !all.contains(Square(base + 2)) &&
            !all.contains(Square(base + 3)) && safe(4, 3, 2))
            out.push({king_home, Square(base + 2), PieceType::None});
    }
}

} // namespace detail

// Every legal move of the side to move.
inline MoveList legal_moves(const Position& pos) {
    MoveList pseudo;
    detail::generate_pseudo_legal(pos, pseudo);
    MoveList legal;
    const Color us = pos.side_to_move();
    for (const Move& m : pseudo) {
        const Position after = pos.made(m);
        if (!after.attacked(after.king_square(us), opposite(us)))
            legal.push(m);
    }
    return legal;
}

// ── How a game ends ───────────────────────────────────────────────
//
// The ends the rules give a position by itself: mate and stalemate (no legal move, in check
// or not), the fifty-move rule (a hundred plies with no capture and no pawn move -- the point
// a player may claim the draw; this file does not tell a claim from the seventy-five-move
// rule's automatic one), and a dead position for want of material. Mate comes first: a mate
// on the hundredth ply is a mate.
//
// The dead positions recognised are the ones no sequence of legal moves can mate in: king
// against king, a king and one knight or one bishop against a king, and kings with bishops
// only, all of them on squares of one colour. Two knights against a king can still be mated
// with the loser's help, so it is not dead, however hopeless.
//
// Repetition needs the game, not the position: `outcome(game)` also ends a game whose last
// position stood twice before since the last capture or pawn move. Two positions are the same
// for this when the same pieces stand on the same squares, the same side is to move, with the
// same castling rights and the same en passant captures possible. An en passant square that
// no pawn can capture to does not count, though `Position` keeps it.

enum class Outcome : std::uint8_t {
    Ongoing,
    Checkmate, // the side to move is mated
    Stalemate,
    FiftyMoves,
    DeadPosition,
    Repetition,
};

// The en passant square, if an en passant capture to it is legal now.
inline std::optional<Square> capturable_en_passant(const Position& pos) {
    const auto ep = pos.en_passant();
    if (!ep)
        return std::nullopt;
    for (const Move& m : legal_moves(pos))
        if (m.to == *ep && pos.piece_on(m.from) == PieceType::Pawn)
            return ep;
    return std::nullopt;
}

// The same position in the sense of the repetition rule.
inline bool same_for_repetition(const Position& a, const Position& b) {
    if (a.side_to_move() != b.side_to_move() || a.castling_rights() != b.castling_rights())
        return false;
    for (Color c : {Color::White, Color::Black})
        for (PieceType t : all_piece_types)
            if (a.pieces(c, t) != b.pieces(c, t))
                return false;
    return capturable_en_passant(a) == capturable_en_passant(b);
}

inline Outcome outcome(const Position& pos) {
    if (legal_moves(pos).size == 0)
        return pos.in_check() ? Outcome::Checkmate : Outcome::Stalemate;
    if (pos.halfmove_clock() >= 100)
        return Outcome::FiftyMoves;
    const auto dead = [&] {
        const auto both = [&](PieceType t) {
            return pos.pieces(Color::White, t) | pos.pieces(Color::Black, t);
        };
        if (!(both(PieceType::Pawn) | both(PieceType::Rook) | both(PieceType::Queen)).empty())
            return false;
        const Bitboard knights = both(PieceType::Knight), bishops = both(PieceType::Bishop);
        if (knights.count() + bishops.count() <= 1)
            return true;
        if (!knights.empty())
            return false;
        constexpr Bitboard dark(0xaa55aa55aa55aa55ull); // a1, c1, ..., b2, d2, ...
        return (bishops & dark).empty() || (bishops & dark) == bishops;
    };
    return dead() ? Outcome::DeadPosition : Outcome::Ongoing;
}

// `game` is every position of the game in order, the current one last. A repetition is looked
// for only back to the last capture or pawn move, since no position before one can recur.
inline Outcome outcome(std::span<const Position> game) {
    const Position& now = game.back();
    if (const Outcome alone = outcome(now); alone != Outcome::Ongoing)
        return alone;
    const std::size_t reach = std::min<std::size_t>(game.size() - 1, static_cast<std::size_t>(now.halfmove_clock()));
    int seen = 1;
    for (std::size_t back = 2; back <= reach; back += 2) // the same side to move: even plies back
        if (same_for_repetition(game[game.size() - 1 - back], now) && ++seen == 3)
            return Outcome::Repetition;
    return Outcome::Ongoing;
}

// ── Counting the tree of legal moves ──────────────────────────────

// The number of positions `depth` moves ahead. The last level is counted, not played.
inline std::uint64_t perft(const Position& pos, int depth) {
    if (depth <= 0)
        return 1;
    const MoveList moves = legal_moves(pos);
    if (depth == 1)
        return static_cast<std::uint64_t>(moves.size);
    std::uint64_t nodes = 0;
    for (const Move& m : moves)
        nodes += perft(pos.made(m), depth - 1);
    return nodes;
}

struct DivideEntry {
    Move move;
    std::uint64_t nodes;
};

// perft split by the first move: the way to find which move a rules error hides behind.
inline std::vector<DivideEntry> perft_divide(const Position& pos, int depth) {
    std::vector<DivideEntry> out;
    for (const Move& m : legal_moves(pos))
        out.push_back({m, perft(pos.made(m), depth - 1)});
    return out;
}

// ── Position: definitions ─────────────────────────────────────────

inline Position Position::made(const Move& m) const {
    Position next = *this;
    const Color us = side_, them = opposite(us);
    const PieceType moving = piece_on(m.from);
    const PieceType captured = piece_on(m.to);
    const bool en_passant_capture = moving == PieceType::Pawn && captured == PieceType::None &&
                                    en_passant_ && m.to == *en_passant_ &&
                                    m.from.file() != m.to.file();

    if (captured != PieceType::None)
        next.remove(them, captured, m.to);
    if (en_passant_capture)
        next.remove(them, PieceType::Pawn, Square::at(m.to.file(), m.from.rank()));

    next.remove(us, moving, m.from);
    next.place(us, m.promotion == PieceType::None ? moving : m.promotion, m.to);

    // The rook goes with a king that moves two files.
    if (moving == PieceType::King &&
        (m.to.file() - m.from.file() == 2 || m.from.file() - m.to.file() == 2)) {
        const int base = m.from.rank() * board_size;
        const bool kingside = m.to.file() > m.from.file();
        next.remove(us, PieceType::Rook, Square(base + (kingside ? 7 : 0)));
        next.place(us, PieceType::Rook, Square(base + (kingside ? 5 : 3)));
    }

    // A right is lost when the king or a rook leaves its square, or a rook is captured on it.
    const auto lost = [](Square s) -> std::uint8_t {
        switch (s.index) {
        case 4:
            return white_kingside | white_queenside;
        case 0:
            return white_queenside;
        case 7:
            return white_kingside;
        case 60:
            return black_kingside | black_queenside;
        case 56:
            return black_queenside;
        case 63:
            return black_kingside;
        default:
            return 0;
        }
    };
    next.castling_ = static_cast<std::uint8_t>(castling_ & ~(lost(m.from) | lost(m.to)));

    const int pushed = m.to.rank() - m.from.rank();
    if (moving == PieceType::Pawn && (pushed == 2 || pushed == -2))
        next.en_passant_ = Square::at(m.from.file(), m.from.rank() + pushed / 2);
    else
        next.en_passant_.reset();

    next.halfmove_ =
        (moving == PieceType::Pawn || captured != PieceType::None || en_passant_capture)
            ? 0
            : halfmove_ + 1;
    if (us == Color::Black)
        ++next.fullmove_;
    next.side_ = them;
    return next;
}

inline Position Position::color_flipped() const {
    Position out;
    for (std::size_t t = 0; t < 6; ++t)
        out.types_[t] = types_[t].flipped_ranks();
    out.colors_[0] = colors_[1].flipped_ranks();
    out.colors_[1] = colors_[0].flipped_ranks();
    out.side_ = opposite(side_);
    out.castling_ = static_cast<std::uint8_t>(((castling_ & 3) << 2) | ((castling_ >> 2) & 3));
    if (en_passant_)
        out.en_passant_ = flip_ranks.apply(*en_passant_);
    out.halfmove_ = halfmove_;
    out.fullmove_ = fullmove_;
    return out;
}

inline Position Position::mirrored_files() const {
    Position out = *this;
    for (std::size_t t = 0; t < 6; ++t)
        out.types_[t] = types_[t].mirrored_files();
    for (std::size_t c = 0; c < 2; ++c)
        out.colors_[c] = colors_[c].mirrored_files();
    // Kingside and queenside swap with the files.
    out.castling_ = static_cast<std::uint8_t>(((castling_ & 5) << 1) | ((castling_ >> 1) & 5));
    if (en_passant_)
        out.en_passant_ = mirror_files.apply(*en_passant_);
    return out;
}

inline Result<Position> Position::from_setup(const PositionSetup& setup) {
    const auto fail = [](std::string message) -> Result<Position> {
        return std::unexpected(Error{ErrorCode::InvalidArgument, "position: " + std::move(message)});
    };
    Position p;
    for (int s = 0; s < square_count; ++s)
        if (const auto& o = setup.squares[static_cast<std::size_t>(s)])
            p.place(o->color, o->type, Square(s));
    p.side_ = setup.side;
    if (setup.castling > 15)
        return fail("the castling rights are not a set of the four");
    p.castling_ = setup.castling;
    if (setup.en_passant) {
        if (setup.en_passant->rank() != (p.side_ == Color::White ? 5 : 2))
            return fail("the en passant square is on the wrong rank for the side to move");
        p.en_passant_ = setup.en_passant;
    }
    if (setup.halfmove < 0)
        return fail("the halfmove clock is negative");
    if (setup.fullmove < 1)
        return fail("the move number is not positive");
    p.halfmove_ = setup.halfmove;
    p.fullmove_ = setup.fullmove;

    for (Color c : {Color::White, Color::Black}) {
        if (p.pieces(c, PieceType::King).count() != 1)
            return fail("each side needs exactly one king");
        if (p.occupied(c).count() > 16)
            return fail("a side has more than sixteen pieces");
    }
    if (!(p.types_[static_cast<std::size_t>(PieceType::Pawn)] & Bitboard(0xff000000000000ffull))
             .empty())
        return fail("a pawn stands on the first or last rank");
    if (p.attacked(p.king_square(opposite(p.side_)), p.side_))
        return fail("the side that is not to move is in check");
    return p;
}

inline Position Position::initial() {
    PositionSetup setup;
    constexpr std::array<PieceType, 8> back{PieceType::Rook,  PieceType::Knight, PieceType::Bishop,
                                            PieceType::Queen, PieceType::King,   PieceType::Bishop,
                                            PieceType::Knight, PieceType::Rook};
    for (int file = 0; file < 8; ++file) {
        const auto at = [&](int rank) { return static_cast<std::size_t>(Square::at(file, rank).index); };
        setup.squares[at(0)] = PositionSetup::Occupant{Color::White, back[static_cast<std::size_t>(file)]};
        setup.squares[at(1)] = PositionSetup::Occupant{Color::White, PieceType::Pawn};
        setup.squares[at(6)] = PositionSetup::Occupant{Color::Black, PieceType::Pawn};
        setup.squares[at(7)] = PositionSetup::Occupant{Color::Black, back[static_cast<std::size_t>(file)]};
    }
    setup.castling = white_kingside | white_queenside | black_kingside | black_queenside;
    return *from_setup(setup);
}

}  // namespace spatium::discrete
