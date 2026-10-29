// A position of the game: FEN in and out, every kind of move, check, and the count of the
// tree of legal moves (perft) against the published numbers for the standard test positions.
//
// A rules error shows up in perft as a wrong count at some depth, and perft_divide says which
// first move it hides behind. The default run goes to the depths under a quarter of a million
// positions; the deeper counts (up to five million) are tagged [.deep] and run on request:
//     test_position "[.deep]"
//
// The check that does not come from a table is the colour symmetry: a position and the same
// game seen from the other side of the table have the same number of moves at every depth,
// castling and en passant included. The left-right mirror is the same kind of check, and it
// is only a symmetry without castling rights, which a test holds as well.

#include <spatium/discrete/position.hpp>
#include <spatium/io/chess.hpp>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <initializer_list>
#include <random>
#include <span>
#include <sstream>
#include <string>
#include <vector>

using namespace spatium;
using namespace spatium::discrete;

namespace {

struct Known {
    const char* name;
    const char* fen;
    std::vector<std::uint64_t> perft; // perft[d - 1]
};

// The standard test positions and their published counts.
const std::vector<Known>& known_positions() {
    static const std::vector<Known> positions{
        {"the starting position",
         "rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1",
         {20, 400, 8902, 197281, 4865609}},
        {"Kiwipete: castling, pins and promotions at once",
         "r3k2r/p1ppqpb1/bn2pnp1/3PN3/1p2P3/2N2Q1p/PPPBBPPP/R3K2R w KQkq - 0 1",
         {48, 2039, 97862, 4085603}},
        {"rooks and pawns in an endgame, with en passant and a pinned pawn",
         "8/2p5/3p4/KP5r/1R3p1k/8/4P1P1/8 w - - 0 1",
         {14, 191, 2812, 43238, 674624}},
        {"promotions with capture, and castling into danger",
         "r3k2r/Pppp1ppp/1b3nbN/nP6/BBP1P3/q4N2/Pp1P2PP/R2Q1RK1 w kq - 0 1",
         {6, 264, 9467, 422333}},
        {"the same, seen from the other side",
         "r2q1rk1/pP1p2pp/Q4n2/bbp1p3/Np6/1B3NBn/pPPP1PPP/R3K2R b KQ - 0 1",
         {6, 264, 9467, 422333}},
        {"a middlegame with a knight about to promote",
         "rnbq1k1r/pp1Pbppp/2p5/8/2B5/8/PPP1NnPP/RNBQK2R w KQ - 1 8",
         {44, 1486, 62379, 2103487}},
        {"a quiet middlegame",
         "r4rk1/1pp1qppp/p1np1n2/2b1p1B1/2B1P1b1/P1NP1N2/1PP1QPPP/R4RK1 w - - 0 10",
         {46, 2079, 89890, 3894594}},
    };
    return positions;
}

Position from(const char* fen) {
    auto p = io::from_fen(fen);
    REQUIRE(p.has_value());
    return *p;
}

std::vector<std::string> uci_of(const Position& p) {
    std::vector<std::string> out;
    for (const Move& m : legal_moves(p))
        out.push_back(io::to_uci(m));
    std::sort(out.begin(), out.end());
    return out;
}

bool has_move(const Position& p, const char* uci) {
    const auto moves = uci_of(p);
    return std::find(moves.begin(), moves.end(), uci) != moves.end();
}

Position after(const Position& p, const char* uci) {
    for (const Move& m : legal_moves(p))
        if (io::to_uci(m) == uci)
            return p.made(m);
    FAIL("no legal move " << uci << " in " << io::to_fen(p));
    return p;
}

// Whether two positions have the same number of positions ahead, at every depth up to `depth`.
bool same_perft(const Position& a, const Position& b, int depth) {
    for (int d = 1; d <= depth; ++d)
        if (perft(a, d) != perft(b, d))
            return false;
    return true;
}

} // namespace

TEST_CASE("a FEN is read and written back unchanged", "[position][fen]") {
    for (const Known& k : known_positions()) {
        INFO(k.name);
        CHECK(io::to_fen(from(k.fen)) == k.fen);
    }
    CHECK(
        io::to_fen(Position::initial()) == "rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1"
    );
    // The two clocks are optional.
    CHECK(io::to_fen(from("4k3/8/8/8/8/8/8/4K3 w - -")) == "4k3/8/8/8/8/8/8/4K3 w - - 0 1");
}

TEST_CASE(
    "a FEN that is not a position the rules can be played from is refused", "[position][fen]"
) {
    // What the text cannot say: io rejects it as a parse error.
    const std::vector<const char*> unreadable{
        "",                                                         // nothing
        "8/8/8/8/8/8/8 w - -",                                      // seven ranks
        "rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBN w KQkq - 0 1",  // a short rank
        "rnbqkbnr/pppppppp/9/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1", // a rank of nine
        "rnbqkbnr/ppxppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1", // an unknown piece
        "4k3/8/8/8/8/8/8/4K3 x - - 0 1",                            // neither side
        "4k3/8/8/8/8/8/8/4K3 w KQx - 0 1",                          // a bad castling field
        "4k3/8/8/8/8/8/8/4K3 w - - x 1",                            // a clock that is not a number
    };
    // What the text says but the rules do not allow: from_setup rejects it.
    const std::vector<const char*> unplayable{
        "8/8/8/8/8/8/8/8 w - - 0 1",                                // no kings
        "4k3/8/8/8/8/8/8/3KK3 w - - 0 1",                           // two white kings
        "P3k3/8/8/8/8/8/8/4K3 w - - 0 1",                           // a pawn on the last rank
        "4k3/8/8/8/8/8/8/4K3 w - e5 0 1",                           // en passant on the wrong rank
        "4k3/8/8/8/8/8/8/4K3 w - - -1 1",                           // a negative clock
        "4k3/8/8/8/8/8/8/4K3 w - - 0 0",                            // move number zero
        "4k3/8/8/4R3/8/8/8/4K3 w - - 0 1", // the side not to move is in check
    };
    for (const char* fen : unreadable) {
        INFO(fen);
        const auto p = io::from_fen(fen);
        REQUIRE_FALSE(p.has_value());
        CHECK(p.error().code == ErrorCode::ParseError);
    }
    for (const char* fen : unplayable) {
        INFO(fen);
        const auto p = io::from_fen(fen);
        REQUIRE_FALSE(p.has_value());
        CHECK(p.error().code == ErrorCode::InvalidArgument);
    }
}

TEST_CASE("the starting position built square by square is the one its FEN reads", "[position][fen]") {
    const auto read = io::from_fen("rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1");
    REQUIRE(read.has_value());
    CHECK(*read == Position::initial());
}

TEST_CASE("the starting position offers the twenty moves", "[position][moves]") {
    const std::vector<std::string> expected{"a2a3", "a2a4", "b1a3", "b1c3", "b2b3", "b2b4", "c2c3",
                                            "c2c4", "d2d3", "d2d4", "e2e3", "e2e4", "f2f3", "f2f4",
                                            "g1f3", "g1h3", "g2g3", "g2g4", "h2h3", "h2h4"};
    CHECK(uci_of(Position::initial()) == expected);
    CHECK_FALSE(Position::initial().in_check());
}

TEST_CASE("a game's opening moves rewrite the position as the notation does", "[position][moves]") {
    // 1. e4 c5 2. Nf3 -- the Sicilian, with the en passant square written after each double push.
    Position p = Position::initial();
    p = after(p, "e2e4");
    CHECK(io::to_fen(p) == "rnbqkbnr/pppppppp/8/8/4P3/8/PPPP1PPP/RNBQKBNR b KQkq e3 0 1");
    p = after(p, "c7c5");
    CHECK(io::to_fen(p) == "rnbqkbnr/pp1ppppp/8/2p5/4P3/8/PPPP1PPP/RNBQKBNR w KQkq c6 0 2");
    p = after(p, "g1f3");
    CHECK(io::to_fen(p) == "rnbqkbnr/pp1ppppp/8/2p5/4P3/5N2/PPPP1PPP/RNBQKB1R b KQkq - 1 2");
}

TEST_CASE(
    "castling needs the right, the room, and a king that is safe the whole way", "[position][moves]"
) {
    const Position both = from("r3k2r/8/8/8/8/8/8/R3K2R w KQkq - 0 1");
    CHECK(has_move(both, "e1g1"));
    CHECK(has_move(both, "e1c1"));
    // The rook goes with the king, and the rights of that side are spent.
    CHECK(io::to_fen(after(both, "e1g1")) == "r3k2r/8/8/8/8/8/8/R4RK1 b kq - 1 1");
    CHECK(io::to_fen(after(both, "e1c1")) == "r3k2r/8/8/8/8/8/8/2KR3R b kq - 1 1");

    // No right, no castling.
    CHECK_FALSE(has_move(from("r3k2r/8/8/8/8/8/8/R3K2R w Qkq - 0 1"), "e1g1"));
    // A piece in the way.
    CHECK_FALSE(has_move(from("r3k2r/8/8/8/8/8/8/R3KB1R w KQkq - 0 1"), "e1g1"));
    // The king crosses an attacked square (f1 is attacked by the rook on f8) -- queenside is still
    // fine.
    const Position crossing = from("4kr2/8/8/8/8/8/8/R3K2R w KQ - 0 1");
    CHECK_FALSE(has_move(crossing, "e1g1"));
    CHECK(has_move(crossing, "e1c1"));
    // The rook, not the king, may be attacked on the queenside: b1 is only passed by the rook.
    CHECK(has_move(from("1r2k3/8/8/8/8/8/8/R3K3 w Q - 0 1"), "e1c1"));
    // Not out of check.
    const Position checked = from("4r1k1/8/8/8/8/8/8/R3K2R w KQ - 0 1");
    CHECK(checked.in_check());
    CHECK_FALSE(has_move(checked, "e1g1"));
    CHECK_FALSE(has_move(checked, "e1c1"));
}

TEST_CASE(
    "a right is lost when the king or a rook moves, or a rook is taken", "[position][moves]"
) {
    const Position start = from("r3k2r/8/8/8/8/8/8/R3K2R w KQkq - 0 1");
    CHECK(
        after(start, "h1g1").castling_rights() ==
        (white_queenside | black_kingside | black_queenside)
    );
    CHECK(after(start, "e1d1").castling_rights() == (black_kingside | black_queenside));
    // The rook on h1 takes the rook on h8: white loses the kingside right by moving, black by being
    // taken.
    const Position taken = after(start, "h1h8");
    CHECK(taken.castling_rights() == (white_queenside | black_queenside));
    CHECK(io::to_fen(taken) == "r3k2R/8/8/8/8/8/8/R3K3 b Qq - 0 1");
}

TEST_CASE("en passant takes the pawn that is not on the target square", "[position][moves]") {
    const Position p = from("rnbqkbnr/ppp1p1pp/8/3pPp2/8/8/PPPP1PPP/RNBQKBNR w KQkq f6 0 3");
    REQUIRE(has_move(p, "e5f6"));
    CHECK(
        io::to_fen(after(p, "e5f6")) == "rnbqkbnr/ppp1p1pp/5P2/3p4/8/8/PPPP1PPP/RNBQKBNR b KQkq - 0 3"
    );

    // The capture that would uncover the king along the rank is not legal: both pawns leave the
    // fifth rank.
    const Position pinned = from("8/8/8/KPp4r/8/8/8/7k w - c6 0 1");
    CHECK_FALSE(has_move(pinned, "b5c6"));
    CHECK(has_move(pinned, "b5b6"));
}

TEST_CASE("a pawn on the last step becomes one of four pieces", "[position][moves]") {
    const Position p = from("8/P6k/8/8/8/8/8/K7 w - - 0 1");
    for (const char* uci : {"a7a8q", "a7a8r", "a7a8b", "a7a8n"})
        CHECK(has_move(p, uci));
    CHECK_FALSE(has_move(p, "a7a8"));
    CHECK(io::to_fen(after(p, "a7a8q")) == "Q7/7k/8/8/8/8/8/K7 b - - 0 1");
    CHECK(io::to_fen(after(p, "a7a8n")) == "N7/7k/8/8/8/8/8/K7 b - - 0 1");
}

TEST_CASE("check, mate and stalemate are told apart", "[position][moves]") {
    // The fool's mate: in check with no move.
    const Position mate = from("rnb1kbnr/pppp1ppp/8/4p3/6Pq/5P2/PPPPP2P/RNBQKBNR w KQkq - 1 3");
    CHECK(mate.in_check());
    CHECK(legal_moves(mate).size == 0);
    // Stalemate: no move, and not in check.
    const Position stalemate = from("7k/5Q2/6K1/8/8/8/8/8 b - - 0 1");
    CHECK_FALSE(stalemate.in_check());
    CHECK(legal_moves(stalemate).size == 0);
    // Check with a way out: the king may step, nothing else is needed.
    CHECK(legal_moves(from("4k3/8/8/8/8/8/4r3/4K3 w - - 0 1")).size > 0);
}

TEST_CASE(
    "the count of the tree of legal moves matches the published numbers", "[position][perft]"
) {
    for (const Known& k : known_positions()) {
        const Position p = from(k.fen);
        for (std::size_t i = 0; i < k.perft.size() && k.perft[i] <= 250000; ++i) {
            INFO(k.name << ", depth " << i + 1);
            CHECK(perft(p, static_cast<int>(i) + 1) == k.perft[i]);
        }
    }
}

TEST_CASE("the deeper counts match as well", "[.deep][position][perft]") {
    for (const Known& k : known_positions()) {
        const Position p = from(k.fen);
        for (std::size_t i = 0; i < k.perft.size(); ++i) {
            if (k.perft[i] <= 250000)
                continue;
            INFO(k.name << ", depth " << i + 1);
            CHECK(perft(p, static_cast<int>(i) + 1) == k.perft[i]);
        }
    }
}

TEST_CASE("perft split by the first move adds up to perft", "[position][perft]") {
    const Position p = from(known_positions()[1].fen);
    for (int depth = 1; depth <= 3; ++depth) {
        std::uint64_t sum = 0;
        for (const DivideEntry& e : perft_divide(p, depth))
            sum += e.nodes;
        CHECK(sum == perft(p, depth));
    }
    CHECK(perft_divide(p, 1).size() == 48);
}

TEST_CASE(
    "the same game seen from the other side has the same number of moves", "[position][symmetry]"
) {
    for (const Known& k : known_positions()) {
        INFO(k.name);
        const Position p = from(k.fen);
        const Position flipped = p.color_flipped();
        CHECK(flipped.color_flipped() == p); // flipping twice is the identity
        CHECK(same_perft(p, flipped, 3));
    }
    // The starting position seen from the other side is itself, with Black to move.
    const Position start = Position::initial().color_flipped();
    CHECK(io::to_fen(start) == "rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR b KQkq - 0 1");
}

TEST_CASE(
    "the left-right mirror is a symmetry when there are no castling rights", "[position][symmetry]"
) {
    // The endgame has none.
    const Position endgame = from(known_positions()[2].fen);
    REQUIRE(endgame.castling_rights() == 0);
    CHECK(endgame.mirrored_files().mirrored_files() == endgame);
    CHECK(same_perft(endgame, endgame.mirrored_files(), 4));
}

TEST_CASE(
    "the mirror is not a symmetry with castling rights, and the check says so",
    "[position][symmetry]"
) {
    // The king of the mirror image stands on d, where no castling starts: it has fewer moves.
    const Position kiwipete = from(known_positions()[1].fen);
    CHECK(perft(kiwipete, 1) == 48);
    CHECK(perft(kiwipete.mirrored_files(), 1) < 48);
    CHECK_FALSE(same_perft(kiwipete, kiwipete.mirrored_files(), 2));
}

TEST_CASE("squares have names and pieces have letters", "[position][fen]") {
    for (int s = 0; s < square_count; ++s) {
        const Square sq(s);
        REQUIRE(io::parse_square(io::square_name(sq)).has_value());
        CHECK(*io::parse_square(io::square_name(sq)) == sq);
    }
    CHECK(io::square_name(Square::at(0, 0)) == "a1");
    CHECK(io::square_name(Square::at(7, 7)) == "h8");
    CHECK(io::square_name(Square::at(4, 3)) == "e4");
    for (const char* bad : {"", "a", "a9", "i1", "A1", "e44"})
        CHECK_FALSE(io::parse_square(bad).has_value());

    const char expected[] = {'p', 'n', 'b', 'r', 'q', 'k'};
    for (std::size_t i = 0; i < all_piece_types.size(); ++i)
        CHECK(io::piece_letter(all_piece_types[i]) == expected[i]);
}

TEST_CASE("the position says what stands on a square", "[position]") {
    const Position p = Position::initial();
    CHECK(p.piece_on(*io::parse_square("e1")) == PieceType::King);
    CHECK(p.piece_on(*io::parse_square("d8")) == PieceType::Queen);
    CHECK(p.piece_on(*io::parse_square("b1")) == PieceType::Knight);
    CHECK(p.piece_on(*io::parse_square("c8")) == PieceType::Bishop);
    CHECK(p.piece_on(*io::parse_square("a1")) == PieceType::Rook);
    CHECK(p.piece_on(*io::parse_square("e2")) == PieceType::Pawn);
    CHECK(p.piece_on(*io::parse_square("e4")) == PieceType::None);
    CHECK(p.color_on(*io::parse_square("e1")) == Color::White);
    CHECK(p.color_on(*io::parse_square("e8")) == Color::Black);
    CHECK(p.pieces(Color::White, PieceType::Pawn).count() == 8);
    CHECK(p.pieces(Color::Black, PieceType::King) == Bitboard::of(*io::parse_square("e8")));
    CHECK(p.king_square(Color::White) == *io::parse_square("e1"));
    CHECK(p.occupied().count() == 32);
}

TEST_CASE("the most moves any position offers fit the move list", "[position][moves]") {
    // The known position with 218 legal moves, the most found in a reachable position.
    const Position p = from("R6R/3Q4/1Q4Q1/4Q3/2Q4Q/Q4Q2/pp1Q4/kBNN1KB1 w - - 0 1");
    const MoveList moves = legal_moves(p);
    CHECK(moves.size == 218);
    int walked = 0;
    for (const Move& m : moves) {
        (void)m;
        ++walked;
    }
    CHECK(walked == moves.size);
    CHECK(moves.size < static_cast<int>(moves.moves.size()));
}

TEST_CASE("a move in standard algebraic notation is read back as the move it was written from",
          "[position][san]") {
    // Random games through every kind of move: the text of each legal move names that move, and
    // no other, in its position.
    std::mt19937_64 rng(11);
    long checked = 0, captures = 0, castles = 0, promotions = 0, checks = 0, mates = 0, tells = 0;
    for (int game = 0; game < 300; ++game) {
        Position p = Position::initial();
        for (int ply = 0; ply < 160; ++ply) {
            const MoveList moves = legal_moves(p);
            if (!moves.size)
                break;
            for (const Move& m : moves) {
                const std::string san = io::to_san(p, m);
                INFO(san << " in " << io::to_fen(p));
                const auto back = io::parse_san(p, san);
                REQUIRE(back.has_value());
                CHECK(*back == m);
                ++checked;
                captures += san.find('x') != std::string::npos;
                castles += san[0] == 'O';
                promotions += san.find('=') != std::string::npos;
                checks += san.back() == '+';
                mates += san.back() == '#';
                // a text for a piece that is told from another by its file or rank
                tells += san.size() > 3 && std::isupper(static_cast<unsigned char>(san[0])) &&
                         san[0] != 'O' && std::islower(static_cast<unsigned char>(san[2])) &&
                         san.find('x') == std::string::npos;
            }
            p = p.made(moves.moves[rng() % static_cast<std::size_t>(moves.size)]);
        }
    }
    CHECK(checked > 100000);
    CHECK(captures > 5000);
    CHECK(castles > 500);
    CHECK(promotions > 100);
    CHECK(checks > 2000);
    CHECK(mates > 20);
    CHECK(tells > 1000);
}

TEST_CASE("a game written as SAN is played out to its mate", "[position][san]") {
    // A game from chess.com that ends 21...Qd1#, castling and a capture on the way.
    const char* game =
        "e4 e5 d3 Bb4+ Nd2 d5 c3 Bc5 Nb3 Bb6 a4 a5 c4 Nf6 Bg5 O-O c5 Ba7 Qc2 dxe4 dxe4 Nbd7 "
        "Nf3 Nxc5 Nxc5 h6 h4 Rb8 Nxe5 Qd6 Qc3 Bxc5 g3 b5 axb5 Rxb5 Rxa5 Rb6 f4 Rd8 Bh3 Qd1#";
    std::istringstream in(game);
    Position p = Position::initial();
    int plies = 0;
    for (std::string token; in >> token; ++plies) {
        INFO(token << " in " << io::to_fen(p));
        const auto m = io::parse_san(p, token);
        REQUIRE(m.has_value());
        CHECK(io::to_san(p, *m) == token);
        p = p.made(*m);
    }
    CHECK(plies == 42);
    CHECK(p.in_check());
    CHECK(legal_moves(p).size == 0);
    CHECK(!io::parse_san(Position::initial(), "e5").has_value());
    CHECK(!io::parse_san(Position::initial(), "Nf6").has_value());
    // two knights to one square: the file tells them, then the rank, then both
    const Position q = from("4k3/8/8/8/8/2N3N1/8/4K3 w - - 0 1");
    CHECK(io::to_san(q, {*io::parse_square("c3"), *io::parse_square("e4")}) == "Nce4");
    const Position r = from("4k3/8/8/8/8/N7/8/N3K3 w - - 0 1");
    CHECK(io::to_san(r, {*io::parse_square("a1"), *io::parse_square("c2")}) == "N1c2");
    const Position t = from("4k3/8/8/8/8/2N3N1/8/4K3 w - - 0 1");
    CHECK(io::to_san(t, {*io::parse_square("g3"), *io::parse_square("e4")}) == "Nge4");
}

TEST_CASE("how a position ends a game by itself", "[position][outcome]") {
    // Fool's mate: the side to move is mated.
    CHECK(outcome(from("rnb1kbnr/pppp1ppp/8/4p3/6Pq/5P2/PPPPP2P/RNBQKBNR w KQkq - 1 3")) ==
          Outcome::Checkmate);
    // The king on h8 has no move and is not in check.
    CHECK(outcome(from("7k/5Q2/6K1/8/8/8/8/8 b - - 0 1")) == Outcome::Stalemate);
    CHECK(outcome(Position::initial()) == Outcome::Ongoing);

    // The hundredth ply with no capture or pawn move, and the one before it.
    CHECK(outcome(from("4k3/8/8/8/8/8/8/R3K3 w - - 100 80")) == Outcome::FiftyMoves);
    CHECK(outcome(from("4k3/8/8/8/8/8/8/R3K3 w - - 99 80")) == Outcome::Ongoing);
    // A mate on the hundredth ply is a mate.
    CHECK(outcome(from("R3k3/8/4K3/8/8/8/8/8 b - - 100 80")) == Outcome::Checkmate);

    // Dead for want of material, and the near misses that are not.
    for (const char* dead : {
             "8/8/8/4k3/8/8/8/4K3 w - - 0 1",    // king against king
             "8/8/8/4k3/8/8/8/2B1K3 w - - 0 1",  // and a bishop
             "8/8/8/4k3/8/8/8/1N2K3 w - - 0 1",  // and a knight
             "5b2/8/8/4k3/8/8/8/2B1K3 w - - 0 1" // a bishop each, both on dark squares (c1, f8)
         }) {
        INFO(dead);
        CHECK(outcome(from(dead)) == Outcome::DeadPosition);
    }
    for (const char* alive : {
             "2b5/8/8/4k3/8/8/8/2B1K3 w - - 0 1", // bishops on squares of both colours (c1, c8)
             "8/8/8/4k3/8/8/8/1NN1K3 w - - 0 1",  // two knights: mate needs the loser's help
             "8/8/8/4k3/8/8/4P3/4K3 w - - 0 1",   // a pawn
             "8/8/8/4k3/8/8/8/R3K3 w - - 0 1",    // a rook
         }) {
        INFO(alive);
        CHECK(outcome(from(alive)) == Outcome::Ongoing);
    }
}

TEST_CASE("a position that stands for the third time ends the game", "[position][outcome]") {
    const auto play = [](std::vector<Position> game, std::initializer_list<const char*> moves) {
        for (const char* uci : moves)
            game.push_back(after(game.back(), uci));
        return game;
    };

    SECTION("the knights out and back twice") {
        auto game = play({Position::initial()}, {"g1f3", "g8f6", "f3g1", "f6g8"});
        CHECK(outcome(std::span<const Position>(game)) == Outcome::Ongoing); // twice
        game = play(game, {"g1f3", "g8f6", "f3g1", "f6g8"});
        CHECK(outcome(std::span<const Position>(game)) == Outcome::Repetition);
        // one ply short of it
        game.pop_back();
        CHECK(outcome(std::span<const Position>(game)) == Outcome::Ongoing);
    }

    SECTION("an en passant square no pawn can take on does not make the position new") {
        // After 1. e4 FEN writes e3, but no black pawn can take there: the position is the
        // one after 3. Ng1 and after 5. Ng1, so 5. Ng1 is the third time.
        auto game = play({Position::initial()},
                         {"e2e4", "g8f6", "g1f3", "f6g8", "f3g1", "g8f6", "g1f3", "f6g8", "f3g1"});
        const Position& first = game[1];
        const Position& third = game.back();
        CHECK(first.en_passant().has_value());
        CHECK_FALSE(capturable_en_passant(first).has_value());
        CHECK_FALSE(first == third); // equal as values they are not
        CHECK(same_for_repetition(first, third));
        CHECK(outcome(std::span<const Position>(game)) == Outcome::Repetition);

        // Compared as values, the clocks differ and no earlier position matches at all.
        // Compared by the first four fields of their FEN -- placement, side, castling, en
        // passant, as a repetition count over FEN strings would -- the first one is missed
        // and there are only two.
        const auto fen_position = [](const Position& p) {
            const std::string fen = io::to_fen(p);
            std::size_t cut = 0;
            for (int field = 0; field < 4; ++field)
                cut = fen.find(' ', cut) + 1;
            return fen.substr(0, cut - 1);
        };
        int by_value = 0, by_fen = 0;
        for (std::size_t i = game.size() % 2 == 0 ? 1 : 0; i < game.size(); i += 2) {
            by_value += game[i] == third;
            by_fen += fen_position(game[i]) == fen_position(third);
        }
        CHECK(by_value == 1);
        CHECK(by_fen == 2);
    }

    SECTION("an en passant capture that is possible makes the position new") {
        // Black's d-pawn on d4 can take e3 after 1. e4 only: here the same placement with
        // the capture available and without it are two positions.
        const Position can = from("4k3/8/8/8/3pP3/8/8/4K3 b - e3 0 1");
        const Position cannot = from("4k3/8/8/8/3pP3/8/8/4K3 b - - 0 1");
        CHECK(capturable_en_passant(can).has_value());
        CHECK_FALSE(same_for_repetition(can, cannot));
    }
}
