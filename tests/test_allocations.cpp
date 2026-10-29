// The hot paths of the library allocate nothing on the heap.
//
// The global operator new is replaced in this translation unit, which replaces it for the
// whole test binary; it only counts, then calls malloc. A count is exact and does not depend
// on the machine's load, so it can be held where a time cannot. Each path is run once first
// (a table built on first use is not a per-call cost) and then counted over many calls.
//
// What allocates on purpose is held to its number: an orbit and the cells a lattice piece
// reaches come back as vectors, one allocation each, reserved to their bound. The counter
// gets a case it must see: an allocation stored where the optimiser cannot drop it.

#include <spatium/algebra/calculus.hpp>
#include <spatium/algebra/dual.hpp>
#include <spatium/algebra/groups/action.hpp>
#include <spatium/algebra/groups/hyperoctahedral.hpp>
#include <spatium/algebra/linear_solve.hpp>
#include <spatium/algebra/matrix.hpp>
#include <spatium/algebra/quadrature.hpp>
#include <spatium/discrete/lattice.hpp>
#include <spatium/discrete/position.hpp>
#include <spatium/io/chess.hpp>
#include <spatium/io/field.hpp>

#include <atomic>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdlib>
#include <new>

namespace {
std::atomic<long> heap_allocations{0};
}

void* operator new(std::size_t n) {
    heap_allocations.fetch_add(1, std::memory_order_relaxed);
    if (void* p = std::malloc(n ? n : 1))
        return p;
    throw std::bad_alloc();
}
void* operator new[](std::size_t n) { return operator new(n); }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }

namespace {

using namespace spatium;

volatile double sink = 0;

// Heap allocations per call of f, over `calls` calls after one warm-up call.
template <class F> double allocations_per_call(F&& f, int calls = 200) {
    f();
    const long before = heap_allocations.load();
    for (int i = 0; i < calls; ++i)
        f();
    return double(heap_allocations.load() - before) / calls;
}

void* volatile kept = nullptr;

} // namespace

TEST_CASE("the counter sees an allocation it cannot be optimised out of", "[allocations]") {
    const double per = allocations_per_call([] {
        void* p = operator new(16);
        operator delete(kept);
        kept = p;
    });
    CHECK(per == 1.0);
    operator delete(kept);
    kept = nullptr;
}

TEST_CASE("algebra's hot paths allocate nothing", "[allocations][algebra]") {
    CHECK(allocations_per_call([] {
              const auto g = gradient([](const auto& x) { return x[0] * x[1] + sin(x[2]) * x[0]; },
                                      Vec3{1, 2, 3});
              sink = sink + g[0];
          }) == 0.0);
    CHECK(allocations_per_call([] {
              Matrix<double, 4, 4> a = Matrix<double, 4, 4>::identity();
              a(0, 1) = 2;
              sink = sink + (a * a)(0, 1);
          }) == 0.0);
    CHECK(allocations_per_call([] {
              Matrix<double, 8, 8> a = Matrix<double, 8, 8>::identity();
              Vec<double, 8> b{};
              b[0] = 1;
              sink = sink + (*solve_direct(a, b))[0];
          }) == 0.0);
    CHECK(allocations_per_call([] {
              sink = sink + integrate([](double x) { return std::sin(x); }, Finite<double>{0, 3}).value;
          }, 20) == 0.0);
    CHECK(allocations_per_call([] {
              sink = sink + integrate([](double x) { return std::sin(x); }, Finite<double>{0, 3},
                                      {.witness = true})
                                .value;
          }, 20) == 0.0);
}

TEST_CASE("a field of up to 32 ops is evaluated without the heap", "[allocations][field]") {
    using F = io::build::Field<double>;
    const F small = sin(F::u()) * F::v() + F::u();
    REQUIRE(small.size() <= 32);
    CHECK(allocations_per_call([&] { sink = sink + small(0.3, 0.4); }) == 0.0);
}

TEST_CASE("the rules of chess allocate nothing", "[allocations][discrete]") {
    const auto kiwipete =
        *io::from_fen("r3k2r/p1ppqpb1/bn2pnp1/3PN3/1p2P3/2N2Q1p/PPPBBPPP/R3K2R w KQkq - 0 1");
    CHECK(allocations_per_call([&] { sink = sink + discrete::legal_moves(kiwipete).size; }) == 0.0);
    CHECK(allocations_per_call([&] { sink = sink + double(discrete::perft(kiwipete, 2)); }, 5) == 0.0);
    CHECK(allocations_per_call([&] { sink = sink + double(discrete::outcome(kiwipete)); }) == 0.0);
    CHECK(allocations_per_call([&] {
              sink = sink + io::from_fen("4k3/8/8/8/8/8/8/4K3 w - - 0 1")->halfmove_clock();
          }) == 0.0);
    const auto first = discrete::legal_moves(kiwipete).moves[0];
    CHECK(allocations_per_call([&] { sink = sink + double(io::to_san(kiwipete, first).size()); }) == 0.0);
}

TEST_CASE("an orbit and a piece's reach are one reserved vector each", "[allocations][discrete]") {
    // at most one: an optimiser may drop an allocation that does not escape
    CHECK(allocations_per_call([] {
              sink = sink + double(orbit(HyperoctahedralGroup<3>{}, SignedPermutationAction<3>{},
                                         IntVec<3>{1, 2, 0})
                                       .size());
          }) <= 1.0);
    const discrete::Shape<2> board{{8, 8}};
    for (const auto& piece : {discrete::knight<2>(), discrete::queen<2>(), discrete::king<2>()})
        CHECK(allocations_per_call([&] {
                  sink = sink + double(discrete::attacks(piece, board, {3, 3},
                                                         [](const discrete::Cell<2>&) { return false; })
                                           .size());
              }) <= 1.0);
}
