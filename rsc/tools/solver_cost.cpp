// Which of two candidate implementations is actually the cheaper one.
//
// RSC's task generators label ground truth by checking candidates
// cheapest-first against a reference-quality last candidate: the cheap one
// wins whenever it is accurate enough, otherwise the expensive one does.
// That is the right notion -- cheapest of the correct -- and the accuracy
// half is genuinely measured, with the two-regime splits probed before
// being written down.
//
// The cost half is not measured. Which candidate is cheapest is fixed by
// its index in the registry, and for the linear domain that index says
// Jacobi. This checks it.
//
// Build: part of the rsc tools target.

#include <spatium/algebra/linear_solve.hpp>
#include <chrono>
#include <print>
#include <random>
#include <vector>
using namespace spatium;
template<std::size_t N>
static void bench(double diag) {
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> off(0.5, 2.0);
    Matrix<double,N,N> A{}; Vec<double,N> b{};
    for (std::size_t i=0;i<N;++i){ for (std::size_t j=0;j<N;++j) A(i,j) = (i==j)? diag : off(rng); b[i]=off(rng); }
    constexpr int reps = 200000;
    auto t0=std::chrono::steady_clock::now();
    volatile double s1=0; for(int i=0;i<reps;++i){ auto x=solve_direct(A,b); if(x) s1+=(*x)[0]; }
    auto t1=std::chrono::steady_clock::now();
    volatile double s2=0; for(int i=0;i<reps;++i){ auto x=solve_jacobi(A,b,20); s2+=x[0]; }
    auto t2=std::chrono::steady_clock::now();
    const double d=std::chrono::duration<double,std::nano>(t1-t0).count()/reps;
    const double j=std::chrono::duration<double,std::nano>(t2-t1).count()/reps;
    std::println("  N={} diag={:.1f} | direct {:>8.1f} ns | jacobi(20) {:>8.1f} ns | jacobi is {:.2f}x direct",
                 N, diag, d, j, j/d);
}
int main() {
    std::println("Which of the two linear solvers is actually cheaper?");
    std::println("RSC labels the cheapest-first candidate as index 0, and that is jacobi.");
    std::println("");
    bench<4>(6.0); bench<4>(20.0); bench<8>(6.0); bench<16>(6.0);
}
