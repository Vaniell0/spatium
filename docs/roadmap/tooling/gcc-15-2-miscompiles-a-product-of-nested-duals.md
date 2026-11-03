---
date: 2026-10-08
title: GCC 15.2 at -O3 miscompiles a product of nested Duals, and the matrix product is unrolled
stage: 7
files: [include/spatium/algebra/matrix.hpp, tests/test_spd.cpp]
open:
  - the reproducer below is not filed with GCC; a person with a bugzilla account should file it, with the flags and the one-flag cure
  - only the matrix product is unrolled; other loops over a Dual of a Dual (a dot product, a trace, a matrix-vector product, the sums in `Vec`) are not known to be affected and are not shown to be safe, and a loop over Dual<Dual<T>> in user code is as exposed as ever
  - the guard is one case in the test binary, so it catches the next compiler or flag that breaks the product only in the build types the tests are run in; CI's nix build is Release, which is what caught this one
---

The previous entry's pull request failed CI on one check, in the `nix build` job and nowhere else: the derivative of the affine-invariant SPD distance through `Dual2` was right in a Debug build and wrong in a Release one (0.418 against a difference quotient of -0.635), for the two `Dual2` cells of the matrix that had just been closed. Bisected with flags: right at -O0 and -O2, wrong at -O3 and -O2 -ftree-vectorize with the default cost model, right again with `-fno-tree-loop-vectorize` and with the sanitizers; `-fno-tree-slp-vectorize`, `-fno-strict-aliasing`, `-ffp-contract=off` and `-fno-inline-functions` change nothing. Narrowed to the product of 2 by 2 matrices whose entries are `Dual<Dual<double>>`: a chain A V A of symmetric matrices came back with (0,1) and (1,0) different. Narrowed further with no library in it at all, to this, which prints a symmetric chain at -O2 and an asymmetric one at -O3 on GCC 15.2.0:

```cpp
#include <cstdio>
template<class T> struct Dl {
    T value{}; T deriv{};
    constexpr Dl() = default;
    constexpr Dl(T v, T d) : value(v), deriv(d) {}
    friend constexpr Dl operator+(const Dl& a, const Dl& b) { return {a.value + b.value, a.deriv + b.deriv}; }
    friend constexpr Dl operator*(const Dl& a, const Dl& b) { return {a.value * b.value, a.deriv * b.value + a.value * b.deriv}; }
};
using D1 = Dl<double>;
using D2 = Dl<D1>;
struct M { D2 m[4]; };                                   // column-major 2 by 2
__attribute__((noinline)) M mul(const M& a, const M& b) {
    M r;
    for (int c = 0; c < 2; ++c) for (int rw = 0; rw < 2; ++rw)
        for (int k = 0; k < 2; ++k) r.m[c*2+rw] = r.m[c*2+rw] + a.m[k*2+rw] * b.m[c*2+k];
    return r;
}
int main() {
    auto mk = [](double a, double b, double d) { M m; m.m[0]=D2{D1{a,0},D1{0,0}}; m.m[1]=m.m[2]=D2{D1{b,0},D1{0,0}}; m.m[3]=D2{D1{d,0},D1{0,0}}; return m; };
    M a = mk(0.79, -0.02, 1.08), v = mk(0.824062, -0.538622, 0.368402);
    D2 t{D1{0,0}, D1{1,0}};
    M vt; for (int i = 0; i < 4; ++i) vt.m[i] = v.m[i] * t;
    M chain = mul(mul(a, vt), a);
    for (int i = 0; i < 4; ++i) std::printf("chain[%d] %.6f\n", i, chain.m[i].deriv.value);   // [1] must equal [2]
}
```

`g++ -O2` prints 0.531465 -0.480745 -0.480745 0.453302; `g++ -O3` prints 0.509572 0.701488 -0.460609 -0.634085. Accumulating in a local, a `#pragma GCC novector` on the inner loop and every other source change tried left it wrong; writing the four sums out by hand fixed it. `Matrix::operator*` is therefore unrolled at compile time for the matrices of up to 64 products (all of this library's), and a test holds a chain of symmetric `Dual2` matrices symmetric in whatever build type the tests run. With it the cell is L3 at -O0, -O2 and -O3 -march=native, and the Release tree passes the SPD and matrix tests.
