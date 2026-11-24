// Small standalone check of Toolbox/fgmres.hpp against the real oneMKL RCI
// routines.  Build with build.sh; no part of suanPan itself is compiled.
#include <Toolbox/fgmres.hpp>

#include <cmath>
#include <cstdio>

namespace {

#define CHECK(condition, message)                                                                                                                          \
    do {                                                                                                                                                   \
        if(!(condition)) {                                                                                                                                 \
            std::fprintf(stderr, "FAIL: %s\n", message);                                                                                                   \
            return 1;                                                                                                                                      \
        }                                                                                                                                                  \
    } while(false)

struct DenseKernel {
    const mat* matrix;
    uword n_rows;

    mat operator*(const mat& x) const { return *matrix * x; }
};

} // namespace

int main() {
    constexpr uword n = 8;
    mat a(n, n, arma::fill::zeros);
    for(uword i = 0; i < n; ++i) {
        a(i, i) = 10. + static_cast<double>(i);
        if(i + 1 < n) {
            a(i, i + 1) = -1.;
            a(i + 1, i) = -1.;
        }
        if(i + 2 < n) {
            a(i, i + 2) = .5;
            a(i + 2, i) = .5;
        }
    }

    vec exact(n);
    for(uword i = 0; i < n; ++i) exact[i] = static_cast<double>(i + 1);
    const vec b = a * exact;
    const vec diagonal = a.diag();
    const DenseKernel kernel{&a, n};

    vec x(n, arma::fill::zeros);
    MKL_INT iterations = -1;
    const auto code = fgmres_solve(kernel, diagonal, x.memptr(), const_cast<double*>(b.memptr()), 1E-12, &iterations);
    std::printf("converged solve: code=%d iterations=%d\n", code, iterations);
    CHECK(0 == code, "convergence expected");
    CHECK(iterations > 0, "iteration count not reported");
    const double relative_residual = arma::norm(a * x - b, 2) / arma::norm(b, 2);
    std::printf("relative residual: %.3e\n", relative_residual);
    CHECK(std::isfinite(relative_residual) && relative_residual < 1E-10, "residual too large");

    // Invalid input is rejected before any oneMKL call.
    CHECK(-1 == fgmres_solve(kernel, diagonal, nullptr, const_cast<double*>(b.memptr()), 1E-12), "null solution accepted");
    CHECK(-1 == fgmres_solve(kernel, vec(3), x.memptr(), const_cast<double*>(b.memptr()), 1E-12), "wrong preconditioner size accepted");
    CHECK(-1 == fgmres_solve(kernel, diagonal, x.memptr(), const_cast<double*>(b.memptr()), 0.), "zero tolerance accepted");
    CHECK(-1 == fgmres_solve(kernel, diagonal, x.memptr(), const_cast<double*>(b.memptr()), std::nan("")), "NaN tolerance accepted");

    // An unattainable tolerance must be reported, not silently called a success.
    x.zeros();
    iterations = -1;
    const auto impossible = fgmres_solve(kernel, diagonal, x.memptr(), const_cast<double*>(b.memptr()), 1E-30, &iterations);
    std::printf("unattainable tolerance: code=%d iterations=%d\n", impossible, iterations);
    CHECK(-7 == impossible, "non-convergence not reported");
    CHECK(iterations > 0, "iteration count missing on non-convergence");

    std::puts("all FGMRES tests passed");
    return 0;
}
