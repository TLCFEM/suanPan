#include "CatchHeader.h"

#include <Domain/Factory.hpp>
#include <Domain/MetaMat/MetaMat>
#include <Toolbox/eigen_solver.h>

namespace {
    template<typename MT> void assign(MT& in_mat, const mat& in_ref) {
        for(uword I{0}; I < in_ref.n_rows; ++I)
            for(uword J{0}; J < in_ref.n_cols; ++J) in_mat.at(I, J) = in_ref(I, J);
    }

    mat tridiagonal(const uword N) {
        mat K(N, N, fill::zeros);

        for(uword I{0}; I < N; ++I) {
            K(I, I) = 2.;
            if(I > 0) K(I, I - 1) = -1.;
            if(I < N - 1) K(I, I + 1) = -1.;
        }

        return K;
    }

    void require_eigen(const vec& eigval, const vec& ref, const double tol = 1E-8) {
        REQUIRE(eigval.n_elem == ref.n_elem);

        const vec sorted = sort(eigval);

        for(uword I{0}; I < ref.n_elem; ++I) REQUIRE(sorted(I) == Approx(ref(I)).epsilon(tol));
    }

    vec tridiagonal_eigen(const uword N) {
        vec eigval;
        mat eigvec;

        eig_sym(eigval, eigvec, tridiagonal(N));

        return sort(eigval);
    }
} // namespace

#ifndef SUANPAN_DISTRIBUTED
TEST_CASE("FEASTDense", "[Utility.FEAST]") {
    constexpr auto N = 8;

    mat K(N, N, fill::zeros), M(N, N, fill::zeros);

    for(uword I{0}; I < N; ++I) {
        K(I, I) = static_cast<double>(I + 1);
        M(I, I) = 1.;
    }

    auto KK = std::make_shared<FullMat<double>>(N, N);
    auto MM = std::make_shared<FullMat<double>>(N, N);

    assign(*KK, K);
    assign(*MM, M);

    vec eigval;
    mat eigvec;

    REQUIRE(feast_solve(eigval, eigvec, KK, MM, StorageScheme::FULL, N, .5 * (N + 1), 1. * N) == SUANPAN_SUCCESS);
    REQUIRE(eigvec.n_cols == eigval.n_elem);

    require_eigen(eigval, regspace(1., 1., static_cast<double>(N)));
}

TEST_CASE("FEASTBand", "[Utility.FEAST]") {
    constexpr auto N = 10;

    const mat K = tridiagonal(N);
    const mat M = eye<mat>(N, N);

    auto KK = std::make_shared<BandMat<double>>(N, 1, 1);
    auto MM = std::make_shared<BandMat<double>>(N, 1, 1);

    assign(*KK, K);
    assign(*MM, M);

    vec eigval;
    mat eigvec;

    REQUIRE(feast_solve(eigval, eigvec, KK, MM, StorageScheme::BAND, N, 2., 1.95, 1, 1) == SUANPAN_SUCCESS);
    REQUIRE(eigvec.n_cols == eigval.n_elem);

    require_eigen(eigval, tridiagonal_eigen(N));
}

TEST_CASE("FEASTBandSymm", "[Utility.FEAST]") {
    constexpr auto N = 10;

    const mat K = tridiagonal(N);
    const mat M = eye<mat>(N, N);

    auto KK = std::make_shared<BandSymmMat<double>>(N, 1);
    auto MM = std::make_shared<BandSymmMat<double>>(N, 1);

    assign(*KK, K);
    assign(*MM, M);

    vec eigval;
    mat eigvec;

    REQUIRE(feast_solve(eigval, eigvec, KK, MM, StorageScheme::BANDSYMM, N, 2., 1.95, 1, 1) == SUANPAN_SUCCESS);
    REQUIRE(eigvec.n_cols == eigval.n_elem);

    require_eigen(eigval, tridiagonal_eigen(N));
}
#endif

TEST_CASE("FEASTSparse", "[Utility.FEAST]") {
    constexpr auto N = 10;

    const mat K = tridiagonal(N);
    const mat M = eye<mat>(N, N);

    auto KK = std::make_shared<SparseMatSuperLU<double>>(N, N);
    auto MM = std::make_shared<SparseMatSuperLU<double>>(N, N);

    assign(*KK, K);
    assign(*MM, M);

    vec eigval;
    mat eigvec;

    REQUIRE(feast_solve(eigval, eigvec, KK, MM, StorageScheme::SPARSE, N, 2., 1.95) == SUANPAN_SUCCESS);
    REQUIRE(eigvec.n_cols == eigval.n_elem);

    require_eigen(eigval, tridiagonal_eigen(N));
}
