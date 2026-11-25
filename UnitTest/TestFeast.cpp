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

    shared_ptr<Factory<double>> make_factory(const shared_ptr<MetaMat<double>>& K, const shared_ptr<MetaMat<double>>& M, const StorageScheme scheme, const unsigned low_band = 0, const unsigned up_band = 0) {
        auto factory = std::make_shared<Factory<double>>(K->n_rows, AnalysisType::EIGEN, scheme);
        factory->set_stiffness(K);
        factory->set_mass(M);
        factory->set_bandwidth(low_band, up_band);
        return factory;
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

    const auto factory = make_factory(KK, MM, StorageScheme::FULL);

    REQUIRE(feast_solve(factory, N, .5 * (N + 1), 1. * N) == SUANPAN_SUCCESS);
    REQUIRE(factory->get_eigenvector().n_cols == factory->get_eigenvalue().n_elem);

    require_eigen(factory->get_eigenvalue(), regspace(1., 1., static_cast<double>(N)));
}

TEST_CASE("FEASTBand", "[Utility.FEAST]") {
    constexpr auto N = 10;

    const mat K = tridiagonal(N);
    const mat M = eye<mat>(N, N);

    auto KK = std::make_shared<BandMat<double>>(N, 1, 1);
    auto MM = std::make_shared<BandMat<double>>(N, 1, 1);

    assign(*KK, K);
    assign(*MM, M);

    const auto factory = make_factory(KK, MM, StorageScheme::BAND, 1, 1);

    REQUIRE(feast_solve(factory, N, 2., 1.95) == SUANPAN_SUCCESS);
    REQUIRE(factory->get_eigenvector().n_cols == factory->get_eigenvalue().n_elem);

    require_eigen(factory->get_eigenvalue(), tridiagonal_eigen(N));
}

TEST_CASE("FEASTSpike", "[Utility.FEAST]") {
    constexpr auto N = 10;

    const mat K = tridiagonal(N);
    const mat M = eye<mat>(N, N);

    auto KK = std::make_shared<BandMatSpike<double>>(N, 1, 1);
    auto MM = std::make_shared<BandMatSpike<double>>(N, 1, 1);

    assign(*KK, K);
    assign(*MM, M);

    const auto factory = make_factory(KK, MM, StorageScheme::BAND, 1, 1);
    factory->set_solver_type(SolverType::SPIKE);

    REQUIRE(feast_solve(factory, N, 2., 1.95) == SUANPAN_SUCCESS);
    REQUIRE(factory->get_eigenvector().n_cols == factory->get_eigenvalue().n_elem);

    require_eigen(factory->get_eigenvalue(), tridiagonal_eigen(N));
}

TEST_CASE("FEASTBandSymm", "[Utility.FEAST]") {
    constexpr auto N = 10;

    const mat K = tridiagonal(N);
    const mat M = eye<mat>(N, N);

    auto KK = std::make_shared<BandSymmMat<double>>(N, 1);
    auto MM = std::make_shared<BandSymmMat<double>>(N, 1);

    assign(*KK, K);
    assign(*MM, M);

    const auto factory = make_factory(KK, MM, StorageScheme::BANDSYMM, 1, 1);

    REQUIRE(feast_solve(factory, N, 2., 1.95) == SUANPAN_SUCCESS);
    REQUIRE(factory->get_eigenvector().n_cols == factory->get_eigenvalue().n_elem);

    require_eigen(factory->get_eigenvalue(), tridiagonal_eigen(N));
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

    const auto factory = make_factory(KK, MM, StorageScheme::SPARSE);

    REQUIRE(feast_solve(factory, N, 2., 1.95) == SUANPAN_SUCCESS);
    REQUIRE(factory->get_eigenvector().n_cols == factory->get_eigenvalue().n_elem);

    require_eigen(factory->get_eigenvalue(), tridiagonal_eigen(N));
}
