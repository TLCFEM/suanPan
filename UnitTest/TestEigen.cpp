#include <Domain/MetaMat/FullMat.hpp>
#include <Toolbox/arpack.h>
#include <Toolbox/utility.h>
#include "CatchHeader.h"

TEST_CASE("Eigensolver", "[Utility.Eigen]") {
    constexpr auto N = 100;
    constexpr auto Q = 6;

    const vec D = regspace(1, 1, N);

    for(auto L = 0; L < N; ++L) {
        const mat P = orth(randn(D.n_elem, D.n_elem));

        mat K = P * diagmat(D) * P.t();

        mat M = 2. * eye(size(K));

        auto KK = make_shared<FullMat<double>>(D.n_elem, D.n_elem);

        for(auto I = 0llu; I < D.n_elem; ++I) for(auto J = 0llu; J < D.n_elem; ++J) KK->at(J, I) = K(J, I);

        auto MM = make_shared<FullMat<double>>(D.n_elem, D.n_elem);

        for(auto I = 0llu; I < D.n_elem; ++I) MM->at(I, I) = M(I, I);

        vec eigval;
        mat eigvec;

        REQUIRE(eig_solve(eigval, eigvec, KK, MM->make_copy(), Q, "SM") == 0);

        for(auto I = 0; I < Q; ++I)
            REQUIRE(Approx(eigval(I)) == .5 * I + .5);

        REQUIRE(eig_solve(eigval, eigvec, KK, MM->make_copy(), Q, "LM") == 0);

        for(auto I = 0; I < Q; ++I)
            REQUIRE(Approx(eigval(Q - 1 - I)) == .5 * (N - I));

        cx_vec cx_eigval;
        cx_mat cx_eigvec;

        REQUIRE(eig_solve(cx_eigval, cx_eigvec, KK->make_copy(), MM->make_copy(), Q, "LM") == 0);

        for(auto I = 0; I < Q; ++I)
            REQUIRE(Approx(cx_eigval(I).real()) == .5 * I + .5);
    }
}
