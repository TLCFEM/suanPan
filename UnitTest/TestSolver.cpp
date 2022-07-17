#include <Toolbox/LBFGS.hpp>
#include <Toolbox/GMRES.hpp>
#include <Toolbox/BiCGSTAB.hpp>
#include "CatchHeader.h"

TEST_CASE("LBFGS Solver", "[Utility.Solver]") {
    Quadratic function;

    vec result{2., 3.};

    LBFGS().optimize(function, result);

    for(const auto I : result)
        REQUIRE(Approx(1.) == I);
}

TEST_CASE("GMRES Solver", "[Utility.Solver]") {
    constexpr auto N = 100;
    for(auto I = 0; I < N; ++I) {
        constexpr int m = 40;
        const mat A = randu(N, N) + 10 * eye(N, N);
        const vec b(N, fill::randu);
        vec x;

        DiagonalPreconditioner preconditioner(A);

        int max_iteration = 200;
        double tolerance = 1E-10;
        GMRES(A, x, b, preconditioner, m, max_iteration, tolerance);

        REQUIRE(norm(solve(A, b) - x) <= 1E2 * N * tolerance);
    }
}

TEST_CASE("BiCGSTAB Solver", "[Utility.Solver]") {
    constexpr auto N = 100;
    for(auto I = 0; I < N; ++I) {
        const mat A = randu(N, N) + 10 * eye(N, N);
        const vec b(N, fill::randu);
        vec x;

        DiagonalPreconditioner preconditioner(A);

        int max_iteration = 500;
        double tolerance = 1E-10;
        BiCGSTAB(A, x, b, preconditioner, max_iteration, tolerance);

        REQUIRE(norm(solve(A, b) - x) <= 1E3 * N * tolerance);
    }
}
