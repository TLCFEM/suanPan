#include "TestSolver.h"
#include <Toolbox/LBFGS.hpp>
#include "CatchHeader.h"
#include "Domain/MetaMat/BandMat.hpp"
#include "Domain/MetaMat/IterativeSolver.hpp"
#include "Domain/MetaMat/SparseMatSuperLU.hpp"

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
        System A(mat(randu(N, N) + 10 * eye(N, N)));
        const vec b(N, fill::randu);
        vec x;

        SolverSetting<double> setting;
        auto preconditioner = Jacobi<double>(A.A);
        setting.preconditioner = &preconditioner;

        GMRES(&A, x, b, setting);

        REQUIRE(norm(solve(A.A, b) - x) <= 1E-12);
    }
}

TEST_CASE("BiCGSTAB Solver", "[Utility.Solver]") {
    constexpr auto N = 100;
    for(auto I = 0; I < N; ++I) {
        System A(mat(randu(N, N) + 10 * eye(N, N)));
        const vec b(N, fill::randu);
        vec x;

        SolverSetting<double> setting;
        auto preconditioner = Jacobi<double>(A.A);
        setting.preconditioner = &preconditioner;

        BiCGSTAB(&A, x, b, setting);

        REQUIRE(norm(solve(A.A, b) - x) <= 1E-12);
    }
}

TEST_CASE("Iterative Solver Sparse", "[Matrix.Solver]") {
    for(auto I = 0; I < 100; ++I) {
        const auto N = randi<uword>(distr_param(100, 200));
        auto A = SparseMatSuperLU<double>(N, N);
        REQUIRE(A.n_rows == N);
        REQUIRE(A.n_cols == N);

        sp_mat B = sprandu(N, N, .02) + speye(N, N) * 1E1;

        const vec C = randu<vec>(N);
        vec x;

        A.zeros();
        for(auto J = B.begin(); J != B.end(); ++J) A.at(J.row(), J.col()) = *J;

        SolverSetting<double> setting;
        auto preconditioner = Jacobi<double>(A);
        setting.preconditioner = &preconditioner;

        BiCGSTAB(&A, x, C, setting);

        REQUIRE(norm(spsolve(B, C) - x) <= 1E-12);

        x.reset();
        setting.max_iteration = 500;
        GMRES(&A, x, C, setting);

        REQUIRE(norm(spsolve(B, C) - x) <= 1E-12);
    }
}

TEST_CASE("Iterative Solver Dense", "[Matrix.Solver]") {
    for(auto I = 0; I < 100; ++I) {
        const auto N = randi<uword>(distr_param(100, 200));
        auto A = BandMat<double>(N, N, 3);
        REQUIRE(A.n_rows == N);
        REQUIRE(A.n_cols == N);

        mat B = randu(N, N) + eye(N, N) * 1E1;

        const vec C = randu<vec>(N);
        vec x;

        A.zeros();
        for(auto i = 0llu; i < N; ++i)
            for(auto j = 0llu; j < N; ++j)
                if(std::abs(static_cast<int>(i) - static_cast<int>(j)) <= 3) A.at(i, j) = B(i, j);
                else B(i, j) = 0.;

        SolverSetting<double> setting;
        setting.iterative_solver = IterativeSolver::BICGSTAB;

        A.set_solver_setting(setting);
        A.iterative_solve(x, C);

        REQUIRE(norm(solve(B, C) - x) <= 1E1 * setting.tolerance);

        setting.iterative_solver = IterativeSolver::GMRES;

        A.set_solver_setting(setting);
        A.iterative_solve(x, C);

        REQUIRE(norm(solve(B, C) - x) <= 1E1 * setting.tolerance);
    }
}

TEST_CASE("Iterative Solver Sparse Mat", "[Matrix.Solver]") {
    for(auto I = 0; I < 100; ++I) {
        const auto N = randi<uword>(distr_param(100, 200));
        auto A = SparseMatSuperLU<double>(N, N);
        REQUIRE(A.n_rows == N);
        REQUIRE(A.n_cols == N);

        sp_mat B = sprandu(N, N, .02) + speye(N, N) * 1E1;

        const mat C = randu<mat>(N, N);
        mat x;

        A.zeros();
        for(auto J = B.begin(); J != B.end(); ++J) A.at(J.row(), J.col()) = *J;

        SolverSetting<double> setting;
        setting.iterative_solver = IterativeSolver::BICGSTAB;
        setting.preconditioner_type = PreconditionerType::ILU;

        A.iterative_solve(x, C);

        REQUIRE(norm(spsolve(B, C) - x) <= 1E-12);
    }
}
