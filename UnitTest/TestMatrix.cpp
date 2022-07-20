#include <Domain/MetaMat/MetaMat>
#include "CatchHeader.h"

template<typename MT, std::invocable T> void test_mat_solve(MT& A, const vec& D, const vec& C, T clear_mat) {
    constexpr double tol = 1E-12;

    vec E(C.n_elem);

    clear_mat();

    // full solve
    A.solve(E, C);
    REQUIRE(norm(E - D) < tol);

    // factored solve
    A.solve(E, C);
    REQUIRE(norm(E - D) < tol);

    clear_mat();

    // r-value full solve
    A.solve(E, mat(C));
    REQUIRE(norm(E - D) < tol);

    // r-value factored solve
    A.solve(E, mat(C));
    REQUIRE(norm(E - D) < tol);

    // mixed precision
    A.get_solver_setting().precision = Precision::MIXED;
    A.get_solver_setting().tolerance = 1E-18;

    clear_mat();

    // full solve
    A.solve(E, C);
    REQUIRE(norm(E - D) < tol);

    // factored solve
    A.solve(E, C);
    REQUIRE(norm(E - D) < tol);

    clear_mat();

    // r-value full solve
    A.solve(E, mat(C));
    REQUIRE(norm(E - D) < tol);

    // r-value factored solve
    A.solve(E, mat(C));
    REQUIRE(norm(E - D) < tol);
}

template<typename MT, std::invocable T> void benchmark_mat_solve(string&& title, MT& A, const vec& C, const vec& E, T&& clear_mat) {
    constexpr double tol = 1E-12;

    vec D;

    BENCHMARK((title + " Full").c_str()) {
        clear_mat();
        A.solve(D, C);
        REQUIRE(norm(E - D) < tol);
    };

    A.get_solver_setting().precision = Precision::MIXED;

    BENCHMARK((title + " Mixed").c_str()) {
        clear_mat();
        A.solve(D, C);
        REQUIRE(norm(E - D) < tol);
    };
}

template<typename T> void benchmark_mat_setup(const int I) {
    sp_mat B = sprandu(I, I, .01);
    B = B + B.t() + speye(I, I) * 1E1;
    const vec C = randu<vec>(I);

    auto A = T(I, I);

    string title;

    if(std::is_same_v<SparseMatSuperLU<double>, T>) title = "SuperLU ";
    else if(std::is_same_v<FullMat<double>, T>) title = "Full ";
#ifdef SUANPAN_CUDA
    else if(std::is_same_v<FullMatCUDA<double>, T>) title = "Full CUDA ";
    else if(std::is_same_v<SparseMatCUDA<double>, T>) title = "Sparse CUDA ";
#endif

    benchmark_mat_solve(std::move(title += std::to_string(I)), A, C, spsolve(B, C), [&] {
        A.zeros();
        for(auto J = B.begin(); J != B.end(); ++J) A.at(J.row(), J.col()) = *J;
    });
}

TEST_CASE("Mixed Precision", "[Matrix.Benchmark]") {
    for(auto I = 32; I < 1024; I *= 2) {
        benchmark_mat_setup<SparseMatSuperLU<double>>(I);
        benchmark_mat_setup<FullMat<double>>(I);
#ifdef SUANPAN_CUDA
        benchmark_mat_setup<FullMatCUDA<double>>(I);
        benchmark_mat_setup<SparseMatCUDA<double>>(I);
#endif
    }
}

TEST_CASE("FullMat", "[Matrix.Dense]") {
    for(auto I = 0; I < 100; ++I) {
        const auto N = randi<uword>(distr_param(100, 200));
        auto A = FullMat<double>(N, N);
        REQUIRE(A.n_rows == N);
        REQUIRE(A.n_cols == N);

        const mat B = randu<mat>(N, N) + eye(N, N) * 1E1;

        auto clear_mat = [&] {
            A.zeros();
            for(auto i = 0llu; i < N; ++i) for(auto j = 0llu; j < N; ++j) A.at(i, j) = B(i, j);
        };

        const vec C = randu<vec>(N);

        clear_mat();

        REQUIRE(norm(A * C - B * C) < 1E-12);
        REQUIRE(norm(A * B - B * B) < 1E-12);

        test_mat_solve(A, solve(B, C), C, clear_mat);
    }
}

TEST_CASE("SymmPackMat", "[Matrix.Dense]") {
    for(auto I = 0; I < 100; ++I) {
        const auto N = randi<uword>(distr_param(100, 200));
        auto A = SymmPackMat<double>(N);
        REQUIRE(A.n_rows == N);

        mat B = diagmat(randu<vec>(N)) + eye(N, N) * 1E1;

        auto clear_mat = [&] {
            A.zeros();
            for(auto i = 0llu; i < N; ++i) for(auto j = 0llu; j < N; ++j) A.at(i, j) = B(i, j);
        };

        const vec C = randu<vec>(N);

        clear_mat();

        REQUIRE(norm(A * C - B * C) < 1E-12);
        REQUIRE(norm(A * B - B * B) < 1E-12);

        test_mat_solve(A, solve(B, C), C, clear_mat);
    }
}

TEST_CASE("BandMat", "[Matrix.Dense]") {
    for(auto I = 0; I < 100; ++I) {
        const auto N = randi<uword>(distr_param(100, 200));
        auto A = BandMat<double>(N, N, 3);
        REQUIRE(A.n_rows == N);
        REQUIRE(A.n_cols == N);

        mat B = randu(N, N) + eye(N, N) * 1E1;

        auto clear_mat = [&] {
            A.zeros();
            for(auto i = 0llu; i < N; ++i)
                for(auto j = 0llu; j < N; ++j)
                    if(std::abs(static_cast<int>(i) - static_cast<int>(j)) <= 3) A.at(i, j) = B(i, j);
                    else B(i, j) = 0.;
        };

        const vec C = randu<vec>(N);

        clear_mat();

        REQUIRE(norm(A * C - B * C) < 1E-12);
        REQUIRE(norm(A * B - B * B) < 1E-12);

        test_mat_solve(A, solve(B, C), C, clear_mat);
    }
}

TEST_CASE("BandMatSpike", "[Matrix.Dense]") {
    for(auto I = 0; I < 100; ++I) {
        const auto N = randi<uword>(distr_param(100, 200));
        auto A = BandMatSpike<double>(N, N, 3);
        REQUIRE(A.n_rows == N);
        REQUIRE(A.n_cols == N);

        mat B = randu(N, N) + eye(N, N) * 1E1;

        auto clear_mat = [&] {
            A.zeros();
            for(auto i = 0llu; i < N; ++i)
                for(auto j = 0llu; j < N; ++j)
                    if(std::abs(static_cast<int>(i) - static_cast<int>(j)) <= 3) A.at(i, j) = B(i, j);
                    else B(i, j) = 0.;
        };

        const vec C = randu<vec>(N);

        clear_mat();

        REQUIRE(norm(A * C - B * C) < 1E-12);
        REQUIRE(norm(A * B - B * B) < 1E-12);

        test_mat_solve(A, solve(B, C), C, clear_mat);
    }
}

TEST_CASE("BandSymmMat", "[Matrix.Dense]") {
    for(auto I = 0; I < 100; ++I) {
        const auto N = randi<uword>(distr_param(100, 200));
        auto A = BandSymmMat<double>(N, 3);
        REQUIRE(A.n_rows == N);
        REQUIRE(A.n_cols == N);

        mat B = randu(N, N);
        B = B + B.t() + eye(N, N) * 1E1;

        auto clear_mat = [&] {
            A.zeros();
            for(auto i = 0llu; i < N; ++i)
                for(auto j = 0llu; j < N; ++j)
                    if(std::abs(static_cast<int>(i) - static_cast<int>(j)) <= 3) A.at(i, j) = B(i, j);
                    else B(i, j) = 0.;
        };

        const vec C = randu<vec>(N);

        clear_mat();

        REQUIRE(norm(A * C - B * C) < 1E-12);
        REQUIRE(norm(A * B - B * B) < 1E-12);

        test_mat_solve(A, solve(B, C), C, clear_mat);
    }
}

TEST_CASE("SparseMatSuperLU", "[Matrix.Sparse]") {
    for(auto I = 0; I < 100; ++I) {
        const auto N = randi<uword>(distr_param(100, 200));
        auto A = SparseMatSuperLU<double>(N, N);
        REQUIRE(A.n_rows == N);
        REQUIRE(A.n_cols == N);

        sp_mat B = sprandu(N, N, .01) + speye(N, N) * 1E1;

        auto clear_mat = [&] {
            A.zeros();
            for(auto J = B.begin(); J != B.end(); ++J) A.at(J.row(), J.col()) = *J;
        };

        const vec C = randu<vec>(N);

        clear_mat();

        REQUIRE(norm(A * C - B * C) < 1E-12);

        test_mat_solve(A, spsolve(B, C), C, clear_mat);
    }
}

TEST_CASE("SparseMatMUMPS", "[Matrix.Sparse]") {
    for(auto I = 0; I < 100; ++I) {
        const auto N = randi<uword>(distr_param(100, 200));
        auto A = SparseMatMUMPS<double>(N, N);
        REQUIRE(A.n_rows == N);
        REQUIRE(A.n_cols == N);

        sp_mat B = sprandu(N, N, .01) + speye(N, N) * 1E1;

        auto clear_mat = [&] {
            A.zeros();
            for(auto J = B.begin(); J != B.end(); ++J) A.at(J.row(), J.col()) = *J;
        };

        const vec C = randu<vec>(N);

        clear_mat();

        REQUIRE(norm(A * C - B * C) < 1E-12);

        test_mat_solve(A, spsolve(B, C), C, clear_mat);
    }
}

#ifdef SUANPAN_MKL
TEST_CASE("SparseMatPARDISO", "[Matrix.Sparse]") {
    for(auto I = 0; I < 100; ++I) {
        const auto N = randi<uword>(distr_param(100, 200));
        auto A = SparseMatPARDISO<double>(N, N);
        REQUIRE(A.n_rows == N);
        REQUIRE(A.n_cols == N);

        sp_mat B = sprandu(N, N, .01) + speye(N, N) * 1E1;

        auto clear_mat = [&] {
            A.zeros();
            for(auto J = B.begin(); J != B.end(); ++J) A.at(J.row(), J.col()) = *J;
        };

        const vec C = randu<vec>(N);

        clear_mat();

        REQUIRE(norm(A * C - B * C) < 1E-12);

        test_mat_solve(A, spsolve(B, C), C, clear_mat);
    }
}
#endif

#ifdef SUANPAN_CUDA
TEST_CASE("SparseMatCUDA", "[Matrix.Sparse]") {
    for(auto I = 0; I < 100; ++I) {
        const auto N = randi<uword>(distr_param(100, 200));
        auto A = SparseMatCUDA<double>(N, N);
        REQUIRE(A.n_rows == N);
        REQUIRE(A.n_cols == N);

        sp_mat B = sprandu(N, N, .01) + speye(N, N) * 1E1;

        auto clear_mat = [&] {
            A.zeros();
            for(auto J = B.begin(); J != B.end(); ++J) A.at(J.row(), J.col()) = *J;
        };

        const vec C = randu<vec>(N);

        clear_mat();

        REQUIRE(norm(A * C - B * C) < 1E-12);

        test_mat_solve(A, spsolve(B, C), C, clear_mat);
    }
}
#endif

TEST_CASE("Triplet/CSR/CSC Sparse", "[Matrix.Sparse]") {
    constexpr auto N = 100;

    triplet_form<double, uword> B(N, N);

    const vec C(N, fill::randn);

    for(auto I = 0; I < 200; ++I) {
        const mat A(sprandu<sp_mat>(N, N, .01));

        B.zeros();

        B.assemble(A, linspace<uvec>(0, N - 1llu, N));

        csr_form<double, uword> D(B);
        csc_form<double, uword> E(B);

        vec F = A * C;

        REQUIRE(norm(F - B * C) <= 1E-13);
        REQUIRE(norm(F - D * C) <= 1E-13);
        REQUIRE(norm(F - E * C) <= 1E-13);

        F = trimatu(A, 1) * C;

        REQUIRE(norm(F - B.strictly_upper() * C) <= 1E-13);

        F = trimatl(A, -1) * C;

        REQUIRE(norm(F - B.strictly_lower() * C) <= 1E-13);
    }
}

TEST_CASE("Benchmark Triplet Assembly", "[Matrix.Sparse]") {
    constexpr unsigned long long N = 1024;
    constexpr unsigned long long REPEAT = 8;
    constexpr unsigned long long NNZ = 1024;

    const triplet_form<double, uword> B(sprandu<sp_mat>(N, N, NNZ * pow(static_cast<double>(N), -2.)));

    REQUIRE(B.n_elem == NNZ);

    for(auto J = 2; J != REPEAT; J *= 2)
        BENCHMARK(string("Assemble " + std::to_string(J)).c_str()) {
            triplet_form<double, uword> C(N + REPEAT, N + REPEAT, B.n_elem * REPEAT);

            for(auto I = 0; I < J; ++I) C.assemble(B, I, I, randu<double>());

            REQUIRE(C.n_elem == NNZ * J);

            C.csc_condense();

            return C;
        };
}

TEST_CASE("Triplet/CSR/CSC Conversion", "[Matrix.Sparse]") {
    constexpr auto N = 128;

    for(auto J = 2; J != N; J *= 2) {
        auto A = mat(sprandu<sp_mat>(N, N, .5));
        const auto INDEX = linspace<uvec>(0, N - 1, N);

        triplet_form<double, uword> B(N, N);

        B.assemble(A, INDEX);
        B.zeros();
        B.assemble(A, INDEX);
        B.assemble(A, INDEX);

        csr_form<double, uword> C(B);
        csc_form<double, uword> D(B);
        csr_form<double, uword> CC(B, SparseBase::ZERO, true);
        csc_form<double, uword> DC(B, SparseBase::ZERO, true);

        const auto E = to_mat(B);
        const auto F = to_mat(C);
        const auto G = to_mat(D);
        const auto FF = to_mat(CC);
        const auto GG = to_mat(DC);

        A *= 2.;

        REQUIRE(norm(A - E) <= 1E-13);
        REQUIRE(norm(A - F) <= 1E-13);
        REQUIRE(norm(A - G) <= 1E-13);
        REQUIRE(norm(A - FF) <= 1E-13);
        REQUIRE(norm(A - GG) <= 1E-13);
    }
}

TEST_CASE("Benchmark Triplet Measure", "[Matrix.Sparse]") {
    constexpr unsigned long long N = 1024;
    constexpr unsigned long long REPEAT = 8;
    constexpr unsigned long long NNZ = 1024;

    const triplet_form<double, uword> B(sprandu<sp_mat>(N, N, NNZ * pow(static_cast<double>(N), -2.)));

    REQUIRE(B.n_elem == NNZ);

    for(auto J = 2; J != REPEAT; J *= 2) {
        constexpr unsigned long long S = 100;
        std::chrono::duration<double> assemble_mean(0);
        std::chrono::duration<double> compress_mean(0);
        for(auto K = 0llu; K < S; ++K) {
            triplet_form<double, uword> C(N + REPEAT, N + REPEAT, B.n_elem * REPEAT);

            auto start = std::chrono::high_resolution_clock::now();
            for(auto I = 0; I < J; ++I) C.assemble(B, I, I, randu<double>());
            auto end = std::chrono::high_resolution_clock::now();
            assemble_mean += end - start;

            REQUIRE(C.n_elem == NNZ * J);

            start = std::chrono::high_resolution_clock::now();
            C.csc_condense();
            end = std::chrono::high_resolution_clock::now();
            compress_mean += end - start;
        }
        // suanpan_info("Assemble: %.3f\n", assemble_mean.count() / static_cast<double>(S));
        // suanpan_info("Compress: %.3f\n", compress_mean.count() / static_cast<double>(S));
    }
}

template<typename T> void test_dense_mat_unify(T A) {
    constexpr auto N = 4;

    constexpr auto V = 2.31212;

    A.at(N, N) = V;
    REQUIRE(Approx(A(N, N)) == V);

    A.unify(N);
    REQUIRE(Approx(A(N, N)) == 1.);

    A.nullify(N);
    REQUIRE(Approx(A(N, N)) == 0.);
}

template<typename T> void test_sparse_mat_unify(T A) {
    constexpr auto N = 4;

    constexpr auto V = 2.31212;

    A.at(N, N) = V;
    REQUIRE(Approx(A(N, N)) == V);

    A.unify(N);
    A.csc_condense();
    REQUIRE(Approx(A(N, N)) == 1.);

    A.nullify(N);
    A.csc_condense();
    REQUIRE(Approx(A(N, N)) == 0.);

    A.unify(N);
    A.csc_condense();
    REQUIRE(Approx(A(N, N)) == 1.);

    A.unify(N);
    A.csr_condense();
    REQUIRE(Approx(A(N, N)) == 1.);

    A.nullify(N);
    A.csr_condense();
    REQUIRE(Approx(A(N, N)) == 0.);

    A.unify(N);
    A.unify(N);
    A.csr_condense();
    REQUIRE(Approx(A(N, N)) == 1.);
}

TEST_CASE("Unify FullMat", "[Matrix.Utility]") { test_dense_mat_unify(FullMat<double>(10, 10)); }

TEST_CASE("Unify BandMat", "[Matrix.Utility]") { test_dense_mat_unify(BandMat<double>(10, 2, 3)); }

TEST_CASE("Unify BandSymmMat", "[Matrix.Utility]") { test_dense_mat_unify(BandSymmMat<double>(10, 2)); }

TEST_CASE("Unify BandMatSpike", "[Matrix.Utility]") { test_dense_mat_unify(BandMatSpike<double>(10, 2, 3)); }

TEST_CASE("Unify SymmPackMat", "[Matrix.Utility]") { test_dense_mat_unify(SymmPackMat<double>(10)); }

TEST_CASE("Unify SparseMatSuperLU", "[Matrix.Utility]") { test_sparse_mat_unify(SparseMatSuperLU<double>(10, 10)); }

TEST_CASE("Unify SparseMatMUMPS", "[Matrix.Utility]") { test_sparse_mat_unify(SparseMatMUMPS<double>(10, 10)); }

#ifdef SUANPAN_MKL
TEST_CASE("Unify SparseMatPARDISO", "[Matrix.Utility]") { test_sparse_mat_unify(SparseMatPARDISO<double>(10, 10)); }
#endif

#ifdef SUANPAN_CUDA
TEST_CASE("Unify SparseMatCUDA", "[Matrix.Utility]") { test_sparse_mat_unify(SparseMatCUDA<double>(10, 10)); }
#endif
