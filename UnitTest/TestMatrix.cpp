#include <Domain/MetaMat/MetaMat>
#include "CatchHeader.h"

template<typename MT, typename ET, std::invocable T> void test_mat_solve(MT& A, const Mat<ET>& D, const Col<ET>& C, T clear_mat) {
    constexpr auto tol = std::numeric_limits<ET>::epsilon() * 100;
    const auto scaled_tol = static_cast<ET>(C.n_elem) * tol;

    Col<ET> E(C.n_elem);

    clear_mat();

    // full solve
    A.solve(E, C);
    REQUIRE(arma::norm<Col<ET>>(E - D) < scaled_tol);

    // factored solve
    A.solve(E, C);
    REQUIRE(arma::norm<Col<ET>>(E - D) < scaled_tol);

    clear_mat();

    // r-value full solve
    A.solve(E, Mat<ET>(C));
    REQUIRE(arma::norm<Col<ET>>(E - D) < scaled_tol);

    // r-value factored solve
    A.solve(E, Mat<ET>(C));
    REQUIRE(arma::norm<Col<ET>>(E - D) < scaled_tol);

    // mixed precision
    A.get_solver_setting().precision = Precision::MIXED;
    A.get_solver_setting().tolerance = tol;

    clear_mat();

    // full solve
    A.solve(E, C);
    REQUIRE(arma::norm<Col<ET>>(E - D) < scaled_tol);

    // factored solve
    A.solve(E, C);
    REQUIRE(arma::norm<Col<ET>>(E - D) < scaled_tol);

    clear_mat();

    // r-value full solve
    A.solve(E, Mat<ET>(C));
    REQUIRE(arma::norm<Col<ET>>(E - D) < scaled_tol);

    // r-value factored solve
    A.solve(E, Mat<ET>(C));
    REQUIRE(arma::norm<Col<ET>>(E - D) < scaled_tol);
}

template<typename ET, std::invocable<u64> F> void test_dense_mat_setup(F new_mat) {
    constexpr auto tol = std::numeric_limits<ET>::epsilon() * 1000;
    for(auto I = 0; I < 100; ++I) {
        const auto N = randi<uword>(distr_param(100, 200));
        auto A = new_mat(N);
        REQUIRE(A.n_rows == N);
        REQUIRE(A.n_cols == N);

        auto B = randu<Mat<ET>>(N, N);
        B = B + B.t() + eye<decltype(B)>(N, N) * 10;

        auto clear_mat = [&] {
            A.zeros();
            for(auto i = 0llu; i < N; ++i)
                for(auto j = 0llu; j < N; ++j)
                    if(std::abs(static_cast<int>(i) - static_cast<int>(j)) <= 3) A.at(i, j) = B(i, j);
                    else B(i, j) = ET(0);
        };

        const auto C = randu<Col<ET>>(N);

        clear_mat();

        REQUIRE(arma::norm<Col<ET>>(A * C - B * C) < tol);
        REQUIRE(arma::norm<Mat<ET>>(A * B - B * B) < tol);

        test_mat_solve(A, solve(B, C).eval(), C, clear_mat);
    }
}

template<typename ET, std::invocable<u64> F> void test_sparse_mat_setup(F new_mat) {
    constexpr auto tol = std::numeric_limits<ET>::epsilon() * 1000;
    for(auto I = 0; I < 100; ++I) {
        const auto N = randi<uword>(distr_param(100, 200));
        auto A = new_mat(N);
        REQUIRE(A.n_rows == N);
        REQUIRE(A.n_cols == N);

        SpMat<ET> B = sprandu<SpMat<ET>>(N, N, .01) + speye<SpMat<ET>>(N, N) * 10;

        auto clear_mat = [&] {
            A.zeros();
            for(auto J = B.begin(); J != B.end(); ++J) A.at(J.row(), J.col()) = *J;
        };

        const auto C = randu<Col<ET>>(N);

        clear_mat();

        REQUIRE(arma::norm<Col<ET>>(A * C - B * C) < tol);

        test_mat_solve(A, spsolve(B, C), C, clear_mat);
    }
}

template<typename MT, typename ET, std::invocable T> void benchmark_mat_solve(string&& title, MT& A, const Col<ET>& C, const Mat<ET>& E, T&& clear_mat) {
    constexpr auto tol = std::numeric_limits<ET>::epsilon() * 1000;
    const auto scaled_tol = static_cast<ET>(C.n_elem) * tol;

    Col<ET> D;

    BENCHMARK((title + " Full").c_str()) {
        clear_mat();
        A.solve(D, C);
        REQUIRE(norm(E - D) < scaled_tol);
    };

    A.get_solver_setting().precision = Precision::MIXED;

    BENCHMARK((title + " Mixed").c_str()) {
        clear_mat();
        A.solve(D, C);
        REQUIRE(norm(E - D) < scaled_tol);
    };
}

template<typename T> T create_new(u64) { throw std::runtime_error("unknown matrix"); }

template<> FullMat<double> create_new(const u64 N) { return {N, N}; }

template<> SymmPackMat<double> create_new(const u64 N) { return SymmPackMat<double>{N}; }

template<> BandMat<double> create_new(const u64 N) { return {N, 3, 3}; }

template<> BandMatSpike<double> create_new(const u64 N) { return {N, N / 10, N / 10}; }

template<> BandSymmMat<double> create_new(const u64 N) { return {N, 3}; }

template<> SparseMatSuperLU<double> create_new(const u64 N) { return {N, N}; }

template<> SparseMatMUMPS<double> create_new(const u64 N) { return {N, N}; }

template<> SparseSymmMatMUMPS<double> create_new(const u64 N) { return {N, N}; }

template<> FullMat<float> create_new(const u64 N) { return {N, N}; }

template<> SymmPackMat<float> create_new(const u64 N) { return SymmPackMat<float>{N}; }

template<> BandMat<float> create_new(const u64 N) { return {N, 3, 3}; }

template<> BandMatSpike<float> create_new(const u64 N) { return {N, N / 10, N / 10}; }

template<> BandSymmMat<float> create_new(const u64 N) { return {N, 3}; }

template<> SparseMatSuperLU<float> create_new(const u64 N) { return {N, N}; }

template<> SparseMatMUMPS<float> create_new(const u64 N) { return {N, N}; }

template<> SparseMatLis<double> create_new(const u64 N) { return {N, N}; }

#ifdef SUANPAN_MKL
template<> SparseMatPARDISO<double> create_new(const u64 N) { return {N, N}; }

template<> SparseMatPARDISO<float> create_new(const u64 N) { return {N, N}; }

#ifdef SUANPAN_MPI
template<> SparseMatMPIPARDISO<double> create_new(const u64 N) { return {N, N}; }

template<> SparseMatMPIPARDISO<float> create_new(const u64 N) { return {N, N}; }
#endif

template<> SparseMatFGMRES<double> create_new(const u64 N) { return {N, N}; }

template<> SparseSymmMatFGMRES<double> create_new(const u64 N) { return {N, N}; }
#endif

#ifdef SUANPAN_CUDA
template<> BandMatCUDA<double> create_new(const u64 N) { return {N, 3, 3}; }

template<> FullMatCUDA<double> create_new(const u64 N) { return {N, N}; }

template<> FullMatCUDA<float> create_new(const u64 N) { return {N, N}; }

template<> SparseMatCUDA<double> create_new(const u64 N) { return {N, N}; }

template<> SparseMatCUDA<float> create_new(const u64 N) { return {N, N}; }

#ifdef SUANPAN_MAGMA
template<> SparseMatMAGMA<double> create_new(const u64 N) {
    istringstream dummy{"--verbose 1"};
    return {N, N, magma_parse_opts<magma_dopts>(dummy)};
}

template<> SparseMatMAGMA<float> create_new(const u64 N) {
    istringstream dummy{"--verbose 1"};
    return {N, N, magma_parse_opts<magma_sopts>(dummy)};
}
#endif
#endif

template<typename T, typename ET> void benchmark_mat_setup(const int I) {
    auto B = sprandu<SpMat<ET>>(I, I, .2);
    B = B + B.t() + speye<SpMat<ET>>(I, I) * 1E1;

    std::vector<std::tuple<uword, uword>> to_erase;
    to_erase.reserve(B.n_nonzero);
    for(auto J = B.begin(); J != B.end(); ++J) if(std::abs(static_cast<int>(J.row()) - static_cast<int>(J.col())) > 3) to_erase.emplace_back(std::tuple{J.row(), J.col()});
    for(const auto& [row, col] : to_erase) B.at(row, col) = 0.;

    const auto C = randu<Col<ET>>(I);

    auto A = create_new<T>(I);

    string title;

    if(std::is_same_v<FullMat<ET>, T>) title = "Full ";
    else if(std::is_same_v<SymmPackMat<ET>, T>) title = "SymmPack ";
    else if(std::is_same_v<BandMat<ET>, T>) title = "Band ";
    else if(std::is_same_v<BandMatSpike<ET>, T>) title = "BandSpike ";
    else if(std::is_same_v<BandSymmMat<ET>, T>) title = "BandSymm ";
    else if(std::is_same_v<SparseMatMUMPS<ET>, T>) title = "MUMPS ";
    else if(std::is_same_v<SparseSymmMatMUMPS<ET>, T>) title = "MUMPS Symm ";
    else if(std::is_same_v<SparseMatLis<ET>, T>) title = "Lis ";
    else if(std::is_same_v<SparseMatSuperLU<ET>, T>) title = "SuperLU ";
#ifdef SUANPAN_MKL
    else if(std::is_same_v<SparseMatPARDISO<ET>, T>) title = "PARDISO ";
#ifdef SUANPAN_MPI
    else if(std::is_same_v<SparseMatMPIPARDISO<ET>, T>) title = "MPI PARDISO ";
#endif
    else if(std::is_same_v<SparseMatFGMRES<ET>, T>) title = "FGMRES ";
    else if(std::is_same_v<SparseSymmMatFGMRES<ET>, T>) title = "FGMRES Symm ";
#endif
#ifdef SUANPAN_CUDA
    else if(std::is_same_v<FullMatCUDA<ET>, T>) title = "Full CUDA ";
    else if(std::is_same_v<SparseMatCUDA<ET>, T>) title = "Sparse CUDA ";
    else if(std::is_same_v<BandMatCUDA<ET>, T>) title = "Band CUDA ";
#ifdef SUANPAN_MAGMA
    else if(std::is_same_v<SparseMatMAGMA<ET>, T>) title = "Magma ";
#endif
#endif

    title += "N=" + std::to_string(I) + " NZ=" + std::to_string(B.n_nonzero) + " NE=" + std::to_string(A.n_elem);

    benchmark_mat_solve(std::move(title), A, C, spsolve(B, C).eval(), [&] {
        A.zeros();
        for(auto J = B.begin(); J != B.end(); ++J) A.at(J.col(), J.row()) = *J;
    });
}

TEST_CASE("Mixed Precision", "[Matrix.Benchmark]") {
    for(auto I = 0x0020; I < 0x0100; I *= 2) {
        benchmark_mat_setup<FullMat<double>, double>(I);
        benchmark_mat_setup<SymmPackMat<double>, double>(I);
        benchmark_mat_setup<BandMat<double>, double>(I);
        benchmark_mat_setup<BandMatSpike<double>, double>(I);
        benchmark_mat_setup<BandSymmMat<double>, double>(I);
        benchmark_mat_setup<SparseMatMUMPS<double>, double>(I);
        benchmark_mat_setup<SparseSymmMatMUMPS<double>, double>(I);
        benchmark_mat_setup<SparseMatSuperLU<double>, double>(I);
#ifdef SUANPAN_MKL
        benchmark_mat_setup<SparseMatPARDISO<double>, double>(I);
        benchmark_mat_setup<SparseMatFGMRES<double>, double>(I);
        benchmark_mat_setup<SparseSymmMatFGMRES<double>, double>(I);
#endif
#ifdef SUANPAN_CUDA
        benchmark_mat_setup<SparseMatCUDA<double>, double>(I);
        benchmark_mat_setup<BandMatCUDA<double>, double>(I);
#endif
    }
#ifdef SUANPAN_CUDA
    for(auto I = 0x0100; I < 0x2000; I *= 2) benchmark_mat_setup<FullMatCUDA<double>, double>(I);
#endif
}

#ifdef SUANPAN_CUDA
TEST_CASE("Large BandMatCUDA", "[Matrix.Benchmark]") {
    benchmark_mat_setup<BandMatCUDA<double>, double>(0x1000);
    benchmark_mat_setup<BandMat<double>, double>(0x1000);
}
#endif

TEST_CASE("Large Mixed Precision", "[Matrix.Benchmark]") {
    for(auto I = 0x400; I < 0x500; I *= 2) {
        benchmark_mat_setup<BandMat<double>, double>(I);
        benchmark_mat_setup<BandMatSpike<double>, double>(I);
        benchmark_mat_setup<BandSymmMat<double>, double>(I);
        benchmark_mat_setup<SparseMatMUMPS<double>, double>(I);
        benchmark_mat_setup<SparseSymmMatMUMPS<double>, double>(I);
        benchmark_mat_setup<SparseMatLis<double>, double>(I);
        benchmark_mat_setup<SparseMatSuperLU<double>, double>(I);
    }
}

TEST_CASE("Large Sparse Solve Type", "[Matrix.Benchmark]") {
    for(auto I = 0x1000; I < 0x5000; I *= 2) {
        benchmark_mat_setup<BandMat<double>, double>(I);
        benchmark_mat_setup<SparseMatMUMPS<double>, double>(I);
        benchmark_mat_setup<SparseMatLis<double>, double>(I);
        benchmark_mat_setup<SparseMatSuperLU<double>, double>(I);
    }
}

TEST_CASE("FullMat", "[Matrix.Dense]") { test_dense_mat_setup<double>(create_new<FullMat<double>>); }

TEST_CASE("SymmPackMat", "[Matrix.Dense]") { test_dense_mat_setup<double>(create_new<SymmPackMat<double>>); }

TEST_CASE("BandMat", "[Matrix.Dense]") { test_dense_mat_setup<double>(create_new<BandMat<double>>); }

TEST_CASE("BandMatSpike", "[Matrix.Dense]") { test_dense_mat_setup<double>(create_new<BandMatSpike<double>>); }

TEST_CASE("BandSymmMat", "[Matrix.Dense]") { test_dense_mat_setup<double>(create_new<BandSymmMat<double>>); }

TEST_CASE("FullMatFloat", "[Matrix.Dense]") { test_dense_mat_setup<float>(create_new<FullMat<float>>); }

TEST_CASE("SymmPackMatFloat", "[Matrix.Dense]") { test_dense_mat_setup<float>(create_new<SymmPackMat<float>>); }

TEST_CASE("BandMatFloat", "[Matrix.Dense]") { test_dense_mat_setup<float>(create_new<BandMat<float>>); }

TEST_CASE("BandMatSpikeFloat", "[Matrix.Dense]") { test_dense_mat_setup<float>(create_new<BandMatSpike<float>>); }

TEST_CASE("BandSymmMatFloat", "[Matrix.Dense]") { test_dense_mat_setup<float>(create_new<BandSymmMat<float>>); }

TEST_CASE("SparseMatSuperLU", "[Matrix.Sparse]") { test_sparse_mat_setup<double>(create_new<SparseMatSuperLU<double>>); }

TEST_CASE("SparseMatSuperLUFloat", "[Matrix.Sparse]") { test_sparse_mat_setup<float>(create_new<SparseMatSuperLU<float>>); }

TEST_CASE("SparseMatMUMPS", "[Matrix.Sparse]") { test_sparse_mat_setup<double>(create_new<SparseMatMUMPS<double>>); }

TEST_CASE("SparseMatMUMPSFloat", "[Matrix.Sparse]") { test_sparse_mat_setup<float>(create_new<SparseMatMUMPS<float>>); }

TEST_CASE("SparseMatLis", "[Matrix.Sparse]") { test_sparse_mat_setup<double>(create_new<SparseMatLis<double>>); }

#ifdef SUANPAN_MKL
TEST_CASE("SparseMatPARDISO", "[Matrix.Sparse]") { test_sparse_mat_setup<double>(create_new<SparseMatPARDISO<double>>); }

TEST_CASE("SparseMatPARDISOFloat", "[Matrix.Sparse]") { test_sparse_mat_setup<float>(create_new<SparseMatPARDISO<float>>); }

TEST_CASE("SparseMatFGMRES", "[Matrix.Sparse]") { test_sparse_mat_setup<double>(create_new<SparseMatFGMRES<double>>); }

#ifdef SUANPAN_MPI
TEST_CASE("SparseMatMPIPARDISO", "[Matrix.Sparse]") { test_sparse_mat_setup<double>(create_new<SparseMatMPIPARDISO<double>>); }
#endif
#endif

#ifdef SUANPAN_CUDA
TEST_CASE("BandMatCUDA", "[Matrix.Dense]") { test_dense_mat_setup<double>(create_new<BandMatCUDA<double>>); }

TEST_CASE("SparseMatCUDA", "[Matrix.Sparse]") { test_sparse_mat_setup<double>(create_new<SparseMatCUDA<double>>); }

TEST_CASE("SparseMatCUDAFloat", "[Matrix.Sparse]") { test_sparse_mat_setup<float>(create_new<SparseMatCUDA<float>>); }

#ifdef SUANPAN_MAGMA
TEST_CASE("SparseMatMAGMA", "[Matrix.Sparse]") { test_sparse_mat_setup<double>(create_new<SparseMatMAGMA<double>>); }

TEST_CASE("SparseMatMAGMAFloat", "[Matrix.Sparse]") { test_sparse_mat_setup<float>(create_new<SparseMatMAGMA<float>>); }

TEST_CASE("Large CUDA Sparse", "[Matrix.Benchmark]") {
    benchmark_mat_setup<SparseMatPARDISO<double>, double>(0x400);
    benchmark_mat_setup<BandSymmMat<double>, double>(0x2976);
    benchmark_mat_setup<SparseMatMAGMA<double>, double>(0x400);
    benchmark_mat_setup<SparseMatMAGMA<float>, float>(0x400);
}
#endif
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
        // suanpan_info("Assemble: {:.3f}\n", assemble_mean.count() / static_cast<double>(S));
        // suanpan_info("Compress: {:.3f}\n", compress_mean.count() / static_cast<double>(S));
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

TEST_CASE("Aligned Round", "[Matrix.Utility]") {
    REQUIRE(0 == round_up<double>(0));
    REQUIRE(0 == round_up<float>(0));
    REQUIRE(8 == round_up<double>(8));
    REQUIRE(16 == round_up<float>(8));
}
