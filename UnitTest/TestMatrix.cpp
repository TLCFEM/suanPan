#include "CatchHeader.h"

#include <Domain/MetaMat/MetaMat>

namespace {
    template<typename MT, typename ET, std::invocable T> void test_mat_solve(MT& A, const Mat<ET>& D, const Col<ET>& C, T clear_mat) {
        constexpr auto tol = std::numeric_limits<ET>::epsilon() * 100;
        const auto scaled_tol = static_cast<ET>(C.n_elem) * tol;

        Col<ET> E(C.n_elem);

        clear_mat();

        // full solve
        A.solve(E, C);
        REQUIRE(arma::norm<Col<ET>>(E - D, "inf") < scaled_tol);

        // factored solve
        A.solve(E, C);
        REQUIRE(arma::norm<Col<ET>>(E - D, "inf") < scaled_tol);

        clear_mat();

        // r-value full solve
        A.solve(E, Mat<ET>(C));
        REQUIRE(arma::norm<Col<ET>>(E - D, "inf") < scaled_tol);

        // r-value factored solve
        A.solve(E, Mat<ET>(C));
        REQUIRE(arma::norm<Col<ET>>(E - D, "inf") < scaled_tol);

        // mixed precision
        A.get_solver_setting().precision = Precision::MIXED;
        A.get_solver_setting().tolerance = tol;

        clear_mat();

        // full solve
        A.solve(E, C);
        REQUIRE(arma::norm<Col<ET>>(E - D, "inf") < scaled_tol);

        // factored solve
        A.solve(E, C);
        REQUIRE(arma::norm<Col<ET>>(E - D, "inf") < scaled_tol);

        clear_mat();

        // r-value full solve
        A.solve(E, Mat<ET>(C));
        REQUIRE(arma::norm<Col<ET>>(E - D, "inf") < scaled_tol);

        // r-value factored solve
        A.solve(E, Mat<ET>(C));
        REQUIRE(arma::norm<Col<ET>>(E - D, "inf") < scaled_tol);
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
                for(uword i{0}; i < N; ++i)
                    for(uword j{0}; j < N; ++j)
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

            test_mat_solve(A, spsolve(B, C,
#ifdef SUANPAN_SUPERLUMT
                                      "lapack"
#else
                                      "superlu"
#endif
                              ),
                           C, clear_mat);
        }
    }

    template<typename MT, typename ET, std::invocable T> void benchmark_mat_solve(std::string&& title, MT& A, const Col<ET>& C, const Mat<ET>& E, T&& clear_mat) {
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

    template<typename T> T create_new(uword) { throw std::runtime_error("unknown matrix"); }

    template<> FullMat<double> create_new(const uword N) { return {N, N}; }

    template<> SymmPackMat<double> create_new(const uword N) { return SymmPackMat<double>{N}; }

    template<> BandMat<double> create_new(const uword N) { return {N, std::max(N / uword{200}, uword{3}), std::max(N / uword{200}, uword{3})}; }

    template<> BandMatSpike<double> create_new(const uword N) { return {N, std::max(N / uword{200}, uword{3}), std::max(N / uword{200}, uword{3})}; }

    template<> BandSymmMat<double> create_new(const uword N) { return {N, std::max(N / uword{200}, uword{3})}; }

    template<> SparseMatSuperLU<double> create_new(const uword N) { return {N, N}; }

    template<> FullMat<float> create_new(const uword N) { return {N, N}; }

    template<> SymmPackMat<float> create_new(const uword N) { return SymmPackMat<float>{N}; }

    template<> BandMat<float> create_new(const uword N) { return {N, std::max(N / uword{200}, uword{3}), std::max(N / uword{200}, uword{3})}; }

    template<> BandMatSpike<float> create_new(const uword N) { return {N, std::max(N / uword{200}, uword{3}), std::max(N / uword{200}, uword{3})}; }

    template<> BandSymmMat<float> create_new(const uword N) { return {N, std::max(N / uword{200}, uword{3})}; }

    template<> SparseMatSuperLU<float> create_new(const uword N) { return {N, N}; }

#ifdef SUANPAN_MKL
    template<> SparseMatPARDISO<double> create_new(const uword N) { return {N, N}; }

    template<> SparseMatPARDISO<float> create_new(const uword N) { return {N, N}; }

#ifdef SUANPAN_DISTRIBUTED
    template<> SparseMatClusterPARDISO<double> create_new(const uword N) { return {N, N}; }

    template<> SparseMatClusterPARDISO<float> create_new(const uword N) { return {N, N}; }
#endif

    template<> SparseMatFGMRES<double> create_new(const uword N) { return {N, N}; }
#endif

#ifdef SUANPAN_CUDA
    template<> FullMatCUDA<double> create_new(const uword N) { return {N, N}; }

    template<> FullMatCUDA<float> create_new(const uword N) { return {N, N}; }

    template<> SparseMatCUDA<double> create_new(const uword N) { return {N, N}; }

    template<> SparseMatCUDA<float> create_new(const uword N) { return {N, N}; }

#ifdef SUANPAN_MAGMA
    template<> BandMatMAGMA<double> create_new(const uword N) { return {N, std::max(N / uword{200}, uword{3}), std::max(N / uword{200}, uword{3})}; }

    template<> BandMatMAGMA<float> create_new(const uword N) { return {N, std::max(N / uword{200}, uword{3}), std::max(N / uword{200}, uword{3})}; }

    template<> SparseMatMAGMA<double> create_new(const uword N) { return {N, N}; }

    template<> SparseMatMAGMA<float> create_new(const uword N) { return {N, N}; }
#endif
#endif

    template<typename T, typename ET> void benchmark_mat_setup(const int I) {
        const auto C = randu<Col<ET>>(I);

        Mat<ET> V(I, 5, fill::ones);
        V.col(2) += 10 * C + 10;

        auto B = spdiags(V, ivec{-2, -1, 0, +1, +2}, I, I);

        auto A = create_new<T>(I);

        std::string title;

        if(std::is_same_v<FullMat<ET>, T>) title = "Full ";
        else if(std::is_same_v<SymmPackMat<ET>, T>) title = "SymmPack ";
        else if(std::is_same_v<BandMat<ET>, T>) title = "Band ";
        else if(std::is_same_v<BandMatSpike<ET>, T>) title = "BandSpike ";
        else if(std::is_same_v<BandSymmMat<ET>, T>) title = "BandSymm ";
        else if(std::is_same_v<SparseMatSuperLU<ET>, T>) title = "SuperLU ";
#ifdef SUANPAN_MKL
        else if(std::is_same_v<SparseMatPARDISO<ET>, T>) title = "PARDISO ";
#ifdef SUANPAN_DISTRIBUTED
        else if(std::is_same_v<SparseMatClusterPARDISO<ET>, T>) title = "Cluster PARDISO ";
        else if(std::is_same_v<SparseSymmMatClusterPARDISO<ET>, T>) title = "Cluster Symm PARDISO ";
        else if(std::is_same_v<SparseSPDMatClusterPARDISO<ET>, T>) title = "Cluster SPD PARDISO ";
#endif
        else if(std::is_same_v<SparseMatFGMRES<ET>, T>) title = "FGMRES ";
#endif
#ifdef SUANPAN_CUDA
        else if(std::is_same_v<FullMatCUDA<ET>, T>) title = "Full CUDA ";
        else if(std::is_same_v<SparseMatCUDA<ET>, T>) title = "Sparse CUDA ";
#ifdef SUANPAN_MAGMA
        else if(std::is_same_v<BandMatMAGMA<ET>, T>) title = "Band Magma ";
        else if(std::is_same_v<SparseMatMAGMA<ET>, T>) title = "Sparse Magma ";
#endif
#endif

        title += "N=" + std::to_string(I) + " NZ=" + std::to_string(B.n_nonzero) + " NE=" + std::to_string(A.n_elem);

        benchmark_mat_solve(std::move(title), A, C, spsolve(B, C,
#ifdef SUANPAN_SUPERLUMT
                                                            "lapack"
#else
                                                            "superlu"
#endif
                                                    ),
                            [&] {
                                A.zeros();
                                for(auto J = B.begin(); J != B.end(); ++J) A.at(J.col(), J.row()) = *J;
                            });
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
} // namespace

TEST_CASE("Mixed Precision", "[Matrix.Benchmark]") {
    for(auto I = 0x0020; I < 0x0100; I *= 2) {
        benchmark_mat_setup<FullMat<double>, double>(I);
        benchmark_mat_setup<SymmPackMat<double>, double>(I);
        benchmark_mat_setup<BandMat<double>, double>(I);
        benchmark_mat_setup<BandMatSpike<double>, double>(I);
        benchmark_mat_setup<BandSymmMat<double>, double>(I);
        benchmark_mat_setup<SparseMatSuperLU<double>, double>(I);
#ifdef SUANPAN_MKL
        benchmark_mat_setup<SparseMatPARDISO<double>, double>(I);
        benchmark_mat_setup<SparseMatFGMRES<double>, double>(I);
#endif
#ifdef SUANPAN_CUDA
        benchmark_mat_setup<SparseMatCUDA<double>, double>(I);
#endif
    }
#ifdef SUANPAN_CUDA
    for(auto I = 0x0100; I < 0x2000; I *= 2) benchmark_mat_setup<FullMatCUDA<double>, double>(I);
#endif
}

TEST_CASE("Large Mixed Precision", "[Matrix.Benchmark]") {
    for(auto I = 0x400; I < 0x500; I *= 2) {
        benchmark_mat_setup<BandMat<double>, double>(I);
        benchmark_mat_setup<BandMatSpike<double>, double>(I);
        benchmark_mat_setup<BandSymmMat<double>, double>(I);
        benchmark_mat_setup<SparseMatSuperLU<double>, double>(I);
    }
}

TEST_CASE("Large Sparse Solve Type", "[Matrix.Benchmark]") {
    for(auto I = 0x1000; I < 0x5000; I *= 2) {
        benchmark_mat_setup<BandMat<double>, double>(I);
        benchmark_mat_setup<SparseMatSuperLU<double>, double>(I);
#ifdef SUANPAN_MKL
        benchmark_mat_setup<SparseMatPARDISO<double>, double>(I);
        benchmark_mat_setup<SparseMatFGMRES<double>, double>(I);
#endif
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

#ifdef SUANPAN_MKL
TEST_CASE("SparseMatPARDISO", "[Matrix.Sparse]") { test_sparse_mat_setup<double>(create_new<SparseMatPARDISO<double>>); }

TEST_CASE("SparseMatPARDISOFloat", "[Matrix.Sparse]") { test_sparse_mat_setup<float>(create_new<SparseMatPARDISO<float>>); }

TEST_CASE("SparseMatFGMRES", "[Matrix.Sparse]") { test_sparse_mat_setup<double>(create_new<SparseMatFGMRES<double>>); }

#ifdef SUANPAN_DISTRIBUTED
TEST_CASE("SparseMatClusterPARDISO", "[Matrix.Sparse]") { test_sparse_mat_setup<double>(create_new<SparseMatClusterPARDISO<double>>); }

TEST_CASE("SparseMatClusterPARDISOFloat", "[Matrix.Sparse]") { test_sparse_mat_setup<float>(create_new<SparseMatClusterPARDISO<float>>); }
#endif
#endif

#ifdef SUANPAN_CUDA
TEST_CASE("SparseMatCUDA", "[Matrix.Sparse]") { test_sparse_mat_setup<double>(create_new<SparseMatCUDA<double>>); }

TEST_CASE("SparseMatCUDAFloat", "[Matrix.Sparse]") { test_sparse_mat_setup<float>(create_new<SparseMatCUDA<float>>); }

#ifdef SUANPAN_MAGMA
TEST_CASE("BandMatMAGMA", "[Matrix.Dense]") { test_dense_mat_setup<double>(create_new<BandMatMAGMA<double>>); }

TEST_CASE("BandMatMAGMAFloat", "[Matrix.Dense]") { test_dense_mat_setup<float>(create_new<BandMatMAGMA<float>>); }

TEST_CASE("SparseMatMAGMA", "[Matrix.Sparse]") { test_sparse_mat_setup<double>(create_new<SparseMatMAGMA<double>>); }

TEST_CASE("SparseMatMAGMAFloat", "[Matrix.Sparse]") { test_sparse_mat_setup<float>(create_new<SparseMatMAGMA<float>>); }

TEST_CASE("Large CUDA Sparse", "[Matrix.Benchmark]") {
    for(auto I = 0x4000; I < 0x10000; I *= 2) {
        benchmark_mat_setup<BandSymmMat<double>, double>(I);
        benchmark_mat_setup<BandMatMAGMA<double>, double>(I);
        benchmark_mat_setup<SparseMatPARDISO<double>, double>(I);
        benchmark_mat_setup<SparseMatMAGMA<double>, double>(I);
        benchmark_mat_setup<SparseMatCUDA<double>, double>(I);
        benchmark_mat_setup<BandSymmMat<float>, float>(I);
        benchmark_mat_setup<BandMatMAGMA<float>, float>(I);
        benchmark_mat_setup<SparseMatPARDISO<float>, float>(I);
        benchmark_mat_setup<SparseMatMAGMA<float>, float>(I);
        benchmark_mat_setup<SparseMatCUDA<float>, float>(I);
    }
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

        B.assemble(A, linspace<uvec>(0, N - 1u, N));

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
        BENCHMARK(std::string("Assemble " + std::to_string(J)).c_str()) {
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

namespace {
    template<typename T> void check_compressed_ptr(const T* ptr, const uword n, const uword n_elem, const uword shift) {
        REQUIRE(static_cast<uword>(ptr[0]) == shift);
        for(uword I = 1; I <= n; ++I) REQUIRE(ptr[I] >= ptr[I - 1]);
        REQUIRE(static_cast<uword>(ptr[n]) == n_elem + shift);
    }

    template<typename F> mat csr_to_dense(const F& in_mat, const uword shift) {
        mat out_mat(static_cast<uword>(in_mat.n_rows), static_cast<uword>(in_mat.n_cols), fill::zeros);
        for(uword I = 0; I < out_mat.n_rows; ++I)
            for(auto J = static_cast<uword>(in_mat.row_mem()[I]) - shift; J < static_cast<uword>(in_mat.row_mem()[I + 1]) - shift; ++J) out_mat(I, static_cast<uword>(in_mat.col_mem()[J]) - shift) += in_mat.val_mem()[J];
        return out_mat;
    }

    template<typename F> mat csc_to_dense(const F& in_mat, const uword shift) {
        mat out_mat(static_cast<uword>(in_mat.n_rows), static_cast<uword>(in_mat.n_cols), fill::zeros);
        for(uword I = 0; I < out_mat.n_cols; ++I)
            for(auto J = static_cast<uword>(in_mat.col_mem()[I]) - shift; J < static_cast<uword>(in_mat.col_mem()[I + 1]) - shift; ++J) out_mat(static_cast<uword>(in_mat.row_mem()[J]) - shift, I) += in_mat.val_mem()[J];
        return out_mat;
    }

    template<typename T> void check_conversion(const mat& A) {
        const auto R = A.n_rows, C = A.n_cols;

        triplet_form<double, uword> B(R, C);

        // insert every entry twice to exercise condensation
        for(uword I = 0; I < R; ++I)
            for(uword J = 0; J < C; ++J)
                if(A(I, J) != 0.) {
                    B.at(I, J) = .5 * A(I, J);
                    B.at(I, J) = .5 * A(I, J);
                }

        const vec X(C, fill::randn);

        for(const auto base : {SparseBase::ZERO, SparseBase::ONE})
            for(const auto full : {false, true}) {
                const auto shift = static_cast<uword>(base);

                auto BB = B;
                csr_form<double, T> D(BB, base, full);
                BB = B;
                csc_form<double, T> E(BB, base, full);

                REQUIRE(static_cast<uword>(D.n_rows) == R);
                REQUIRE(static_cast<uword>(D.n_cols) == C);
                REQUIRE(static_cast<uword>(E.n_rows) == R);
                REQUIRE(static_cast<uword>(E.n_cols) == C);

                check_compressed_ptr(D.row_mem(), R, static_cast<uword>(D.n_elem), shift);
                check_compressed_ptr(E.col_mem(), C, static_cast<uword>(E.n_elem), shift);

                REQUIRE(norm(A - csr_to_dense(D, shift), "inf") == 0.);
                REQUIRE(norm(A - csc_to_dense(E, shift), "inf") == 0.);

                if(SparseBase::ZERO == base) {
                    REQUIRE(norm(A - to_mat(D), "inf") == 0.);
                    REQUIRE(norm(A - to_mat(E), "inf") == 0.);
                    for(uword I = 0; I < R; ++I)
                        for(uword J = 0; J < C; ++J) {
                            REQUIRE(D(static_cast<T>(I), static_cast<T>(J)) == A(I, J));
                            REQUIRE(E(static_cast<T>(I), static_cast<T>(J)) == A(I, J));
                        }
                    if(R == C) {
                        REQUIRE(norm(A * X - D * X, "inf") <= 1E-13);
                        REQUIRE(norm(A * X - E * X, "inf") <= 1E-13);
                    }
                }
            }

        // assignment operator path
        auto BB = B;
        csr_form<double, T> D;
        D = BB;
        BB = B;
        csc_form<double, T> E;
        E = BB;

        check_compressed_ptr(D.row_mem(), R, static_cast<uword>(D.n_elem), 0);
        check_compressed_ptr(E.col_mem(), C, static_cast<uword>(E.n_elem), 0);

        REQUIRE(norm(A - to_mat(D), "inf") == 0.);
        REQUIRE(norm(A - to_mat(E), "inf") == 0.);
    }
} // namespace

TEST_CASE("Triplet/CSR/CSC Conversion Empty Rows/Columns", "[Matrix.Sparse]") {
    // leading, trailing and interior empty rows/columns in a rectangular matrix
    mat A(7, 5, fill::zeros);
    A(1, 1) = 1.;
    A(1, 3) = 2.;
    A(3, 1) = 3.;
    A(4, 3) = 4.;

    check_conversion<uword>(A);
    check_conversion<int>(A);

    // transposed shape so columns outnumber rows
    check_conversion<uword>(mat(A.t()));
    check_conversion<int>(mat(A.t()));

    // only the first entry is populated
    mat F(4, 6, fill::zeros);
    F(0, 0) = 1.;

    check_conversion<uword>(F);
    check_conversion<int>(F);

    // only the last entry is populated
    mat L(6, 4, fill::zeros);
    L(5, 3) = 1.;

    check_conversion<uword>(L);
    check_conversion<int>(L);

    // fully empty matrix
    check_conversion<uword>(mat(5, 3, fill::zeros));
    check_conversion<int>(mat(3, 5, fill::zeros));

    // random sparse matrices with many empty rows/columns
    for(auto I = 0; I < 50; ++I) {
        check_conversion<uword>(mat(sprandu<sp_mat>(50, 50, .02)));
        check_conversion<int>(mat(sprandu<sp_mat>(40, 60, .02)));
        check_conversion<uword>(mat(sprandu<sp_mat>(60, 40, .02)));
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

TEST_CASE("Unify FullMat", "[Matrix.Utility]") { test_dense_mat_unify(FullMat<double>(10, 10)); }

TEST_CASE("Unify BandMat", "[Matrix.Utility]") { test_dense_mat_unify(BandMat<double>(10, 2, 3)); }

TEST_CASE("Unify BandSymmMat", "[Matrix.Utility]") { test_dense_mat_unify(BandSymmMat<double>(10, 2)); }

TEST_CASE("Unify BandMatSpike", "[Matrix.Utility]") { test_dense_mat_unify(BandMatSpike<double>(10, 2, 3)); }

TEST_CASE("Unify SymmPackMat", "[Matrix.Utility]") { test_dense_mat_unify(SymmPackMat<double>(10)); }

TEST_CASE("Unify SparseMatSuperLU", "[Matrix.Utility]") { test_sparse_mat_unify(SparseMatSuperLU<double>(10, 10)); }

#ifdef SUANPAN_MKL
TEST_CASE("Unify SparseMatPARDISO", "[Matrix.Utility]") { test_sparse_mat_unify(SparseMatPARDISO<double>(10, 10)); }
#endif

#ifdef SUANPAN_CUDA
TEST_CASE("Unify SparseMatCUDA", "[Matrix.Utility]") { test_sparse_mat_unify(SparseMatCUDA<double>(10, 10)); }

#ifdef SUANPAN_MAGMA
TEST_CASE("Unify BandMatMAGMA", "[Matrix.Utility]") { test_dense_mat_unify(BandMatMAGMA<double>(10, 2, 3)); }
#endif
#endif

TEST_CASE("Aligned Round", "[Matrix.Utility]") {
    REQUIRE(0 == round_up<double>(0));
    REQUIRE(0 == round_up<float>(0));
    REQUIRE(8 == round_up<double>(8));
    REQUIRE(16 == round_up<float>(8));
}
