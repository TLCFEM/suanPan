class glue_solve_symm {
public:
	template<typename T1, typename T2> static void apply(Mat<typename T1::elem_type>& X, const Glue<T1, T2, glue_solve_symm>& S) {
		arma_extra_debug_sigprint();

		if(glue_solve_symm::apply(X, S.A, S.B) == false) arma_stop_runtime_error("solve(): solution not found");
	}

	template<typename eT, typename T1, typename T2> static bool apply(Mat<eT>& X, const T1& A, const T2& B) {
		auto UPLO = 'U';
		auto N = static_cast<int>(A.n_size);
		auto NRHS = static_cast<int>(B.n_cols);
		const auto IPIV = new int[N];
		auto LDB = N;
		auto INFO = 0;

		X = Mat<eT>(B.memptr(), N, NRHS);

		if(is_float<eT>::value) {
			using T = float;
			arma_fortran(arma_sspsv)(&UPLO, &N, &NRHS, (T*)A.memptr(), IPIV, (T*)X.memptr(), &LDB, &INFO);
		}
		else if(is_double<eT>::value) {
			using T = double;
			arma_fortran(arma_dspsv)(&UPLO, &N, &NRHS, (T*)A.memptr(), IPIV, (T*)X.memptr(), &LDB, &INFO);
		}

		delete[] IPIV;

		return INFO == 0;
	}
};
