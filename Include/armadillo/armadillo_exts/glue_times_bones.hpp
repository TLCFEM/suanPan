class glue_times_symm {
public:
	template<typename T1, typename T2> 
	arma_hot
	static

	void apply(Mat<typename T1::elem_type>& out, const Glue<T1, T2, glue_times_symm>& I) {
		typedef typename T1::elem_type eT;
		auto& A = I.A;
		auto& X = I.B;

		out = X;

		auto UPLO = 'U';
		auto N = static_cast<int>(A.n_size);
		eT ALPHA = 1.;
		auto INC = 1;
		eT BETA = 0.;

		if(is_float<eT>::value) {
			using T = float;
			arma_fortran(arma_sspmv)(&UPLO, &N, (T*)&ALPHA, (T*)A.memptr(), (T*)X.memptr(), &INC, (T*)&BETA, (T*)out.memptr(), &INC);
		}
		else if(is_double<eT>::value) {
			using T = double;
			arma_fortran(arma_dspmv)(&UPLO, &N, (T*)&ALPHA, (T*)A.memptr(), (T*)X.memptr(), &INC, (T*)&BETA, (T*)out.memptr(), &INC);
		}
	}
};

class glue_mixed_times_symm {
public:
	template<typename T1, typename T2> inline static void apply(Mat<typename eT_promoter<T1, T2>::eT>& out, const mtGlue<typename eT_promoter<T1, T2>::eT, T1, T2, glue_mixed_times_symm>& X) {
		arma_extra_debug_sigprint();

		typedef typename T1::elem_type eT1;
		typedef typename T2::elem_type eT2;

		const unwrap_check_mixed<T1> tmp1(X.A, out);
		const unwrap_check_mixed<T2> tmp2(X.B, out);

		const Mat<eT1>& A = tmp1.M;
		const Mat<eT2>& B = tmp2.M;

		arma_debug_assert_mul_size(A, B, "matrix multiplication");

		out.set_size(A.n_rows, B.n_cols);

		gemm_mixed<>::apply(out, A, B);
	}
};

template<typename eT> inline Col<eT> sp_mv(const SymmMat<eT>& A, const Col<eT>& X) {
	auto Y = X;

	auto UPLO = 'U';
	auto N = static_cast<int>(A.n_size);
	eT ALPHA = 1.;
	auto INC = 1;
	eT BETA = 0.;

	if(is_float<eT>::value) {
		using T = float;
		arma_fortran(arma_sspmv)(&UPLO, &N, (T*)&ALPHA, (T*)A.memptr(), (T*)X.memptr(), &INC, (T*)&BETA, (T*)Y.memptr(), &INC);
	}
	else if(is_double<eT>::value) {
		using T = double;
		arma_fortran(arma_dspmv)(&UPLO, &N, (T*)&ALPHA, (T*)A.memptr(), (T*)X.memptr(), &INC, (T*)&BETA, (T*)Y.memptr(), &INC);
	}

	return Y;
};
