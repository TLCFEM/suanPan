class smop_inv {
public:
	template<typename T1> static void apply(T1& out, const SmOp<T1, smop_inv>& in) {
		arma_extra_debug_sigprint();

		out = in.m;

		if(sp_inv(out) != 0) {
			out.reset();
			arma_stop_runtime_error("sp_inv(): matrix seems singular");
		}
	}
};

template<typename eT> int sp_inv(SymmMat<eT>& A) {
	auto UPLO = 'U';
	auto N = static_cast<int>(A.n_size);
	const auto IPIV = new int[N];
	auto INFO = 0;

	if(is_float<eT>::value) {
		using T = float;
		arma_fortran(arma_ssptrf)(&UPLO, &N, (T*)A.memptr(), IPIV, &INFO);
	}
	else if(is_double<eT>::value) {
		using T = double;
		arma_fortran(arma_dsptrf)(&UPLO, &N, (T*)A.memptr(), IPIV, &INFO);
	}

	if(INFO != 0) return INFO;

	const auto WORK = new eT[N];

	if(is_float<eT>::value) {
		using T = float;
		arma_fortran(arma_ssptri)(&UPLO, &N, (T*)A.memptr(), IPIV, (T*)WORK, &INFO);
	}
	else if(is_double<eT>::value) {
		using T = double;
		arma_fortran(arma_dsptri)(&UPLO, &N, (T*)A.memptr(), IPIV, (T*)WORK, &INFO);
	}

	delete[] WORK;
	delete[] IPIV;

	return INFO;
}
