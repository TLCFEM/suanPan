template<typename T1, typename T2> typename enable_if2<is_supported_blas_type < typename T1::elem_type>::value &&is_SymmMat<T1>::value,


const Glue<T1, T2, glue_solve_symm>
>
::result solve(const T1& A, const Base<typename T1::elem_type, T2>& B) {
	arma_extra_debug_sigprint();

	return Glue<T1, T2, glue_solve_symm>(A.get_ref(), B.get_ref());
}
