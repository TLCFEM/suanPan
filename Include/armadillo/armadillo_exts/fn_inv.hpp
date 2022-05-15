template<typename T1> arma_warn_unused

typename enable_if2<is_supported_blas_type<typename T1::elem_type>::value && is_SymmMat<T1>::value, const SmOp<T1, smop_inv>>::result inv(const Base<typename T1::elem_type, T1>& X) {
	arma_extra_debug_sigprint();

	return SmOp<T1, smop_inv>(X.get_ref());
}
