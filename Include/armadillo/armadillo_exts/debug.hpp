template<typename eT1, typename eT2> 
arma_hot

void arma_assert_same_size(const SymmMat<eT1>& A, const SymmMat<eT2>& B, const char* x) {
	const uword A_n_size = A.n_size;
	const uword B_n_size = B.n_size;

	if(A_n_size != B_n_size) arma_stop_logic_error(arma_incompat_size_string(A_n_size, A_n_size, B_n_size, B_n_size, x));
}
