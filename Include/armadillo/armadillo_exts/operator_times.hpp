template<typename T1> typename enable_if2<is_SymmMat<T1>::value, const eOp<T1, eop_scalar_times>>::result operator*(const T1& X, const typename T1::elem_type k) {
	arma_extra_debug_sigprint();

	return eOp<T1, eop_scalar_times>(X, k);
}

template<typename T1> typename enable_if2<is_SymmMat<T1>::value, const eOp<T1, eop_scalar_times>>::result operator*(const typename T1::elem_type k, const T1& X) {
	arma_extra_debug_sigprint();

	return eOp<T1, eop_scalar_times>(X, k);
}

template<typename T1, typename T2> typename enable_if2<is_SymmMat<T1>::value && is_Mat<T2>::value && is_same_type < typename T1::elem_type, typename T2::elem_type>::value
,
const Glue<T1, T2, glue_times_symm>
>
::result operator*(const T1& X, const T2& Y) {
	arma_extra_debug_sigprint();

	return Glue<T1, T2, glue_times_symm>(X, Y);
}

template<typename T1, typename T2> typename enable_if2<is_Mat<T1>::value && is_SymmMat<T2>::value && is_same_type < typename T1::elem_type, typename T2::elem_type>::value
,
const Glue<T1, T2, glue_times_symm>
>
::result operator*(const T1& X, const T2& Y) {
	arma_extra_debug_sigprint();

	return Glue<T1, T2, glue_times_symm>(X, Y);
}

template<typename T1, typename T2> inline typename enable_if2<(is_SymmMat<T1>::value && is_arma_type<T2>::value && (is_same_type < typename T1::elem_type, typename T2::elem_type > ::no)), const mtGlue<typename promote_type < typename T1::elem_type, typename T2::elem_type>::result, T1, T2, glue_mixed_times_symm>
>
::result operator*(const T1& X, const T2& Y) {
	arma_extra_debug_sigprint();

	typedef typename T1::elem_type eT1;
	typedef typename T2::elem_type eT2;

	typedef typename promote_type<eT1, eT2>::result out_eT;

	promote_type<eT1, eT2>::check();

	return mtGlue<out_eT, T1, T2, glue_mixed_times_symm>(X, Y);
}

template<typename T1, typename T2> inline typename enable_if2<(is_arma_type<T1>::value && is_SymmMat<T2>::value && (is_same_type < typename T1::elem_type, typename T2::elem_type > ::no)), const mtGlue<typename promote_type < typename T1::elem_type, typename T2::elem_type>::result, T1, T2, glue_mixed_times_symm>
>
::result operator*(const T1& X, const T2& Y) {
	arma_extra_debug_sigprint();

	typedef typename T1::elem_type eT1;
	typedef typename T2::elem_type eT2;

	typedef typename promote_type<eT1, eT2>::result out_eT;

	promote_type<eT1, eT2>::check();

	return mtGlue<out_eT, T1, T2, glue_mixed_times_symm>(X, Y);
}
