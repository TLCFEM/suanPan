template<typename T1, typename bdop_type> BdOp<T1, bdop_type>::BdOp(const T1& in_m)
	: m(in_m) { arma_extra_debug_sigprint(); }

template<typename T1, typename bdop_type> BdOp<T1, bdop_type>::BdOp(const T1& in_m, const typename T1::elem_type in_aux)
	: m(in_m)
	, aux(in_aux) { arma_extra_debug_sigprint(); }

template<typename T1, typename bdop_type> BdOp<T1, bdop_type>::BdOp(const T1& in_m, const uword in_aux_uword_a, const uword in_aux_uword_b)
	: m(in_m)
	, aux_uword_a(in_aux_uword_a)
	, aux_uword_b(in_aux_uword_b) { arma_extra_debug_sigprint(); }

template<typename T1, typename bdop_type> BdOp<T1, bdop_type>::~BdOp() { arma_extra_debug_sigprint(); }
