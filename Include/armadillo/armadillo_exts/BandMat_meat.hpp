template<typename eT> BandMat<eT>::~BandMat() {}

template<typename eT> BandMat<eT>::BandMat()
	: n_cols(0)
	, n_l(0)
	, n_u(0)
	, n_elem(0)
	, mem_state(0)
	, mem() {}

template<typename eT> BandMat<eT>::BandMat(const uword& in_size, const uword& in_l, const uword& in_u)
	: n_cols(in_size)
	, n_l(in_l)
	, n_u(in_u)
	, n_elem((in_u + 2 * in_l + 1) * in_size)
	, mem_state(0)
	, mem() {
	arma_extra_debug_sigprint_this(this);
	init_cold();
}

template<typename eT> template<typename fill_type> BandMat<eT>::BandMat(const uword& in_size, const uword& in_l, const uword& in_u, const fill::fill_class<fill_type>& f)
	: n_cols(in_size)
	, n_l(in_l)
	, n_u(in_u)
	, n_elem((in_u + 2 * in_l + 1) * in_size)
	, mem_state(0)
	, mem() {
	arma_extra_debug_sigprint_this(this);

	init_cold();

	(*this).fill(f);
}

template<typename eT> BandMat<eT>& BandMat<eT>::operator=(const eT& val) {
	init_warm(1, 0, 0);
	access::rw(mem[0]) = val;
	return *this;
}

template<typename eT> BandMat<eT>& BandMat<eT>::operator+=(const eT& val) {
	arma_extra_debug_sigprint();

	arrayops::inplace_plus(memptr(), val, n_elem);

	return *this;
}

template<typename eT> BandMat<eT>& BandMat<eT>::operator-=(const eT& val) {
	arma_extra_debug_sigprint();

	arrayops::inplace_minus(memptr(), val, n_elem);

	return *this;
}

template<typename eT> BandMat<eT>& BandMat<eT>::operator*=(const eT& val) {
	arma_extra_debug_sigprint();

	arrayops::inplace_mul(memptr(), val, n_elem);

	return *this;
}

template<typename eT> BandMat<eT>& BandMat<eT>::operator/=(const eT& val) {
	arma_extra_debug_sigprint();

	arrayops::inplace_div(memptr(), val, n_elem);

	return *this;
}

template<typename eT> BandMat<eT>::BandMat(const BandMat& m)
	: n_cols(m.n_cols)
	, n_elem(m.n_elem)
	, mem_state(0)
	, mem() {
	arma_extra_debug_sigprint(arma_str::format("this = %x   in_mat = %x") % this % &m);

	init_cold();

	arrayops::copy(memptr(), m.mem, m.n_elem);
}

template<typename eT> BandMat<eT>& BandMat<eT>::operator=(const BandMat& m) {
	arma_extra_debug_sigprint(arma_str::format("this = %x   in_mat = %x") % this % &m);

	if(this != &m) {
		init_warm(m.n_cols, m.n_l, m.n_u);

		arrayops::copy(memptr(), m.mem, m.n_elem);
	}

	return *this;
}

template<typename eT> BandMat<eT>& BandMat<eT>::operator+=(const BandMat& m) {
	arma_extra_debug_sigprint();

	arma_debug_assert_same_size(*this, m, "addition");

	arrayops::inplace_plus(memptr(), m.memptr(), n_elem);

	return *this;
}

template<typename eT> BandMat<eT>& BandMat<eT>::operator-=(const BandMat& m) {
	arma_extra_debug_sigprint();

	arma_debug_assert_same_size(*this, m, "subtraction");

	arrayops::inplace_minus(memptr(), m.memptr(), n_elem);

	return *this;
}

template<typename eT> BandMat<eT>& BandMat<eT>::operator%=(const BandMat& m) {
	arma_extra_debug_sigprint();

	arma_debug_assert_same_size(*this, m, "element-wise multiplication");

	arrayops::inplace_mul(memptr(), m.memptr(), n_elem);

	return *this;
}

template<typename eT> BandMat<eT>& BandMat<eT>::operator/=(const BandMat& m) {
	arma_extra_debug_sigprint();

	arma_debug_assert_same_size(*this, m, "element-wise division");

	arrayops::inplace_div(memptr(), m.memptr(), n_elem);

	return *this;
}

template<typename eT> template<typename T1, typename bdop_type> BandMat<eT>::BandMat(const BdOp<T1, bdop_type>& X)
	: n_cols(0)
	, n_l(0)
	, n_u(0)
	, n_elem(0)
	, mem_state(0)
	, mem() {
	arma_extra_debug_sigprint_this(this);

	arma_type_check((is_same_type<eT, typename T1::elem_type>::no));

	bdop_type::apply(*this, X);
}

template<typename eT> template<typename T1, typename bdop_type> BandMat<eT>& BandMat<eT>::operator=(const BdOp<T1, bdop_type>& X) {
	arma_extra_debug_sigprint();

	arma_type_check((is_same_type<eT, typename T1::elem_type>::no));

	bdop_type::apply(*this, X);

	return *this;
}

template<typename eT> eT& BandMat<eT>::at(const uword& in_row, const uword& in_col) { return access::rw(mem[n_a + n_b * in_col + in_row]); }

template<typename eT> const eT& BandMat<eT>::at(const uword& in_row, const uword& in_col) const { return mem[n_a + n_b * in_col + in_row]; }

template<typename eT> eT& BandMat<eT>::operator()(const uword& in_row, const uword& in_col) {
	arma_debug_check(in_row >= n_cols || in_col >= n_cols, "BandMat::operator(): index out of bounds");
	return at(in_row, in_col);
}

template<typename eT> const eT& BandMat<eT>::operator()(const uword& in_row, const uword& in_col) const {
	arma_debug_check(in_row >= n_cols || in_col >= n_cols, "BandMat::operator(): index out of bounds");
	return at(in_row, in_col);
}

template<typename eT> void BandMat<eT>::init_cold() {
	arma_extra_debug_sigprint(arma_str::format("n_size = %d") % n_cols);

#if(defined(ARMA_USE_CXX11) || defined(ARMA_64BIT_WORD))
    auto error_message = "BandMat::init(): requested size is too large";
#else
	const char* error_message = "BandMat::init(): requested size is too large; suggest to "
		"compile in C++11 mode or enable ARMA_64BIT_WORD";
#endif

	arma_debug_check(n_cols > ARMA_MAX_UHWORD ? n_elem > ARMA_MAX_UWORD : false, error_message);

	if(n_elem <= arma_config::mat_prealloc)
		if(n_elem == 0) access::rw(mem) = NULL;
		else {
			arma_extra_debug_print("BandMat::init(): using local memory");
			access::rw(mem) = mem_local;
		}
	else {
		arma_extra_debug_print("BandMat::init(): acquiring memory");
		access::rw(mem) = memory::acquire<eT>(n_elem);
	}
}

template<typename eT> void BandMat<eT>::init_warm(const uword& in_size, const uword& in_l, const uword& in_u) {
	arma_extra_debug_sigprint(arma_str::format("in_n_size = %d") % in_size);

	if(n_cols == in_size && n_l == in_l && n_u == in_u) return;

	auto err_state = false;
	char* err_msg = nullptr;

	const uhword t_mem_state = mem_state;

	arma_debug_set_error(err_state, err_msg, t_mem_state == 3, "BandMat::init(): size is fixed and hence cannot be changed");

#if(defined(ARMA_USE_CXX11) || defined(ARMA_64BIT_WORD))
    auto error_message = "BandMat::init(): requested size is too large";
#else
	const char* error_message = "BandMat::init(): requested size is too large; suggest "
		"to compile in C++11 mode or enable ARMA_64BIT_WORD";
#endif

	arma_debug_set_error(err_state, err_msg, in_size > ARMA_MAX_UHWORD ? (2 * in_l + in_u + 1) * in_size > ARMA_MAX_UWORD : false, error_message);

	arma_debug_check(err_state, err_msg);

	const auto old_n_elem = n_elem;
	const auto new_n_elem = (2 * in_l + in_u + 1) * in_size;

	if(old_n_elem == new_n_elem) { arma_extra_debug_print("BandMat::init(): reusing memory"); }
	else {
		arma_debug_check(t_mem_state == 2,
		                 "BandMat::init(): mismatch between size of "
		                 "auxiliary memory and requested size");

		if(new_n_elem < old_n_elem) {
			if(t_mem_state == 0 && new_n_elem <= arma_config::mat_prealloc) {
				if(old_n_elem > arma_config::mat_prealloc) {
					arma_extra_debug_print("BandMat::init(): releasing memory");
					memory::release(access::rw(mem));
				}

				access::rw(mem) = new_n_elem == 0 ? NULL : mem_local;
			}
			else { arma_extra_debug_print("BandMat::init(): reusing memory"); }
		}
		else {
			if(t_mem_state == 0 && old_n_elem > arma_config::mat_prealloc) {
				arma_extra_debug_print("BandMat::init(): releasing memory");
				memory::release(access::rw(mem));
			}

			access::rw(mem) = new_n_elem <= arma_config::mat_prealloc ? mem_local : memory::acquire<eT>(new_n_elem);

			access::rw(mem_state) = 0;
		}

		access::rw(n_cols) = in_size;
		access::rw(n_l) = in_l;
		access::rw(n_u) = in_u;
		access::rw(n_s) = n_l + n_u;
		access::rw(n_rows) = 2 * n_l + n_u + 1;
		access::rw(n_a) = n_rows - n_l;
		access::rw(n_b) = n_rows - 1;
		access::rw(n_elem) = new_n_elem;
	}
}

template<typename eT> eT* BandMat<eT>::memptr() { return const_cast<eT*>(mem); }

template<typename eT> const eT* BandMat<eT>::memptr() const { return mem; }

template<typename eT> void BandMat<eT>::set_size(const uword in_size, const uword& in_l, const uword& in_u) {
	arma_extra_debug_sigprint();

	init_warm(in_size, in_l, in_u);
}

template<typename eT> const BandMat<eT>& BandMat<eT>::fill(const eT val) {
	arma_extra_debug_sigprint();

	arrayops::inplace_set(memptr(), val, n_elem);

	return *this;
}

template<typename eT> template<typename fill_type> const BandMat<eT>& BandMat<eT>::fill(const fill::fill_class<fill_type>&) {
	arma_extra_debug_sigprint();

	if(is_same_type<fill_type, fill::fill_zeros>::yes) (*this).zeros();
	if(is_same_type<fill_type, fill::fill_ones>::yes) (*this).ones();
	if(is_same_type<fill_type, fill::fill_eye>::yes) (*this).eye();

	return *this;
}

template<typename eT> const BandMat<eT>& BandMat<eT>::zeros() {
	arma_extra_debug_sigprint();

	arrayops::fill_zeros(memptr(), n_elem);

	return *this;
}

template<typename eT> const BandMat<eT>& BandMat<eT>::zeros(const uword in_size, const uword& in_l, const uword& in_u) {
	arma_extra_debug_sigprint();

	set_size(in_size, in_l, in_u);

	return (*this).zeros();
}

template<typename eT> const BandMat<eT>& BandMat<eT>::ones() {
	arma_extra_debug_sigprint();

	return fill(eT(1));
}

template<typename eT> const BandMat<eT>& BandMat<eT>::ones(const uword in_size, const uword& in_l, const uword& in_u) {
	arma_extra_debug_sigprint();

	set_size(in_size, in_l, in_u);

	return fill(eT(1));
}

template<typename eT> const BandMat<eT>& BandMat<eT>::eye() {
	arma_extra_debug_sigprint();

	(*this).zeros();

	for(uword ii = 0; ii < n_cols; ++ii) at(ii, ii) = eT(1);

	return *this;
}

template<typename eT> const BandMat<eT>& BandMat<eT>::eye(const uword in_size, const uword& in_l, const uword& in_u) {
	arma_extra_debug_sigprint();

	set_size(in_size, in_l, in_u);

	return (*this).eye();
}

template<typename eT> void BandMat<eT>::reset() {
	arma_extra_debug_sigprint();

	init_warm(0, 0, 0);
}
