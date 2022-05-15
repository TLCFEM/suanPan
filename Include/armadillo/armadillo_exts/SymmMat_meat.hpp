template<typename eT> SymmMat<eT>::~SymmMat() {}

template<typename eT> SymmMat<eT>::SymmMat()
	: n_size(0)
	, n_elem(0)
	, mem_state(0)
	, mem() {}

template<typename eT> SymmMat<eT>::SymmMat(const uword& in_size)
	: n_size(in_size)
	, n_elem((in_size + 1) * in_size / 2)
	, mem_state(0)
	, mem() {
	arma_extra_debug_sigprint_this(this);
	init_cold();
}

template<typename eT> SymmMat<eT>::SymmMat(const SizeMat& s)
	: n_size(s.n_rows)
	, n_elem((s.n_rows + 1) * s.n_cols / 2)
	, mem_state(0)
	, mem() {
	arma_debug_check(s.n_rows != s.n_cols, "SymmMat() only accepts sqaure matrix.");

	init_cold();
}

template<typename eT> template<typename fill_type> SymmMat<eT>::SymmMat(const uword& in_size, const fill::fill_class<fill_type>& f)
	: n_size(in_size)
	, n_elem((in_size + 1) * in_size / 2)
	, mem_state(0)
	, mem() {
	arma_extra_debug_sigprint_this(this);

	init_cold();

	(*this).fill(f);
}

template<typename eT> template<typename fill_type> SymmMat<eT>::SymmMat(const SizeMat& s, const fill::fill_class<fill_type>& f)
	: n_size(s.n_rows)
	, n_elem((s.n_rows + 1) * s.n_cols / 2)
	, mem_state(0)
	, mem() {
	arma_extra_debug_sigprint_this(this);

	arma_debug_check(s.n_rows != s.n_cols, "SymmMat() only accepts sqaure matrix.");

	init_cold();

	(*this).fill(f);
}

template<typename eT> SymmMat<eT>::SymmMat(const Mat<eT>& in_mat)
	: n_size(in_mat.n_rows)
	, n_elem((in_mat.n_rows + 1) * in_mat.n_cols / 2)
	, mem_state(0)
	, mem() {
	arma_debug_check(in_mat.n_rows != in_mat.n_cols, "SymmMat() only accepts sqaure matrix.");

	init_cold();

	auto tmp_ptr = const_cast<eT*>(mem);
	for(auto j = 0; j < n_size; ++j) for(auto i = 0; i <= j; ++i) *tmp_ptr++ = in_mat(i, j);
}

template<typename eT> SymmMat<eT>& SymmMat<eT>::operator=(const eT& val) {
	init_warm(1);
	access::rw(mem[0]) = val;
	return *this;
}

template<typename eT> SymmMat<eT>& SymmMat<eT>::operator+=(const eT& val) {
	arma_extra_debug_sigprint();

	arrayops::inplace_plus(memptr(), val, n_elem);

	return *this;
}

template<typename eT> SymmMat<eT>& SymmMat<eT>::operator-=(const eT& val) {
	arma_extra_debug_sigprint();

	arrayops::inplace_minus(memptr(), val, n_elem);

	return *this;
}

template<typename eT> SymmMat<eT>& SymmMat<eT>::operator*=(const eT& val) {
	arma_extra_debug_sigprint();

	arrayops::inplace_mul(memptr(), val, n_elem);

	return *this;
}

template<typename eT> SymmMat<eT>& SymmMat<eT>::operator/=(const eT& val) {
	arma_extra_debug_sigprint();

	arrayops::inplace_div(memptr(), val, n_elem);

	return *this;
}

template<typename eT> SymmMat<eT>::SymmMat(const SymmMat& m)
	: n_size(m.n_size)
	, n_elem(m.n_elem)
	, mem_state(0)
	, mem() {
	arma_extra_debug_sigprint(arma_str::format("this = %x   in_mat = %x") % this % &m);

	init_cold();

	arrayops::copy(memptr(), m.mem, m.n_elem);
}

template<typename eT> SymmMat<eT>& SymmMat<eT>::operator=(const SymmMat& m) {
	arma_extra_debug_sigprint(arma_str::format("this = %x   in_mat = %x") % this % &m);

	if(this != &m) {
		init_warm(m.n_size);

		arrayops::copy(memptr(), m.mem, m.n_elem);
	}

	return *this;
}

template<typename eT> SymmMat<eT>& SymmMat<eT>::operator+=(const SymmMat& m) {
	arma_extra_debug_sigprint();

	arma_debug_assert_same_size(*this, m, "addition");

	arrayops::inplace_plus(memptr(), m.memptr(), n_elem);

	return *this;
}

template<typename eT> SymmMat<eT>& SymmMat<eT>::operator-=(const SymmMat& m) {
	arma_extra_debug_sigprint();

	arma_debug_assert_same_size(*this, m, "subtraction");

	arrayops::inplace_minus(memptr(), m.memptr(), n_elem);

	return *this;
}

template<typename eT> SymmMat<eT>& SymmMat<eT>::operator%=(const SymmMat& m) {
	arma_extra_debug_sigprint();

	arma_debug_assert_same_size(*this, m, "element-wise multiplication");

	arrayops::inplace_mul(memptr(), m.memptr(), n_elem);

	return *this;
}

template<typename eT> SymmMat<eT>& SymmMat<eT>::operator/=(const SymmMat& m) {
	arma_extra_debug_sigprint();

	arma_debug_assert_same_size(*this, m, "element-wise division");

	arrayops::inplace_div(memptr(), m.memptr(), n_elem);

	return *this;
}

template<typename eT> template<typename T1, typename smop_type> SymmMat<eT>::SymmMat(const SmOp<T1, smop_type>& X)
	: n_size(0)
	, n_elem(0)
	, mem_state(0)
	, mem() {
	arma_extra_debug_sigprint_this(this);

	arma_type_check((is_same_type<eT, typename T1::elem_type>::no));

	smop_type::apply(*this, X);
}

template<typename eT> template<typename T1, typename T2, typename glue_type> SymmMat<eT>& SymmMat<eT>::operator=(const Glue<T1, T2, glue_type>& X) {
	arma_extra_debug_sigprint();

	arma_type_check((is_same_type<eT, typename T1::elem_type>::no));
	arma_type_check((is_same_type<eT, typename T2::elem_type>::no));

	glue_type::apply(*this, X);

	return *this;
}

template<typename eT> template<typename T1, typename smop_type> SymmMat<eT>& SymmMat<eT>::operator=(const SmOp<T1, smop_type>& X) {
	arma_extra_debug_sigprint();

	arma_type_check((is_same_type<eT, typename T1::elem_type>::no));

	smop_type::apply(*this, X);

	return *this;
}

template<typename eT> eT& SymmMat<eT>::at(const uword& in_row, const uword& in_col) {
	const auto tmp_loc = in_col > in_row ? (in_col * in_col + in_col) / 2 + in_row : (in_row * in_row + in_row) / 2 + in_col;

	return access::rw(mem[tmp_loc]);
}

template<typename eT> const eT& SymmMat<eT>::at(const uword& in_row, const uword& in_col) const {
	const auto tmp_loc = in_col > in_row ? (in_col * in_col + in_col) / 2 + in_row : (in_row * in_row + in_row) / 2 + in_col;

	return mem[tmp_loc];
}

template<typename eT> eT& SymmMat<eT>::operator()(const uword& in_row, const uword& in_col) {
	arma_debug_check(in_row >= n_size || in_col >= n_size, "SymmMat::operator(): index out of bounds");
	return at(in_row, in_col);
}

template<typename eT> const eT& SymmMat<eT>::operator()(const uword& in_row, const uword& in_col) const {
	arma_debug_check(in_row >= n_size || in_col >= n_size, "SymmMat::operator(): index out of bounds");
	return at(in_row, in_col);
}

template<typename eT> void SymmMat<eT>::init_cold() {
	arma_extra_debug_sigprint(arma_str::format("n_size = %d") % n_size);

#if(defined(ARMA_USE_CXX11) || defined(ARMA_64BIT_WORD))
    auto error_message = "SymmMat::init(): requested size is too large";
#else
	const char* error_message = "SymmMat::init(): requested size is too large; suggest to compile in C++11 mode or enable ARMA_64BIT_WORD";
#endif

	arma_debug_check(n_size > ARMA_MAX_UHWORD ? n_elem > ARMA_MAX_UWORD : false, error_message);

	if(n_elem <= arma_config::mat_prealloc)
		if(n_elem == 0) access::rw(mem) = NULL;
		else {
			arma_extra_debug_print("SymmMat::init(): using local memory");
			access::rw(mem) = mem_local;
		}
	else {
		arma_extra_debug_print("SymmMat::init(): acquiring memory");
		access::rw(mem) = memory::acquire<eT>(n_elem);
	}
}

template<typename eT> void SymmMat<eT>::init_warm(const uword& in_size) {
	arma_extra_debug_sigprint(arma_str::format("in_n_size = %d") % in_size);

	if(n_size == in_size) return;

	auto err_state = false;
	char* err_msg = nullptr;

	const uhword t_mem_state = mem_state;

	arma_debug_set_error(err_state, err_msg, t_mem_state == 3, "SymmMat::init(): size is fixed and hence cannot be changed");

#if(defined(ARMA_USE_CXX11) || defined(ARMA_64BIT_WORD))
    auto error_message = "SymmMat::init(): requested size is too large";
#else
	const char* error_message = "SymmMat::init(): requested size is too large; suggest to compile in C++11 mode or enable ARMA_64BIT_WORD";
#endif

	arma_debug_set_error(err_state, err_msg, in_size > ARMA_MAX_UHWORD ? (in_size + 1) * in_size / 2 > ARMA_MAX_UWORD : false, error_message);

	arma_debug_check(err_state, err_msg);

	const auto old_n_elem = n_elem;
	const auto new_n_elem = (in_size + 1) * in_size / 2;

	if(old_n_elem == new_n_elem) { arma_extra_debug_print("SymmMat::init(): reusing memory"); }
	else {
		arma_debug_check(t_mem_state == 2, "SymmMat::init(): mismatch between size of auxiliary memory and requested size");

		if(new_n_elem < old_n_elem) {
			if(t_mem_state == 0 && new_n_elem <= arma_config::mat_prealloc) {
				if(old_n_elem > arma_config::mat_prealloc) {
					arma_extra_debug_print("SymmMat::init(): releasing memory");
					memory::release(access::rw(mem));
				}

				access::rw(mem) = new_n_elem == 0 ? NULL : mem_local;
			}
			else { arma_extra_debug_print("SymmMat::init(): reusing memory"); }
		}
		else {
			if(t_mem_state == 0 && old_n_elem > arma_config::mat_prealloc) {
				arma_extra_debug_print("SymmMat::init(): releasing memory");
				memory::release(access::rw(mem));
			}

			access::rw(mem) = new_n_elem <= arma_config::mat_prealloc ? mem_local : memory::acquire<eT>(new_n_elem);

			access::rw(mem_state) = 0;
		}

		access::rw(n_size) = in_size;
		access::rw(n_elem) = new_n_elem;
	}
}

template<typename eT> void SymmMat<eT>::print() const {
	auto& o = std::cout;

	const auto save_flags = o.flags();
	const auto save_precision = o.precision();

	o.unsetf(ios::scientific);
	o.setf(ios::fixed);
	o.precision(4);

	for(auto i = 0; i < n_size; i++) {
		for(auto j = 0; j < n_size; j++) {
			o.width(8);
			o << at(i, j) << " ";
		}
		o << endl;
	}

	o.flags(save_flags);
	o.precision(save_precision);
}

template<typename eT> eT* SymmMat<eT>::memptr() { return const_cast<eT*>(mem); }

template<typename eT> const eT* SymmMat<eT>::memptr() const { return mem; }

template<typename eT> void SymmMat<eT>::set_size(const uword in_size) {
	arma_extra_debug_sigprint();

	init_warm(in_size);
}

template<typename eT> const SymmMat<eT>& SymmMat<eT>::fill(const eT val) {
	arma_extra_debug_sigprint();

	arrayops::inplace_set(memptr(), val, n_elem);

	return *this;
}

template<typename eT> template<typename fill_type> const SymmMat<eT>& SymmMat<eT>::fill(const fill::fill_class<fill_type>&) {
	arma_extra_debug_sigprint();

	if(is_same_type<fill_type, fill::fill_zeros>::yes) (*this).zeros();
	else if(is_same_type<fill_type, fill::fill_ones>::yes) (*this).ones();
	else if(is_same_type<fill_type, fill::fill_eye>::yes) (*this).eye();

	return *this;
}

template<typename eT> const SymmMat<eT>& SymmMat<eT>::zeros() {
	arma_extra_debug_sigprint();

	arrayops::fill_zeros(memptr(), n_elem);

	return *this;
}

template<typename eT> const SymmMat<eT>& SymmMat<eT>::zeros(const uword in_size) {
	arma_extra_debug_sigprint();

	set_size(in_size);

	return (*this).zeros();
}

template<typename eT> const SymmMat<eT>& SymmMat<eT>::ones() {
	arma_extra_debug_sigprint();

	return fill(eT(1));
}

template<typename eT> const SymmMat<eT>& SymmMat<eT>::ones(const uword in_size) {
	arma_extra_debug_sigprint();

	set_size(in_size);

	return fill(eT(1));
}

template<typename eT> const SymmMat<eT>& SymmMat<eT>::eye() {
	arma_extra_debug_sigprint();

	(*this).zeros();

	for(uword ii = 0; ii < n_size; ++ii) at(ii, ii) = eT(1);

	return *this;
}

template<typename eT> const SymmMat<eT>& SymmMat<eT>::eye(const uword in_size) {
	arma_extra_debug_sigprint();

	set_size(in_size);

	return (*this).eye();
}

template<typename eT> void SymmMat<eT>::reset() {
	arma_extra_debug_sigprint();

	init_warm(0);
}

template<typename eT> int sp_solve(Col<eT>& X, SymmMat<eT>& A, const Col<eT>& B) {
	X = B;

	auto UPLO = 'U';
	auto N = static_cast<int>(A.n_size);
	auto NRHS = 1;
	const auto IPIV = new int[N];
	auto LDB = N;
	auto INFO = 0;

	if(is_float<eT>::value) {
		using T = float;
		arma_fortran(arma_sspsv)(&UPLO, &N, &NRHS, (T*)A.memptr(), IPIV, (T*)X.memptr(), &LDB, &INFO);
	}
	else if(is_double<eT>::value) {
		using T = double;
		arma_fortran(arma_dspsv)(&UPLO, &N, &NRHS, (T*)A.memptr(), IPIV, (T*)X.memptr(), &LDB, &INFO);
	}

	delete[] IPIV;

	return INFO;
}
