template<typename eT> class SymmMat : public Base<eT, SymmMat<eT>> {
public:
	typedef eT elem_type;
	typedef typename get_pod_type<eT>::result pod_type;

	static const bool is_row = false;
	static const bool is_col = false;

	const uword n_size;
	const uword n_elem;
	const uhword mem_state;

	arma_aligned const eT*
	const mem;

	arma_align_mem eT mem_local[arma_config::mat_prealloc];

	~SymmMat();
	SymmMat();

	explicit SymmMat(const uword& in_size);
	explicit SymmMat(const SizeMat& s);

	template<typename fill_type> SymmMat(const uword& in_size, const fill::fill_class<fill_type>& f);
	template<typename fill_type> SymmMat(const SizeMat& s, const fill::fill_class<fill_type>& f);

	SymmMat(const Mat<eT>& in_mat);

	SymmMat& operator=(const eT& val);
	SymmMat& operator+=(const eT& val);
	SymmMat& operator-=(const eT& val);
	SymmMat& operator*=(const eT& val);
	SymmMat& operator/=(const eT& val);

	SymmMat(const SymmMat& m);
	SymmMat& operator=(const SymmMat& m);
	SymmMat& operator+=(const SymmMat& m);
	SymmMat& operator-=(const SymmMat& m);
	SymmMat& operator*=(const SymmMat& m) = delete;
	SymmMat& operator%=(const SymmMat& m);
	SymmMat& operator/=(const SymmMat& m);

	template<typename T1, typename op_type> SymmMat(const SmOp<T1, op_type>& X);

	template<typename T1, typename T2, typename glue_type> SymmMat& operator=(const Glue<T1, T2, glue_type>& X);
	template<typename T1, typename op_type> SymmMat& operator=(const SmOp<T1, op_type>& X);

	eT& at(const uword& in_row, const uword& in_col);
	const eT& at(const uword& in_row, const uword& in_col) const;
	eT& operator()(const uword& in_row, const uword& in_col);
	const eT& operator()(const uword& in_row, const uword& in_col) const;

	eT* memptr();
	const eT* memptr() const;

	void set_size(const uword in_size);

	arma_hot const SymmMat& fill(const eT val);
	template<typename fill_type> arma_hot const SymmMat& fill(const fill::fill_class<fill_type>&);

	const SymmMat& zeros();
	const SymmMat& zeros(const uword in_size);

	const SymmMat& ones();
	const SymmMat& ones(const uword in_size);

	const SymmMat& eye();
	const SymmMat& eye(const uword in_size);

	void reset();

	void init_cold();
	void init_warm(const uword& in_size);

	void print() const;
};
