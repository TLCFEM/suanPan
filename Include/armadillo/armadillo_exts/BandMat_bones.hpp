template<typename eT> class BandMat : public Base<eT, BandMat<eT>> {
public:
	typedef eT elem_type;
	typedef typename get_pod_type<eT>::result pod_type;

	static const bool is_row = false;
	static const bool is_col = false;

	const uword n_cols;
	const uword n_l;
	const uword n_u;
	const uword n_rows = 2 * n_l + n_u + 1;

private:
	const uword n_s = n_l + n_u;
	const uword n_a = n_rows - n_l;
	const uword n_b = n_rows - 1;

public:
	const uword n_elem;
	const uhword mem_state;

	arma_aligned const eT*
	const mem;

	arma_align_mem eT mem_local[arma_config::mat_prealloc];

	~BandMat();
	BandMat();

	explicit BandMat(const uword& in_size, const uword& in_l, const uword& in_u);

	template<typename fill_type> BandMat(const uword& in_size, const uword& in_l, const uword& in_u, const fill::fill_class<fill_type>& f);

	BandMat& operator=(const eT& val);
	BandMat& operator+=(const eT& val);
	BandMat& operator-=(const eT& val);
	BandMat& operator*=(const eT& val);
	BandMat& operator/=(const eT& val);

	BandMat(const BandMat& m);
	BandMat& operator=(const BandMat& m);
	BandMat& operator+=(const BandMat& m);
	BandMat& operator-=(const BandMat& m);
	BandMat& operator*=(const BandMat& m) = delete;
	BandMat& operator%=(const BandMat& m);
	BandMat& operator/=(const BandMat& m);

	template<typename T1, typename bdop_type> BandMat(const BdOp<T1, bdop_type>& X);

	template<typename T1, typename bdop_type> BandMat& operator=(const BdOp<T1, bdop_type>& X);

	eT& at(const uword& in_row, const uword& in_col);
	const eT& at(const uword& in_row, const uword& in_col) const;
	eT& operator()(const uword& in_row, const uword& in_col);
	const eT& operator()(const uword& in_row, const uword& in_col) const;

	eT* memptr();
	const eT* memptr() const;

	void set_size(const uword in_size, const uword& in_l, const uword& in_u);

	arma_hot const BandMat& fill(const eT val);
	template<typename fill_type> arma_hot const BandMat& fill(const fill::fill_class<fill_type>&);

	const BandMat& zeros();
	const BandMat& zeros(const uword in_size, const uword& in_l, const uword& in_u);

	const BandMat& ones();
	const BandMat& ones(const uword in_size, const uword& in_l, const uword& in_u);

	const BandMat& eye();
	const BandMat& eye(const uword in_size, const uword& in_l, const uword& in_u);

	void reset();

	void init_cold();
	void init_warm(const uword& in_size, const uword& in_l, const uword& in_u);
};
