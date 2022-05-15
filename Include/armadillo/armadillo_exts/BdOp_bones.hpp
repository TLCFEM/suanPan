template<typename T1, typename bdop_type> class BdOp : public Base<typename T1::elem_type, BdOp<T1, bdop_type>> {
public:
	typedef typename T1::elem_type elem_type;
	typedef typename get_pod_type<elem_type>::result pod_type;

	static const bool is_row = false;
	static const bool is_col = false;

	explicit BdOp(const T1& in_m);
	BdOp(const T1& in_m, const elem_type in_aux);
	BdOp(const T1& in_m, const uword in_aux_uword_a, const uword in_aux_uword_b);
	~BdOp();

	arma_aligned const T1& m;           //!< storage of reference to the operand (eg. a matrix)
	arma_aligned elem_type aux;         //!< storage of auxiliary data, user defined format
	arma_aligned uword aux_uword_a = 0; //!< storage of auxiliary data, uword format
	arma_aligned uword aux_uword_b = 0; //!< storage of auxiliary data, uword format
};
