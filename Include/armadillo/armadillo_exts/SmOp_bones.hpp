template<typename T1, typename smop_type> class SmOp : public Base<typename T1::elem_type, SmOp<T1, smop_type>> {
public:
	typedef typename T1::elem_type elem_type;
	typedef typename get_pod_type<elem_type>::result pod_type;

	static const bool is_row = false;
	static const bool is_col = false;

	explicit SmOp(const T1& in_m);
	SmOp(const T1& in_m, const elem_type in_aux);
	SmOp(const T1& in_m, const uword in_aux_uword_a, const uword in_aux_uword_b);
	~SmOp();

	arma_aligned const T1& m;           //!< storage of reference to the operand (eg. a matrix)
	arma_aligned elem_type aux;         //!< storage of auxiliary data, user defined format
	arma_aligned uword aux_uword_a = 0; //!< storage of auxiliary data, uword format
	arma_aligned uword aux_uword_b = 0; //!< storage of auxiliary data, uword format
};
