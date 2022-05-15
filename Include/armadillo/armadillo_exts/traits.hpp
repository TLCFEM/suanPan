template<typename T> struct is_SymmMat {
	static const bool value = false;
};

template<typename eT> struct is_SymmMat<SymmMat<eT>> {
	static const bool value = true;
};

// template <typename T>
// struct is_BandMat {
//    static const bool value = false;
//};
//
// template <typename eT>
// struct is_BandMat<BandMat<eT>> {
//    static const bool value = true;
//};
