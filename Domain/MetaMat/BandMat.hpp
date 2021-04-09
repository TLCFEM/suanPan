/*******************************************************************************
 * Copyright (C) 2017-2021 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
/**
 * @class BandMat
 * @brief A BandMat class that holds matrices.
 *
 * @author tlc
 * @date 06/09/2017
 * @version 0.1.0
 * @file BandMat.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef BANDMAT_HPP
#define BANDMAT_HPP

template<typename T> class BandMat final : public MetaMat<T> {
	using MetaMat<T>::TRAN;

	static T bin;

	const uword l_band;
	const uword u_band;
	const uword s_band;
	const uword m_rows;       // memory block layout
	podarray<float> s_memory; // float storage used in mixed precision algorithm
protected:
	unique_ptr<MetaMat<T>> factorize() override;
public:
	using MetaMat<T>::IPIV;
	using MetaMat<T>::factored;
	using MetaMat<T>::n_cols;
	using MetaMat<T>::n_rows;
	using MetaMat<T>::n_elem;
	using MetaMat<T>::memory;
	using MetaMat<T>::precision;
	using MetaMat<T>::tolerance;

	BandMat();
	BandMat(uword, uword, uword);

	unique_ptr<MetaMat<T>> make_copy() override;

	const T& operator()(uword, uword) const override;
	T& at(uword, uword) override;

	Mat<T> operator*(const Mat<T>&) override;

	int solve(Mat<T>&, const Mat<T>&) override;
	int solve_trs(Mat<T>&, const Mat<T>&) override;
};

template<typename T> T BandMat<T>::bin = 0.;

template<typename T> BandMat<T>::BandMat()
	: MetaMat<T>()
	, l_band(0)
	, u_band(0)
	, s_band(0)
	, m_rows(0) {}

template<typename T> BandMat<T>::BandMat(const uword in_size, const uword in_l, const uword in_u)
	: MetaMat<T>(in_size, in_size, (2 * in_l + in_u + 1) * in_size)
	, l_band(in_l)
	, u_band(in_u)
	, s_band(in_l + in_u)
	, m_rows(2 * in_l + in_u + 1) {}

template<typename T> unique_ptr<MetaMat<T>> BandMat<T>::make_copy() { return make_unique<BandMat<T>>(*this); }

template<typename T> const T& BandMat<T>::operator()(const uword in_row, const uword in_col) const {
	// const auto n_bw = static_cast<long long>(in_row) - static_cast<long long>(in_col);
	// if(n_bw > static_cast<long long>(l_band) || n_bw < -static_cast<long long>(u_band)) return bin = 0.;

	if(in_row > in_col + l_band || in_row + u_band < in_col) return bin = 0.;

	return memory[in_row + s_band + in_col * (m_rows - 1)];
}

template<typename T> T& BandMat<T>::at(const uword in_row, const uword in_col) {
	// const auto n_bw = static_cast<long long>(in_row) - static_cast<long long>(in_col);
	// if(n_bw > static_cast<long long>(l_band) || n_bw < -static_cast<long long>(u_band)) return bin = 0.;

	if(in_row > in_col + l_band || in_row + u_band < in_col) return bin = 0.;

	return access::rw(memory[in_row + s_band + in_col * (m_rows - 1)]);
}

template<typename T> Mat<T> BandMat<T>::operator*(const Mat<T>& X) {
	Mat<T> Y(size(X));

	auto M = static_cast<int>(n_rows);
	auto N = static_cast<int>(n_cols);
	auto KL = static_cast<int>(l_band);
	auto KU = static_cast<int>(u_band);
	T ALPHA = 1.;
	auto LDA = static_cast<int>(m_rows);
	auto INC = 1;
	T BETA = 0.;

#ifdef SUANPAN_MT
	tbb::parallel_for(0llu, X.n_cols, [&](const uword I) {
		if(std::is_same<T, float>::value) {
			using E = float;
			arma_fortran(arma_sgbmv)(&TRAN, &M, &N, &KL, &KU, (E*)&ALPHA, (E*)(this->memptr() + l_band), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC);
		}
		else if(std::is_same<T, double>::value) {
			using E = double;
			arma_fortran(arma_dgbmv)(&TRAN, &M, &N, &KL, &KU, (E*)&ALPHA, (E*)(this->memptr() + l_band), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC);
		}
	});
#else
	for(uword I = 0; I < X.n_cols; ++I)
		if(std::is_same<T, float>::value) {
			using E = float;
			arma_fortran(arma_sgbmv)(&TRAN, &M, &N, &KL, &KU, (E*)&ALPHA, (E*)(this->memptr() + l_band), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC);
		}
		else if(std::is_same<T, double>::value) {
			using E = double;
			arma_fortran(arma_dgbmv)(&TRAN, &M, &N, &KL, &KU, (E*)&ALPHA, (E*)(this->memptr() + l_band), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC);
		}
#endif

	return Y;
}

template<typename T> int BandMat<T>::solve(Mat<T>& X, const Mat<T>& B) {
	auto INFO = 0;

	//
	// double precision
	//
	if(Precision::DOUBLE == precision) {
		if(factored) {
			suanpan_warning("the matrix is factored.\n");
			return this->solve_trs(X, B);
		}

		suanpan_debug([&]() { if(n_rows != n_cols) throw invalid_argument("requires a square matrix"); });

		X = B;

		auto N = static_cast<int>(n_rows);
		auto KL = static_cast<int>(l_band);
		auto KU = static_cast<int>(u_band);
		auto NRHS = static_cast<int>(B.n_cols);
		auto LDAB = static_cast<int>(m_rows);
		auto LDB = static_cast<int>(B.n_rows);
		IPIV.zeros(N);

		if(std::is_same<T, float>::value) {
			using E = float;
			arma_fortran(arma_sgbsv)(&N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, IPIV.memptr(), (E*)X.memptr(), &LDB, &INFO);
		}
		else if(std::is_same<T, double>::value) {
			using E = double;
			arma_fortran(arma_dgbsv)(&N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, IPIV.memptr(), (E*)X.memptr(), &LDB, &INFO);
		}

		if(INFO != 0) suanpan_error("solve() receives error code %u from base driver, the matrix is probably singular.\n", INFO);
		else factored = true;

		return INFO;
	}

	//
	// single precision (mixed precision)
	//

	s_memory.set_size(this->n_elem);

#ifdef SUANPAN_MT
	tbb::parallel_for(0llu, this->n_elem, [&](const uword I) { s_memory(I) = static_cast<float>(this->memory[I]); });
#else
	for(uword I = 0; I < this->n_elem; ++I) s_memory(I) = static_cast<float>(this->memory[I]);
#endif

	auto M = static_cast<int>(n_rows);
	auto N = static_cast<int>(n_cols);
	auto KL = static_cast<int>(l_band);
	auto KU = static_cast<int>(u_band);
	auto LDAB = static_cast<int>(m_rows);
	IPIV.zeros(N);

	arma_fortran(arma_sgbtrf)(&M, &N, &KL, &KU, s_memory.memptr(), &LDAB, this->IPIV.memptr(), &INFO);

	if(INFO != 0) return INFO;

	factored = true;

	X = arma::zeros(B.n_rows, B.n_cols);

	mat full_residual = B;
	fmat residual(B.n_rows, B.n_cols);
#ifdef SUANPAN_MT
	tbb::parallel_for(0llu, full_residual.n_elem, [&](const uword I) { residual(I) = static_cast<float>(full_residual(I)); });
#else
	for(uword I = 0; I < full_residual.n_elem; ++I) residual(I) = static_cast<float>(full_residual(I));
#endif

	auto multiplier = 1.;

	auto NRHS = static_cast<int>(B.n_cols);
	auto LDB = static_cast<int>(B.n_rows);
	auto ALPHA = -1.;
	auto BETA = 0.;
	auto INCXY = 1;

	auto counter = 0;
	while(++counter < 10) {
		arma_fortran(arma_sgbtrs)(&TRAN, &N, &KL, &KU, &NRHS, s_memory.memptr(), &LDAB, this->IPIV.memptr(), residual.memptr(), &LDB, &INFO);
		if(0 != INFO) break;

#ifdef SUANPAN_MT
		tbb::parallel_for(0llu, X.n_elem, [&](const uword I) { X(I) += multiplier * residual(I); });
#else
		for(uword I = 0; I < X.n_elem; ++I) X(I) += multiplier * residual(I);
#endif

		for(uword I = 0; I < B.n_cols; ++I)
			arma_fortran(arma_dgbmv)(&TRAN, &M, &N, &KL, &KU, &ALPHA, this->memptr() + l_band, &LDAB, X.colptr(I), &INCXY, &BETA, full_residual.colptr(I), &INCXY);

		full_residual += B;

		multiplier = norm(full_residual);

		suanpan_debug("mixed precision algorithm multiplier: %.5E\n", multiplier);

		if(multiplier < tolerance) break;

#ifdef SUANPAN_MT
		tbb::parallel_for(0llu, full_residual.n_elem, [&](const uword I) { residual(I) = static_cast<float>(full_residual(I) / multiplier); });
#else
		for(uword I = 0; I < full_residual.n_elem; ++I) residual(I) = static_cast<float>(full_residual(I) / multiplier);
#endif
	}

	return INFO;
}

template<typename T> int BandMat<T>::solve_trs(Mat<T>& X, const Mat<T>& B) {
	if(!factored) return this->solve(X, B);

	if(IPIV.is_empty()) return SUANPAN_FAIL;

	suanpan_debug([&]() { if(n_rows != n_cols) throw invalid_argument("requires a square matrix"); });

	auto INFO = 0;

	auto N = static_cast<int>(n_rows);
	auto KL = static_cast<int>(l_band);
	auto KU = static_cast<int>(u_band);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDAB = static_cast<int>(m_rows);
	auto LDB = static_cast<int>(B.n_rows);

	//
	// double precision
	//
	if(Precision::DOUBLE == precision) {
		X = B;

		if(std::is_same<T, float>::value) {
			using E = float;
			arma_fortran(arma_sgbtrs)(&TRAN, &N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, IPIV.memptr(), (E*)X.memptr(), &LDB, &INFO);
		}
		else if(std::is_same<T, double>::value) {
			using E = double;
			arma_fortran(arma_dgbtrs)(&TRAN, &N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, IPIV.memptr(), (E*)X.memptr(), &LDB, &INFO);
		}

		if(INFO != 0) suanpan_error("solve() receives error code %u from base driver, the matrix is probably singular.\n", INFO);

		return INFO;
	}

	//
	// single precision (mixed precision)
	//

	X = arma::zeros(B.n_rows, B.n_cols);

	mat full_residual = B;
	fmat residual(B.n_rows, B.n_cols);
#ifdef SUANPAN_MT
	tbb::parallel_for(0llu, full_residual.n_elem, [&](const uword I) { residual(I) = static_cast<float>(full_residual(I)); });
#else
	for(uword I = 0; I < full_residual.n_elem; ++I) residual(I) = static_cast<float>(full_residual(I));
#endif

	auto multiplier = 1.;

	auto ALPHA = -1.;
	auto BETA = 0.;
	auto INCXY = 1;

	auto counter = 0;
	while(++counter < 10) {
		arma_fortran(arma_sgbtrs)(&TRAN, &N, &KL, &KU, &NRHS, s_memory.memptr(), &LDAB, this->IPIV.memptr(), residual.memptr(), &LDB, &INFO);
		if(0 != INFO) break;

#ifdef SUANPAN_MT
		tbb::parallel_for(0llu, X.n_elem, [&](const uword I) { X(I) += multiplier * residual(I); });
#else
		for(uword I = 0; I < X.n_elem; ++I) X(I) += multiplier * residual(I);
#endif

		for(uword I = 0; I < B.n_cols; ++I)
			arma_fortran(arma_dgbmv)(&TRAN, &N, &N, &KL, &KU, &ALPHA, this->memptr() + l_band, &LDAB, X.colptr(I), &INCXY, &BETA, full_residual.colptr(I), &INCXY);

		full_residual += B;

		multiplier = norm(full_residual);

		suanpan_debug("mixed precision algorithm multiplier: %.5E\n", multiplier);

		if(multiplier < tolerance) break;

#ifdef SUANPAN_MT
		tbb::parallel_for(0llu, full_residual.n_elem, [&](const uword I) { residual(I) = static_cast<float>(full_residual(I) / multiplier); });
#else
		for(uword I = 0; I < full_residual.n_elem; ++I) residual(I) = static_cast<float>(full_residual(I) / multiplier);
#endif
	}

	return INFO;
}

template<typename T> unique_ptr<MetaMat<T>> BandMat<T>::factorize() {
	auto X = make_unique<BandMat<T>>(*this);

	if(factored) {
		suanpan_warning("the matrix is factored.\n");
		return X;
	}

	suanpan_debug([&]() { if(n_rows != n_cols) throw invalid_argument("requires a square matrix"); });

	auto M = static_cast<int>(n_rows);
	auto N = static_cast<int>(n_cols);
	auto KL = static_cast<int>(l_band);
	auto KU = static_cast<int>(u_band);
	auto LDAB = static_cast<int>(m_rows);
	X->IPIV.zeros(N);
	auto INFO = 0;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_sgbtrf)(&M, &N, &KL, &KU, (E*)X->memptr(), &LDAB, X->IPIV.memptr(), &INFO);
	}
	else if(std::is_same<T, double>::value) {
		using E = double;
		arma_fortran(arma_dgbtrf)(&M, &N, &KL, &KU, (E*)X->memptr(), &LDAB, X->IPIV.memptr(), &INFO);
	}

	if(INFO != 0) {
		suanpan_error("factorize() fails.\n");
		X->reset();
	}
	else X->factored = true;

	return X;
}

#endif

//! @}
