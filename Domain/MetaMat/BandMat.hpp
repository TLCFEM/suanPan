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
	static const char TRAN;

	static T bin;

	const uword l_band;
	const uword u_band;
	const uword s_band;
	const uword m_rows;       // memory block layout
	podarray<float> s_memory; // float storage used in mixed precision algorithm
protected:
	int solve_trs(Mat<T>&, const Mat<T>&) override;
	int solve_trs(Mat<T>&, Mat<T>&&) override;
public:
	BandMat();
	BandMat(uword, uword, uword);

	unique_ptr<MetaMat<T>> make_copy() override;

	void unify(uword) override;

	const T& operator()(uword, uword) const override;
	T& at(uword, uword) override;

	Mat<T> operator*(const Mat<T>&) override;

	int solve(Mat<T>&, const Mat<T>&) override;
	int solve(Mat<T>&, Mat<T>&&) override;
};

template<typename T> const char BandMat<T>::TRAN = 'N';

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

template<typename T> void BandMat<T>::unify(const uword idx) {
#ifdef SUANPAN_MT
	tbb::parallel_for(std::max(idx, u_band) - u_band, std::min(this->n_rows, idx + l_band + 1), [&](const uword I) { access::rw(this->memory[I + s_band + idx * (m_rows - 1)]) = 0.; });
	tbb::parallel_for(std::max(idx, l_band) - l_band, std::min(this->n_cols, idx + u_band + 1), [&](const uword I) { access::rw(this->memory[idx + s_band + I * (m_rows - 1)]) = 0.; });
#else
	for(auto I = std::max(idx, u_band) - u_band; I < std::min(this->n_rows, idx + l_band + 1); ++I) access::rw(this->memory[I + s_band + idx * (m_rows - 1)]) = 0.;
	for(auto I = std::max(idx, l_band) - l_band; I < std::min(this->n_cols, idx + u_band + 1); ++I) access::rw(this->memory[idx + s_band + I * (m_rows - 1)]) = 0.;
#endif
	access::rw(this->memory[s_band + idx * m_rows]) = 1.;
}

template<typename T> const T& BandMat<T>::operator()(const uword in_row, const uword in_col) const {
	if(in_row > in_col + l_band || in_row + u_band < in_col) return bin = 0.;

	return this->memory[in_row + s_band + in_col * (m_rows - 1)];
}

template<typename T> T& BandMat<T>::at(const uword in_row, const uword in_col) {
	if(in_row > in_col + l_band || in_row + u_band < in_col) return bin = 0.;

	return access::rw(this->memory[in_row + s_band + in_col * (m_rows - 1)]);
}

template<typename T> Mat<T> BandMat<T>::operator*(const Mat<T>& X) {
	Mat<T> Y(size(X));

	auto M = static_cast<int>(this->n_rows);
	auto N = static_cast<int>(this->n_cols);
	auto KL = static_cast<int>(l_band);
	auto KU = static_cast<int>(u_band);
	T ALPHA = 1.;
	auto LDA = static_cast<int>(m_rows);
	auto INC = 1;
	T BETA = 0.;

	if(std::is_same<T, float>::value) {
		using E = float;
#ifdef SUANPAN_MT
		tbb::parallel_for(0llu, X.n_cols, [&](const uword I) { arma_fortran(arma_sgbmv)(&TRAN, &M, &N, &KL, &KU, (E*)&ALPHA, (E*)(this->memptr() + l_band), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
#else
		for(uword I = 0; I < X.n_cols; ++I)
			arma_fortran(arma_sgbmv)(&TRAN, &M, &N, &KL, &KU, (E*)&ALPHA, (E*)(this->memptr() + l_band), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC);
#endif
	}
	else if(std::is_same<T, double>::value) {
		using E = double;
#ifdef SUANPAN_MT
		tbb::parallel_for(0llu, X.n_cols, [&](const uword I) { arma_fortran(arma_dgbmv)(&TRAN, &M, &N, &KL, &KU, (E*)&ALPHA, (E*)(this->memptr() + l_band), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
#else
		for(uword I = 0; I < X.n_cols; ++I)
			arma_fortran(arma_dgbmv)(&TRAN, &M, &N, &KL, &KU, (E*)&ALPHA, (E*)(this->memptr() + l_band), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC);
#endif
	}

	return Y;
}

template<typename T> int BandMat<T>::solve(Mat<T>& X, const Mat<T>& B) {
	if(this->factored) return this->solve_trs(X, B);

	suanpan_debug([&]() { if(this->n_rows != this->n_cols) throw invalid_argument("requires a square matrix"); });

	auto INFO = 0;

	auto N = static_cast<int>(this->n_rows);
	auto KL = static_cast<int>(l_band);
	auto KU = static_cast<int>(u_band);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDAB = static_cast<int>(m_rows);
	auto LDB = static_cast<int>(B.n_rows);
	this->IPIV.zeros(N);

	this->factored = true;

	if(std::is_same<T, float>::value) {
		using E = float;

		X = B;
		arma_fortran(arma_sgbsv)(&N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, this->IPIV.memptr(), (E*)X.memptr(), &LDB, &INFO);
	}
	else if(Precision::FULL == this->precision) {
		using E = double;

		X = B;
		arma_fortran(arma_dgbsv)(&N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, this->IPIV.memptr(), (E*)X.memptr(), &LDB, &INFO);
	}
	else {
		s_memory = this->to_float();

		arma_fortran(arma_sgbtrf)(&N, &N, &KL, &KU, s_memory.memptr(), &LDAB, this->IPIV.memptr(), &INFO);

		if(0 == INFO) INFO = this->solve_trs(X, B);
	}

	if(0 != INFO) suanpan_error("solve() receives error code %u from base driver, the matrix is probably singular.\n", INFO);

	return INFO;
}

template<typename T> int BandMat<T>::solve_trs(Mat<T>& X, const Mat<T>& B) {
	auto INFO = 0;

	auto N = static_cast<int>(this->n_rows);
	auto KL = static_cast<int>(l_band);
	auto KU = static_cast<int>(u_band);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDAB = static_cast<int>(m_rows);
	auto LDB = static_cast<int>(B.n_rows);

	if(std::is_same<T, float>::value) {
		using E = float;

		X = B;
		arma_fortran(arma_sgbtrs)(&TRAN, &N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, this->IPIV.memptr(), (E*)X.memptr(), &LDB, &INFO);
	}
	else if(Precision::FULL == this->precision) {
		using E = double;

		X = B;
		arma_fortran(arma_dgbtrs)(&TRAN, &N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, this->IPIV.memptr(), (E*)X.memptr(), &LDB, &INFO);
	}
	else {
		X = arma::zeros(B.n_rows, B.n_cols);

		mat full_residual = B;

		auto multiplier = 1.;

		auto counter = 0;
		while(++counter < 20) {
			auto residual = conv_to<fmat>::from(full_residual / multiplier);

			arma_fortran(arma_sgbtrs)(&TRAN, &N, &KL, &KU, &NRHS, s_memory.memptr(), &LDAB, this->IPIV.memptr(), residual.memptr(), &LDB, &INFO);
			if(0 != INFO) break;

			const auto incre = multiplier * conv_to<mat>::from(residual);

			X += incre;

			multiplier = norm(full_residual -= this->operator*(incre));

			suanpan_debug("mixed precision algorithm multiplier: %.5E\n", multiplier);

			if(multiplier < this->tolerance) break;
		}
	}

	if(INFO != 0) suanpan_error("solve() receives error code %u from base driver, the matrix is probably singular.\n", INFO);

	return INFO;
}

template<typename T> int BandMat<T>::solve(Mat<T>& X, Mat<T>&& B) {
	if(this->factored) return this->solve_trs(X, std::forward<Mat<T>>(B));

	suanpan_debug([&]() { if(this->n_rows != this->n_cols) throw invalid_argument("requires a square matrix"); });

	auto INFO = 0;

	auto N = static_cast<int>(this->n_rows);
	auto KL = static_cast<int>(l_band);
	auto KU = static_cast<int>(u_band);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDAB = static_cast<int>(m_rows);
	auto LDB = static_cast<int>(B.n_rows);
	this->IPIV.zeros(N);

	this->factored = true;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_sgbsv)(&N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, this->IPIV.memptr(), (E*)B.memptr(), &LDB, &INFO);
		X = std::move(B);
	}
	else if(Precision::FULL == this->precision) {
		using E = double;
		arma_fortran(arma_dgbsv)(&N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, this->IPIV.memptr(), (E*)B.memptr(), &LDB, &INFO);
		X = std::move(B);
	}
	else {
		s_memory = this->to_float();

		arma_fortran(arma_sgbtrf)(&N, &N, &KL, &KU, s_memory.memptr(), &LDAB, this->IPIV.memptr(), &INFO);

		if(0 == INFO) INFO = this->solve_trs(X, std::forward<Mat<T>>(B));
	}

	if(0 != INFO) suanpan_error("solve() receives error code %u from base driver, the matrix is probably singular.\n", INFO);

	return INFO;
}

template<typename T> int BandMat<T>::solve_trs(Mat<T>& X, Mat<T>&& B) {
	auto INFO = 0;

	auto N = static_cast<int>(this->n_rows);
	auto KL = static_cast<int>(l_band);
	auto KU = static_cast<int>(u_band);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDAB = static_cast<int>(m_rows);
	auto LDB = static_cast<int>(B.n_rows);

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_sgbtrs)(&TRAN, &N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, this->IPIV.memptr(), (E*)B.memptr(), &LDB, &INFO);
		X = std::move(B);
	}
	else if(Precision::FULL == this->precision) {
		using E = double;
		arma_fortran(arma_dgbtrs)(&TRAN, &N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, this->IPIV.memptr(), (E*)B.memptr(), &LDB, &INFO);
		X = std::move(B);
	}
	else {
		X = arma::zeros(B.n_rows, B.n_cols);

		auto multiplier = 1.;

		auto counter = 0;
		while(++counter < 20) {
			auto residual = conv_to<fmat>::from(B / multiplier);

			arma_fortran(arma_sgbtrs)(&TRAN, &N, &KL, &KU, &NRHS, s_memory.memptr(), &LDAB, this->IPIV.memptr(), residual.memptr(), &LDB, &INFO);
			if(0 != INFO) break;

			const auto incre = multiplier * conv_to<mat>::from(residual);

			X += incre;

			multiplier = norm(B -= this->operator*(incre));

			suanpan_debug("mixed precision algorithm multiplier: %.5E\n", multiplier);

			if(multiplier < this->tolerance) break;
		}
	}

	if(INFO != 0) suanpan_error("solve() receives error code %u from base driver, the matrix is probably singular.\n", INFO);

	return INFO;
}
#endif

//! @}
