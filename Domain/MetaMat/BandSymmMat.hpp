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
 * @class BandSymmMat
 * @brief A BandSymmMat class that holds matrices.
 *
 * @author tlc
 * @date 06/09/2017
 * @version 0.1.0
 * @file BandSymmMat.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef BANDSYMMMAT_HPP
#define BANDSYMMMAT_HPP

#include "DenseMat.hpp"

template<typename T> class BandSymmMat final : public DenseMat<T> {
	static constexpr char UPLO = 'L';

	static T bin;

	const uword band;
	const uword m_rows; // memory block layout

	int solve_trs(Mat<T>&, Mat<T>&&);
	int solve_trs(Mat<T>&, const Mat<T>&);
public:
	BandSymmMat(uword, uword);

	unique_ptr<MetaMat<T>> make_copy() override;

	void unify(uword) override;

	const T& operator()(uword, uword) const override;
	T& at(uword, uword) override;

	Mat<T> operator*(const Mat<T>&) override;

	int solve(Mat<T>&, Mat<T>&&) override;
	int solve(Mat<T>&, const Mat<T>&) override;
};

template<typename T> T BandSymmMat<T>::bin = 0.;

template<typename T> BandSymmMat<T>::BandSymmMat(const uword in_size, const uword in_bandwidth)
	: DenseMat<T>(in_size, in_size, (in_bandwidth + 1) * in_size)
	, band(in_bandwidth)
	, m_rows(in_bandwidth + 1) {}

template<typename T> unique_ptr<MetaMat<T>> BandSymmMat<T>::make_copy() { return std::make_unique<BandSymmMat<T>>(*this); }

template<typename T> void BandSymmMat<T>::unify(const uword K) {
#ifdef SUANPAN_MT
	tbb::parallel_for(std::max(band, K) - band, K, [&](const uword I) { access::rw(this->memory[K - I + I * m_rows]) = 0.; });
	tbb::parallel_for(K + 1, std::min(this->n_rows, K + band + 1), [&](const uword I) { access::rw(this->memory[I - K + K * m_rows]) = 0.; });
#else
	for(auto I = std::max(band, K) - band; I < K; ++I) access::rw(this->memory[K - I + I * m_rows]) = 0.;
	for(auto I = K + 1; I < std::min(this->n_rows, K + band + 1); ++I) access::rw(this->memory[I - K + K * m_rows]) = 0.;
#endif
	access::rw(this->memory[K * m_rows]) = 1.;
}

template<typename T> const T& BandSymmMat<T>::operator()(const uword in_row, const uword in_col) const {
	if(in_row > band + in_col) return bin = 0.;
	return this->memory[in_row > in_col ? in_row - in_col + in_col * m_rows : in_col - in_row + in_row * m_rows];
}

template<typename T> T& BandSymmMat<T>::at(const uword in_row, const uword in_col) {
	if(in_row > band + in_col || in_row < in_col) return bin = 0.;
	return access::rw(this->memory[in_row - in_col + in_col * m_rows]);
}

template<typename T> Mat<T> BandSymmMat<T>::operator*(const Mat<T>& X) {
	Mat<T> Y(arma::size(X));

	auto N = static_cast<int>(this->n_cols);
	auto K = static_cast<int>(band);
	auto LDA = static_cast<int>(m_rows);
	auto INC = 1;
	T ALPHA = 1.;
	T BETA = 0.;

#ifdef SUANPAN_MT
	if(std::is_same<T, float>::value) {
		using E = float;
		tbb::parallel_for(0llu, X.n_cols, [&](const uword I) { arma_fortran(arma_ssbmv)(&UPLO, &N, &K, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
	}
	else if(std::is_same<T, double>::value) {
		using E = double;
		tbb::parallel_for(0llu, X.n_cols, [&](const uword I) { arma_fortran(arma_dsbmv)(&UPLO, &N, &K, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
	}
#else
	if(std::is_same<T, float>::value) {
		using E = float;
		for(uword I = 0; I < X.n_cols; ++I)
			arma_fortran(arma_ssbmv)(&UPLO, &N, &K, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC);
	}
	else if(std::is_same<T, double>::value) {
		using E = double;
		for(uword I = 0; I < X.n_cols; ++I)
			arma_fortran(arma_dsbmv)(&UPLO, &N, &K, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC);
	}
#endif

	return Y;
}

template<typename T> int BandSymmMat<T>::solve(Mat<T>& X, const Mat<T>& B) {
	if(this->factored) return this->solve_trs(X, B);

	auto N = static_cast<int>(this->n_rows);
	auto KD = static_cast<int>(band);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDAB = static_cast<int>(m_rows);
	auto LDB = static_cast<int>(B.n_rows);
	auto INFO = 0;

	this->factored = true;

	if(std::is_same<T, float>::value) {
		using E = float;
		X = B;
		arma_fortran(arma_spbsv)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)X.memptr(), &LDB, &INFO);
	}
	else if(Precision::FULL == this->precision) {
		using E = double;
		X = B;
		arma_fortran(arma_dpbsv)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)X.memptr(), &LDB, &INFO);
	}
	else {
		this->s_memory = this->to_float();
		arma_fortran(arma_spbtrf)(&UPLO, &N, &KD, this->s_memory.memptr(), &LDAB, &INFO);
		if(0 == INFO) INFO = this->solve_trs(X, B);
	}

	if(0 != INFO) suanpan_error("solve() receives error code %u from the base driver, the matrix is probably singular.\n", INFO);

	return INFO;
}

template<typename T> int BandSymmMat<T>::solve_trs(Mat<T>& X, const Mat<T>& B) {
	auto N = static_cast<int>(this->n_rows);
	auto KD = static_cast<int>(band);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDAB = static_cast<int>(m_rows);
	auto LDB = static_cast<int>(B.n_rows);
	auto INFO = 0;

	if(std::is_same<T, float>::value) {
		using E = float;
		X = B;
		arma_fortran(arma_spbtrs)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)X.memptr(), &LDB, &INFO);
	}
	else if(Precision::FULL == this->precision) {
		using E = double;
		X = B;
		arma_fortran(arma_dpbtrs)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)X.memptr(), &LDB, &INFO);
	}
	else {
		X = arma::zeros(B.n_rows, B.n_cols);

		mat full_residual = B;

		auto multiplier = norm(full_residual);

		auto counter = 0u;
		while(counter++ < this->refinement) {
			if(multiplier < this->tolerance) break;

			auto residual = conv_to<fmat>::from(full_residual / multiplier);

			arma_fortran(arma_spbtrs)(&UPLO, &N, &KD, &NRHS, this->s_memory.memptr(), &LDAB, residual.memptr(), &LDB, &INFO);
			if(0 != INFO) break;

			const mat incre = multiplier * conv_to<mat>::from(residual);

			X += incre;

			suanpan_debug("mixed precision algorithm multiplier: %.5E\n", multiplier = arma::norm(full_residual -= this->operator*(incre)));
		}
	}

	if(0 != INFO) suanpan_error("solve() receives error code %u from the base driver, the matrix is probably singular.\n", INFO);

	return INFO;
}

template<typename T> int BandSymmMat<T>::solve(Mat<T>& X, Mat<T>&& B) {
	if(this->factored) return this->solve_trs(X, std::forward<Mat<T>>(B));

	auto N = static_cast<int>(this->n_rows);
	auto KD = static_cast<int>(band);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDAB = static_cast<int>(m_rows);
	auto LDB = static_cast<int>(B.n_rows);
	auto INFO = 0;

	this->factored = true;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_spbsv)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)B.memptr(), &LDB, &INFO);
		X = std::move(B);
	}
	else if(Precision::FULL == this->precision) {
		using E = double;
		arma_fortran(arma_dpbsv)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)B.memptr(), &LDB, &INFO);
		X = std::move(B);
	}
	else {
		this->s_memory = this->to_float();
		arma_fortran(arma_spbtrf)(&UPLO, &N, &KD, this->s_memory.memptr(), &LDAB, &INFO);
		if(0 == INFO) INFO = this->solve_trs(X, std::forward<Mat<T>>(B));
	}

	if(0 != INFO) suanpan_error("solve() receives error code %u from the base driver, the matrix is probably singular.\n", INFO);

	return INFO;
}

template<typename T> int BandSymmMat<T>::solve_trs(Mat<T>& X, Mat<T>&& B) {
	auto N = static_cast<int>(this->n_rows);
	auto KD = static_cast<int>(band);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDAB = static_cast<int>(m_rows);
	auto LDB = static_cast<int>(B.n_rows);
	auto INFO = 0;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_spbtrs)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)B.memptr(), &LDB, &INFO);
		X = std::move(B);
	}
	else if(Precision::FULL == this->precision) {
		using E = double;
		arma_fortran(arma_dpbtrs)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)B.memptr(), &LDB, &INFO);
		X = std::move(B);
	}
	else {
		X = arma::zeros(B.n_rows, B.n_cols);

		auto multiplier = norm(B);

		auto counter = 0u;
		while(counter++ < this->refinement) {
			if(multiplier < this->tolerance) break;

			auto residual = conv_to<fmat>::from(B / multiplier);

			arma_fortran(arma_spbtrs)(&UPLO, &N, &KD, &NRHS, this->s_memory.memptr(), &LDAB, residual.memptr(), &LDB, &INFO);
			if(0 != INFO) break;

			const mat incre = multiplier * conv_to<mat>::from(residual);

			X += incre;

			suanpan_debug("mixed precision algorithm multiplier: %.5E\n", multiplier = arma::norm(B -= this->operator*(incre)));
		}
	}

	if(0 != INFO) suanpan_error("solve() receives error code %u from the base driver, the matrix is probably singular.\n", INFO);

	return INFO;
}

#endif

//! @}
