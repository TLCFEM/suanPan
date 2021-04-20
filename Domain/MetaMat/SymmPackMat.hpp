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
 * @class SymmPackMat
 * @brief A SymmPackMat class that holds matrices.
 *
 * @author tlc
 * @date 06/09/2017
 * @version 0.1.0
 * @file SymmPackMat.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef SYMMPACKMAT_HPP
#define SYMMPACKMAT_HPP

#include "DenseMat.hpp"

template<typename T> class SymmPackMat final : public DenseMat<T> {
	static constexpr char UPLO = 'U';

	static T bin;

	int solve_trs(Mat<T>&, Mat<T>&&);
	int solve_trs(Mat<T>&, const Mat<T>&);
public:
	explicit SymmPackMat(uword);

	unique_ptr<MetaMat<T>> make_copy() override;

	void unify(uword) override;

	const T& operator()(uword, uword) const override;
	T& at(uword, uword) override;

	Mat<T> operator*(const Mat<T>&) override;

	int solve(Mat<T>&, Mat<T>&&) override;
	int solve(Mat<T>&, const Mat<T>&) override;
};

template<typename T> T SymmPackMat<T>::bin = 0.;

template<typename T> SymmPackMat<T>::SymmPackMat(const uword in_size)
	: DenseMat<T>(in_size, in_size, (in_size + 1) * in_size / 2) {}

template<typename T> unique_ptr<MetaMat<T>> SymmPackMat<T>::make_copy() { return std::make_unique<SymmPackMat<T>>(*this); }

template<typename T> void SymmPackMat<T>::unify(const uword K) {
#ifdef SUANPAN_MT
	tbb::parallel_for(0llu, K, [&](const uword I) { access::rw(this->memory[(K * K + K) / 2 + I]) = 0.; });
	tbb::parallel_for(K + 1llu, this->n_rows, [&](const uword I) { access::rw(this->memory[(I * I + I) / 2 + K]) = 0.; });
#else
	for(auto I = 0llu; I < K; ++I) access::rw(this->memory[(K * K + K) / 2 + I]) = 0.;
	for(auto I = K + 1llu; I < this->n_rows; ++I) access::rw(this->memory[(I * I + I) / 2 + K]) = 0.;
#endif
	access::rw(this->memory[(K * K + 3 * K) / 2]) = 1.;
}

template<typename T> const T& SymmPackMat<T>::operator()(const uword in_row, const uword in_col) const { return this->memory[in_col > in_row ? (in_col * in_col + in_col) / 2 + in_row : (in_row * in_row + in_row) / 2 + in_col]; }

template<typename T> T& SymmPackMat<T>::at(const uword in_row, const uword in_col) { return in_col < in_row ? bin : access::rw(this->memory[(in_col * in_col + in_col) / 2 + in_row]); }

template<const char S, const char T, typename T1> Mat<T1> spmm(const SymmPackMat<T1>& A, const Mat<T1>& B);

template<typename T> Mat<T> SymmPackMat<T>::operator*(const Mat<T>& X) {
	if(!X.is_colvec()) return spmm<'R', 'N'>(*this, X);

	auto Y = X;

	auto N = static_cast<int>(this->n_rows);
	auto INC = 1;
	T ALPHA = 1.;
	T BETA = 0.;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_sspmv)(&UPLO, &N, (E*)&ALPHA, (E*)this->memptr(), (E*)X.memptr(), &INC, (E*)&BETA, (E*)Y.memptr(), &INC);
	}
	else if(std::is_same<T, double>::value) {
		using E = double;
		arma_fortran(arma_dspmv)(&UPLO, &N, (E*)&ALPHA, (E*)this->memptr(), (E*)X.memptr(), &INC, (E*)&BETA, (E*)Y.memptr(), &INC);
	}

	return Y;
}

template<typename T> int SymmPackMat<T>::solve(Mat<T>& X, const Mat<T>& B) {
	if(this->factored) return this->solve_trs(X, B);

	auto N = static_cast<int>(this->n_rows);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDB = static_cast<int>(B.n_rows);
	auto INFO = 0;

	this->factored = true;

	if(std::is_same<T, float>::value) {
		using E = float;
		X = B;
		arma_fortran(arma_sppsv)(&UPLO, &N, &NRHS, (E*)this->memptr(), (E*)X.memptr(), &LDB, &INFO);
	}
	else if(Precision::FULL == this->precision) {
		using E = double;
		X = B;
		arma_fortran(arma_dppsv)(&UPLO, &N, &NRHS, (E*)this->memptr(), (E*)X.memptr(), &LDB, &INFO);
	}
	else {
		this->s_memory = this->to_float();
		arma_fortran(arma_spptrf)(&UPLO, &N, this->s_memory.memptr(), &INFO);
		if(0 == INFO) INFO = this->solve_trs(X, B);
	}

	if(0 != INFO) suanpan_error("solve() receives error code %u from the base driver, the matrix is probably singular.\n", INFO);

	return INFO;
}

template<typename T> int SymmPackMat<T>::solve_trs(Mat<T>& X, const Mat<T>& B) {
	auto N = static_cast<int>(this->n_rows);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDB = static_cast<int>(B.n_rows);
	auto INFO = 0;

	if(std::is_same<T, float>::value) {
		using E = float;
		X = B;
		arma_fortran(arma_spptrs)(&UPLO, &N, &NRHS, (E*)this->memptr(), (E*)X.memptr(), &LDB, &INFO);
	}
	else if(Precision::FULL == this->precision) {
		using E = double;
		X = B;
		arma_fortran(arma_dpptrs)(&UPLO, &N, &NRHS, (E*)this->memptr(), (E*)X.memptr(), &LDB, &INFO);
	}
	else {
		X = arma::zeros(B.n_rows, B.n_cols);

		mat full_residual = B;

		auto multiplier = 1.;

		auto counter = 0;
		while(++counter < 20) {
			auto residual = conv_to<fmat>::from(full_residual / multiplier);

			arma_fortran(arma_spptrs)(&UPLO, &N, &NRHS, this->s_memory.memptr(), residual.memptr(), &LDB, &INFO);
			if(0 != INFO) break;

			const mat incre = multiplier * conv_to<mat>::from(residual);

			X += incre;

			suanpan_debug("mixed precision algorithm multiplier: %.5E\n", multiplier = norm(full_residual -= this->operator*(incre)));

			if(multiplier < this->tolerance) break;
		}
	}

	if(INFO != 0) suanpan_error("solve() receives error code %u from base driver, the matrix is probably singular.\n", INFO);

	return INFO;
}

template<typename T> int SymmPackMat<T>::solve(Mat<T>& X, Mat<T>&& B) {
	if(this->factored) return this->solve_trs(X, std::forward<Mat<T>>(B));

	auto N = static_cast<int>(this->n_rows);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDB = static_cast<int>(B.n_rows);
	auto INFO = 0;

	this->factored = true;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_sppsv)(&UPLO, &N, &NRHS, (E*)this->memptr(), (E*)B.memptr(), &LDB, &INFO);
		X = std::move(B);
	}
	else if(Precision::FULL == this->precision) {
		using E = double;
		arma_fortran(arma_dppsv)(&UPLO, &N, &NRHS, (E*)this->memptr(), (E*)B.memptr(), &LDB, &INFO);
		X = std::move(B);
	}
	else {
		this->s_memory = this->to_float();
		arma_fortran(arma_spptrf)(&UPLO, &N, this->s_memory.memptr(), &INFO);
		if(0 == INFO) INFO = this->solve_trs(X, std::forward<Mat<T>>(B));
	}

	if(0 != INFO) suanpan_error("solve() receives error code %u from the base driver, the matrix is probably singular.\n", INFO);

	return INFO;
}

template<typename T> int SymmPackMat<T>::solve_trs(Mat<T>& X, Mat<T>&& B) {
	auto N = static_cast<int>(this->n_rows);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDB = static_cast<int>(B.n_rows);
	auto INFO = 0;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_spptrs)(&UPLO, &N, &NRHS, (E*)this->memptr(), (E*)B.memptr(), &LDB, &INFO);
		X = std::move(B);
	}
	else if(Precision::FULL == this->precision) {
		using E = double;
		arma_fortran(arma_dpptrs)(&UPLO, &N, &NRHS, (E*)this->memptr(), (E*)B.memptr(), &LDB, &INFO);
		X = std::move(B);
	}
	else {
		X = arma::zeros(B.n_rows, B.n_cols);

		auto multiplier = 1.;

		auto counter = 0;
		while(++counter < 20) {
			auto residual = conv_to<fmat>::from(B / multiplier);

			arma_fortran(arma_spptrs)(&UPLO, &N, &NRHS, this->s_memory.memptr(), residual.memptr(), &LDB, &INFO);
			if(0 != INFO) break;

			const mat incre = multiplier * conv_to<mat>::from(residual);

			X += incre;

			suanpan_debug("mixed precision algorithm multiplier: %.5E\n", multiplier = norm(B -= this->operator*(incre)));

			if(multiplier < this->tolerance) break;
		}
	}

	if(INFO != 0) suanpan_error("solve() receives error code %u from base driver, the matrix is probably singular.\n", INFO);

	return INFO;
}

#endif

//! @}
