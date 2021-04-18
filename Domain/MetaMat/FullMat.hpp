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
 * @class FullMat
 * @brief A FullMat class that holds matrices.
 *
 * @author tlc
 * @date 06/09/2017
 * @version 0.1.0
 * @file FullMat.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef FULLMAT_HPP
#define FULLMAT_HPP

template<typename T> class FullMat : public MetaMat<T> {
	static const char TRAN;

	podarray<float> s_memory; // float storage used in mixed precision algorithm
protected:
	int solve_trs(Mat<T>&, const Mat<T>&) override;
public:
	FullMat();
	FullMat(uword, uword);

	unique_ptr<MetaMat<T>> make_copy() override;

	void unify(uword) override;

	const T& operator()(uword, uword) const override;
	T& at(uword, uword) override;

	Mat<T> operator*(const Mat<T>&) override;

	int solve(Mat<T>&, const Mat<T>&) override;

	void save(const char*) override;
};

template<typename T> const char FullMat<T>::TRAN = 'N';

template<typename T> FullMat<T>::FullMat()
	: MetaMat<T>() {}

template<typename T> FullMat<T>::FullMat(const uword in_rows, const uword in_cols)
	: MetaMat<T>(in_rows, in_cols, in_rows * in_cols) {}

template<typename T> unique_ptr<MetaMat<T>> FullMat<T>::make_copy() { return make_unique<FullMat<T>>(*this); }

template<typename T> void FullMat<T>::unify(const uword idx) {
#ifdef SUANPAN_MT
	tbb::parallel_for(0llu, this->n_rows, [&](const uword I) { at(I, idx) = 0.; });
	tbb::parallel_for(0llu, this->n_cols, [&](const uword I) { at(idx, I) = 0.; });
#else
	for(uword I = 0; I < this->n_rows; ++I) at(I, idx) = 0.;
	for(uword I = 0; I < this->n_cols; ++I) at(idx, I) = 0.;
#endif
	at(idx, idx) = 1.;
}

template<typename T> const T& FullMat<T>::operator()(const uword in_row, const uword in_col) const { return this->memory[in_row + in_col * this->n_rows]; }

template<typename T> T& FullMat<T>::at(const uword in_row, const uword in_col) { return access::rw(this->memory[in_row + in_col * this->n_rows]); }

template<typename T> Mat<T> FullMat<T>::operator*(const Mat<T>& B) {
	Mat<T> C(size(B));

	auto M = static_cast<int>(this->n_rows);
	auto N = static_cast<int>(this->n_cols);
	auto LDA = M;

	T ALPHA = 1.;
	T BETA = 0.;

	if(1 == B.n_cols) {
		auto INCX = 1;
		auto INCY = 1;

		if(std::is_same<T, float>::value) {
			using E = float;
			arma_fortran(arma_sgemv)(&TRAN, &M, &N, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)B.memptr(), &INCX, (E*)&BETA, (E*)C.memptr(), &INCY);
		}
		else if(std::is_same<T, double>::value) {
			using E = double;
			arma_fortran(arma_dgemv)(&TRAN, &M, &N, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)B.memptr(), &INCX, (E*)&BETA, (E*)C.memptr(), &INCY);
		}
	}
	else {
		auto K = static_cast<int>(B.n_cols);
		auto LDB = N;
		auto LDC = M;

		if(std::is_same<T, float>::value) {
			using E = float;
			arma_fortran(arma_sgemm)(&TRAN, &TRAN, &M, &K, &N, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)B.memptr(), &LDB, (E*)&BETA, (E*)C.memptr(), &LDC);
		}
		else if(std::is_same<T, double>::value) {
			using E = double;
			arma_fortran(arma_dgemm)(&TRAN, &TRAN, &M, &K, &N, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)B.memptr(), &LDB, (E*)&BETA, (E*)C.memptr(), &LDC);
		}
	}

	return C;
}

template<typename T> int FullMat<T>::solve(Mat<T>& X, const Mat<T>& B) {
	if(this->factored) return this->solve_trs(X, B);

	auto N = static_cast<int>(this->n_rows);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDB = static_cast<int>(B.n_rows);
	auto INFO = 0;

	this->IPIV.zeros(N);

	this->factored = true;

	if(std::is_same<T, float>::value) {
		using E = float;

		X = B;
		arma_fortran(arma_sgesv)(&N, &NRHS, (E*)this->memptr(), &N, this->IPIV.memptr(), (E*)X.memptr(), &LDB, &INFO);
	}
	else if(Precision::FULL == this->precision) {
		using E = double;

		X = B;
		arma_fortran(arma_dgesv)(&N, &NRHS, (E*)this->memptr(), &N, this->IPIV.memptr(), (E*)X.memptr(), &LDB, &INFO);
	}
	else {
		s_memory = this->to_float();

		arma_fortran(arma_sgetrf)(&N, &N, s_memory.memptr(), &N, this->IPIV.memptr(), &INFO);

		if(0 == INFO) INFO = solve_trs(X, B);
	}

	if(0 != INFO) suanpan_error("solve() receives error code %u from the base driver, the matrix is probably singular.\n", INFO);

	return INFO;
}

template<typename T> int FullMat<T>::solve_trs(Mat<T>& X, const Mat<T>& B) {
	auto N = static_cast<int>(this->n_rows);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDB = static_cast<int>(B.n_rows);
	auto INFO = 0;

	if(std::is_same<T, float>::value) {
		using E = float;

		X = B;
		arma_fortran(arma_sgetrs)(&TRAN, &N, &NRHS, (E*)this->memptr(), &N, this->IPIV.memptr(), (E*)X.memptr(), &LDB, &INFO);
	}
	else if(Precision::FULL == this->precision) {
		using E = double;

		X = B;
		arma_fortran(arma_dgetrs)(&TRAN, &N, &NRHS, (E*)this->memptr(), &N, this->IPIV.memptr(), (E*)X.memptr(), &LDB, &INFO);
	}
	else {
		X = arma::zeros(B.n_rows, B.n_cols);

		mat full_residual = B;

		auto multiplier = 1.;

		auto counter = 0;
		while(++counter < 20) {
			auto residual = conv_to<fmat>::from(full_residual / multiplier);

			arma_fortran(arma_sgetrs)(&TRAN, &N, &NRHS, s_memory.memptr(), &N, this->IPIV.memptr(), residual.memptr(), &LDB, &INFO);
			if(0 != INFO) break;

			multiplier = norm(full_residual = B - this->operator*(X += multiplier * conv_to<mat>::from(residual)));

			suanpan_debug("mixed precision algorithm multiplier: %.5E\n", multiplier);

			if(multiplier < this->tolerance) break;
		}
	}

	return INFO;
}

template<typename T> void FullMat<T>::save(const char* name) {
	Mat<T> aux(access::rwp(this->memory), this->n_rows, this->n_cols, false, false);
	aux.save(name, raw_ascii);
}

#endif

//! @}
