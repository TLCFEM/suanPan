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
protected:
	int solve_trs(Mat<T>&, const Mat<T>&) override;
public:
	FullMat();
	FullMat(uword, uword);

	unique_ptr<MetaMat<T>> make_copy() override;

	Mat<T> operator*(const Mat<T>&) override;

	int solve(Mat<T>&, const Mat<T>&) override;

	void save(const char*) override;
};

template<typename T> FullMat<T>::FullMat()
	: MetaMat<T>() {}

template<typename T> FullMat<T>::FullMat(const uword in_rows, const uword in_cols)
	: MetaMat<T>(in_rows, in_cols, in_rows * in_cols) {}

template<typename T> unique_ptr<MetaMat<T>> FullMat<T>::make_copy() { return make_unique<FullMat<T>>(*this); }

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
			arma_fortran(arma_sgemv)(&this->TRAN, &M, &N, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)B.memptr(), &INCX, (E*)&BETA, (E*)C.memptr(), &INCY);
		}
		else if(std::is_same<T, double>::value) {
			using E = double;
			arma_fortran(arma_dgemv)(&this->TRAN, &M, &N, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)B.memptr(), &INCX, (E*)&BETA, (E*)C.memptr(), &INCY);
		}
	}
	else {
		auto K = static_cast<int>(B.n_cols);
		auto LDB = N;
		auto LDC = M;

		if(std::is_same<T, float>::value) {
			using E = float;
			arma_fortran(arma_sgemm)(&this->TRAN, &this->TRAN, &M, &K, &N, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)B.memptr(), &LDB, (E*)&BETA, (E*)C.memptr(), &LDC);
		}
		else if(std::is_same<T, double>::value) {
			using E = double;
			arma_fortran(arma_dgemm)(&this->TRAN, &this->TRAN, &M, &K, &N, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)B.memptr(), &LDB, (E*)&BETA, (E*)C.memptr(), &LDC);
		}
	}

	return C;
}

template<typename T> int FullMat<T>::solve(Mat<T>& X, const Mat<T>& B) {
	if(this->factored) return this->solve_trs(X, B);

	auto N = static_cast<int>(this->n_rows);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDA = N;
	auto LDB = static_cast<int>(B.n_rows);
	auto INFO = 0;

	this->IPIV.zeros(N);

	X = B;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_sgesv)(&N, &NRHS, (E*)this->memptr(), &LDA, this->IPIV.memptr(), (E*)X.memptr(), &LDB, &INFO);
	}
	else if(std::is_same<T, double>::value) {
		using E = double;
		arma_fortran(arma_dgesv)(&N, &NRHS, (E*)this->memptr(), &LDA, this->IPIV.memptr(), (E*)X.memptr(), &LDB, &INFO);
	}

	if(0 == INFO) this->factored = true;
	else suanpan_error("solve() receives error code %u from the base driver, the matrix is probably singular.\n", INFO);

	return INFO;
}

template<typename T> int FullMat<T>::solve_trs(Mat<T>& X, const Mat<T>& B) {
	X = B;

	auto N = static_cast<int>(this->n_rows);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDA = N;
	auto LDB = static_cast<int>(B.n_rows);
	auto INFO = 0;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_sgetrs)(const_cast<char*>(&this->TRAN), &N, &NRHS, (E*)this->memptr(), &LDA, this->IPIV.memptr(), (E*)X.memptr(), &LDB, &INFO);
	}
	else if(std::is_same<T, double>::value) {
		using E = double;
		arma_fortran(arma_dgetrs)(const_cast<char*>(&this->TRAN), &N, &NRHS, (E*)this->memptr(), &LDA, this->IPIV.memptr(), (E*)X.memptr(), &LDB, &INFO);
	}

	return INFO;
}

template<typename T> void FullMat<T>::save(const char* name) {
	Mat<T> aux(access::rwp(this->memory), this->n_rows, this->n_cols, false, false);
	aux.save(name, raw_ascii);
}

#endif

//! @}
