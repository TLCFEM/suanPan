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
	unique_ptr<MetaMat<T>> factorize() override;

	unique_ptr<MetaMat<T>> i() override;
public:
	using MetaMat<T>::IPIV;
	using MetaMat<T>::factored;
	using MetaMat<T>::n_cols;
	using MetaMat<T>::n_rows;
	using MetaMat<T>::memory;

	FullMat();
	FullMat(uword, uword);

	unique_ptr<MetaMat<T>> make_copy() override;

	Mat<T> operator*(const Mat<T>&) override;

	int solve(Mat<T>&, const Mat<T>&) override;

	int solve_trs(Mat<T>&, const Mat<T>&) override;

	void save(const char*) override;
};

template<typename T> FullMat<T>::FullMat()
	: MetaMat<T>() {}

template<typename T> FullMat<T>::FullMat(const uword in_rows, const uword in_cols)
	: MetaMat<T>(in_rows, in_cols, in_rows * in_cols) {}

template<typename T> unique_ptr<MetaMat<T>> FullMat<T>::make_copy() { return make_unique<FullMat<T>>(*this); }

template<typename T> Mat<T> FullMat<T>::operator*(const Mat<T>& B) {
	Mat<T> C(size(B));

	if(1 == B.n_cols) {
		auto M = static_cast<int>(n_rows);
		auto N = static_cast<int>(n_cols);
		T ALPHA = 1.;
		auto LDA = M;
		auto INCX = 1;
		T BETA = 0.;
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
		auto M = static_cast<int>(n_rows);
		auto N = static_cast<int>(B.n_cols);
		auto K = static_cast<int>(n_cols);
		T ALPHA = 1.;
		auto LDA = M;
		auto LDB = K;
		T BETA = 0.;
		auto LDC = M;

		if(std::is_same<T, float>::value) {
			using E = float;
			arma_fortran(arma_sgemm)(&this->TRAN, &this->TRAN, &M, &N, &K, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)B.memptr(), &LDB, (E*)&BETA, (E*)C.memptr(), &LDC);
		}
		else if(std::is_same<T, double>::value) {
			using E = double;
			arma_fortran(arma_dgemm)(&this->TRAN, &this->TRAN, &M, &N, &K, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)B.memptr(), &LDB, (E*)&BETA, (E*)C.memptr(), &LDC);
		}
	}

	return C;
}

template<typename T> int FullMat<T>::solve(Mat<T>& X, const Mat<T>& B) {
	if(factored) {
		suanpan_debug("the matrix is factored.\n");
		return this->solve_trs(X, B);
	}

	auto INFO = 0;

	X = B;

	auto N = static_cast<int>(n_rows);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDA = N;
	auto LDB = static_cast<int>(B.n_rows);
	IPIV.zeros(N);

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_sgesv)(&N, &NRHS, (E*)this->memptr(), &LDA, IPIV.memptr(), (E*)X.memptr(), &LDB, &INFO);
	}
	else if(std::is_same<T, double>::value) {
		using E = double;
		arma_fortran(arma_dgesv)(&N, &NRHS, (E*)this->memptr(), &LDA, IPIV.memptr(), (E*)X.memptr(), &LDB, &INFO);
	}

	if(INFO == 0) factored = true;
	else suanpan_error("solve() receives error code %u from the base driver, the matrix is probably singular.\n", INFO);

	return INFO;

}

template<typename T> int FullMat<T>::solve_trs(Mat<T>& X, const Mat<T>& B) {
	if(!factored) {
		suanpan_debug("the matrix is not factored.\n");
		return this->solve(X, B);
	}

	if(IPIV.is_empty()) return -1;

	X = B;

	auto N = static_cast<int>(n_rows);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDA = N;
	auto LDB = static_cast<int>(B.n_rows);
	auto INFO = 0;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_sgetrs)(const_cast<char*>(&this->TRAN), &N, &NRHS, (E*)this->memptr(), &LDA, IPIV.memptr(), (E*)X.memptr(), &LDB, &INFO);
	}
	else if(std::is_same<T, double>::value) {
		using E = double;
		arma_fortran(arma_dgetrs)(const_cast<char*>(&this->TRAN), &N, &NRHS, (E*)this->memptr(), &LDA, IPIV.memptr(), (E*)X.memptr(), &LDB, &INFO);
	}

	return INFO;
}

template<typename T> unique_ptr<MetaMat<T>> FullMat<T>::factorize() {
	auto X = make_unique<FullMat<T>>(*this);

	if(factored) {
		suanpan_warning("the matrix is factored.\n");
		return X;
	}

	arma_debug_check(X->n_rows != X->n_cols, "i() only accepts sqaure matrix.");

	auto M = static_cast<int>(X->n_rows);
	auto N = M;
	auto LDA = M;
	X->IPIV.zeros(N);
	auto INFO = 0;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_sgetrf)(&M, &N, (E*)X->memptr(), &LDA, X->IPIV.memptr(), &INFO);
	}
	else if(std::is_same<T, double>::value) {
		using E = double;
		arma_fortran(arma_dgetrf)(&M, &N, (E*)X->memptr(), &LDA, X->IPIV.memptr(), &INFO);
	}

	if(INFO != 0) {
		suanpan_error("factorize() fails.\n");
		X->reset();
	}
	else X->factored = true;

	return X;
}

template<typename T> unique_ptr<MetaMat<T>> FullMat<T>::i() {
	auto X = make_unique<FullMat<T>>(*this);

	arma_debug_check(X->n_rows != X->n_cols, "i() only accepts sqaure matrix.");

	auto M = static_cast<int>(X->n_rows);
	auto N = M;
	auto LDA = M;
	X->IPIV.zeros(N);
	auto INFO = 0;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_sgetrf)(&M, &N, (E*)X->memptr(), &LDA, X->IPIV.memptr(), &INFO);
	}
	else if(std::is_same<T, double>::value) {
		using E = double;
		arma_fortran(arma_dgetrf)(&M, &N, (E*)X->memptr(), &LDA, X->IPIV.memptr(), &INFO);
	}

	if(INFO != 0) {
		X->reset();
		return X;
	}

	auto LWORK = 8 * M;
	const auto WORK = new T[LWORK];

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_sgetri)(&N, (E*)X->memptr(), &LDA, X->IPIV.memptr(), (E*)WORK, &LWORK, &INFO);
	}
	else if(std::is_same<T, double>::value) {
		using E = double;
		arma_fortran(arma_dgetri)(&N, (E*)X->memptr(), &LDA, X->IPIV.memptr(), (E*)WORK, &LWORK, &INFO);
	}

	delete[] WORK;

	if(INFO != 0) X->reset();

	return X;
}

template<typename T> void FullMat<T>::save(const char* name) {
	Mat<T> aux(access::rwp(this->memory), this->n_rows, this->n_cols, false, false);
	aux.save(name, raw_ascii);
}

#endif

//! @}
