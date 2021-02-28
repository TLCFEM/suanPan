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

template<typename T> class BandSymmMat final : public MetaMat<T> {
	static const char UPLO;

	static T bin;

	const uword band;
	const uword m_rows; // memory block layout
protected:
	unique_ptr<MetaMat<T>> factorize() override;
public:
	using MetaMat<T>::factored;
	using MetaMat<T>::n_cols;
	using MetaMat<T>::n_rows;
	using MetaMat<T>::n_elem;
	using MetaMat<T>::memory;

	BandSymmMat();
	BandSymmMat(uword, uword);

	unique_ptr<MetaMat<T>> make_copy() override;

	void unify(uword) override;

	const T& operator()(uword, uword) const override;
	T& at(uword, uword) override;

	Mat<T> operator*(const Mat<T>&) override;

	int solve(Mat<T>&, const Mat<T>&) override;
	int solve_trs(Mat<T>&, const Mat<T>&) override;
};

template<typename T> const char BandSymmMat<T>::UPLO = 'L';

template<typename T> T BandSymmMat<T>::bin = 0.;

template<typename T> BandSymmMat<T>::BandSymmMat()
	: MetaMat<T>()
	, band(0)
	, m_rows(0) {}

template<typename T> BandSymmMat<T>::BandSymmMat(const uword in_size, const uword in_bandwidth)
	: MetaMat<T>(in_size, in_size, (in_bandwidth + 1llu) * in_size)
	, band(in_bandwidth)
	, m_rows(in_bandwidth + 1llu) {}

template<typename T> unique_ptr<MetaMat<T>> BandSymmMat<T>::make_copy() { return make_unique<BandSymmMat<T>>(*this); }

template<typename T> void BandSymmMat<T>::unify(const uword idx) {
#ifdef SUANPAN_MT
	tbb::parallel_for(0llu, idx, [&](const uword I) { at(idx, I) = 0.; });
	tbb::parallel_for(static_cast<uword>(idx + 1llu), n_rows, [&](const uword I) { at(I, idx) = 0.; });
#else
	for(uword I = 0; I < idx; ++I) at(idx, I) = 0.;
	for(auto I = idx + 1llu; I < n_rows; ++I) at(I, idx) = 0.;
#endif
	at(idx, idx) = 1.;
}

template<typename T> const T& BandSymmMat<T>::operator()(const uword in_row, const uword in_col) const {
	if(in_row > band + in_col) return bin = 0.;
	return memory[in_row > in_col ? in_row - in_col + in_col * m_rows : in_col - in_row + in_row * m_rows];
}

template<typename T> T& BandSymmMat<T>::at(const uword in_row, const uword in_col) {
	if(in_row > band + in_col || in_row < in_col) return bin = 0.;

	return access::rw(memory[in_row - in_col + in_col * m_rows]);
}

template<typename T> Mat<T> BandSymmMat<T>::operator*(const Mat<T>& X) {
	Mat<T> Y(size(X));

	auto N = static_cast<int>(n_cols);
	auto K = static_cast<int>(band);
	T ALPHA = 1.;
	auto LDA = static_cast<int>(m_rows);
	auto INC = 1;
	T BETA = 0.;

#ifdef SUANPAN_MT
	tbb::parallel_for(0llu, X.n_cols, [&](const uword I) {
		if(std::is_same<T, float>::value) {
			using E = float;
			arma_fortran(arma_ssbmv)(&UPLO, &N, &K, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC);
		}
		else if(std::is_same<T, double>::value) {
			using E = double;
			arma_fortran(arma_dsbmv)(&UPLO, &N, &K, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC);
		}
	});
#else
	for(uword I = 0; I < X.n_cols; ++I)
		if(std::is_same<T, float>::value) {
			using E = float;
			arma_fortran(arma_ssbmv)(&UPLO, &N, &K, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC);
		}
		else if(std::is_same<T, double>::value) {
			using E = double;
			arma_fortran(arma_dsbmv)(&UPLO, &N, &K, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC);
		}
#endif

	return Y;
}

template<typename T> int BandSymmMat<T>::solve(Mat<T>& X, const Mat<T>& B) {
	if(factored) {
		suanpan_debug("the matrix is factored.\n");
		return this->solve_trs(X, B);
	}

	X = B;

	auto N = static_cast<int>(n_rows);
	auto KD = static_cast<int>(band);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDAB = static_cast<int>(m_rows);
	auto LDB = static_cast<int>(B.n_rows);
	auto INFO = 0;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_spbsv)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)X.memptr(), &LDB, &INFO);
	}
	else if(std::is_same<T, double>::value) {
		using E = double;
		arma_fortran(arma_dpbsv)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)X.memptr(), &LDB, &INFO);
	}

	if(INFO != 0) suanpan_error("solve() receives error code %u from the base driver, the matrix is probably singular.\n", INFO);
	else factored = true;

	return INFO;
}

template<typename T> int BandSymmMat<T>::solve_trs(Mat<T>& X, const Mat<T>& B) {
	if(!factored) {
		suanpan_debug("the matrix is not factored.\n");
		return this->solve(X, B);
	}

	X = B;

	auto N = static_cast<int>(n_rows);
	auto KD = static_cast<int>(band);
	auto NRHS = static_cast<int>(B.n_cols);
	auto LDAB = static_cast<int>(m_rows);
	auto LDB = static_cast<int>(B.n_rows);
	auto INFO = 0;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_spbtrs)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)X.memptr(), &LDB, &INFO);
	}
	else if(std::is_same<T, double>::value) {
		using E = double;
		arma_fortran(arma_dpbtrs)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)X.memptr(), &LDB, &INFO);
	}

	if(INFO != 0) suanpan_error("solve() receives error code %u from the base driver, the matrix is probably singular.\n", INFO);

	return INFO;
}

template<typename T> unique_ptr<MetaMat<T>> BandSymmMat<T>::factorize() {
	auto X = make_unique<BandSymmMat<T>>(*this);

	if(factored) {
		suanpan_warning("the matrix is factored.\n");
		return X;
	}

	auto N = static_cast<int>(n_rows);
	auto KD = static_cast<int>(band);
	auto LDAB = static_cast<int>(m_rows);
	auto INFO = 0;

	if(std::is_same<T, float>::value) {
		using E = float;
		arma_fortran(arma_spbtrf)(&UPLO, &N, &KD, (E*)X->memptr(), &LDAB, &INFO);
	}
	else if(std::is_same<T, double>::value) {
		using E = double;
		arma_fortran(arma_dpbtrf)(&UPLO, &N, &KD, (E*)X->memptr(), &LDAB, &INFO);
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
