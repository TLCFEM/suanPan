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
 * @class MetaMat
 * @brief A MetaMat class that holds matrices.
 *
 * @author tlc
 * @date 08/09/2017
 * @version 0.2.0
 * @file MetaMat.hpp
 * @addtogroup MetaMat
 * @{
 */

#ifndef METAMAT_HPP
#define METAMAT_HPP

#include <Toolbox/debug.h>
#include "triplet_form.hpp"

enum class Precision { MIXED, FULL };

template<typename T> class MetaMat {
protected:
	bool factored = false;

	double tolerance = 1E-13;

	Precision precision = Precision::FULL;

	unsigned refinement = 10;
public:
	triplet_form<T, uword> triplet_mat;

	const uword n_rows;
	const uword n_cols;
	const uword n_elem;

	MetaMat(uword, uword, uword);
	MetaMat(const MetaMat&) = default;
	MetaMat(MetaMat&&) noexcept = delete;
	MetaMat& operator=(const MetaMat&) = delete;
	MetaMat& operator=(MetaMat&&) noexcept = delete;
	virtual ~MetaMat() = default;

	void set_tolerance(double);
	void set_precision(Precision);
	void set_refinement(unsigned);

	[[nodiscard]] virtual bool is_empty() const = 0;
	virtual void zeros() = 0;

	virtual unique_ptr<MetaMat> make_copy() = 0;

	virtual void unify(uword) = 0;

	[[nodiscard]] virtual T max() const = 0;

	virtual const T& operator()(uword, uword) const = 0;
	virtual T& at(uword, uword) = 0;

	[[nodiscard]] virtual const T* memptr() const = 0;
	virtual T* memptr() = 0;

	virtual void operator+=(const shared_ptr<MetaMat>&) = 0;
	virtual void operator-=(const shared_ptr<MetaMat>&) = 0;

	virtual Mat<T> operator*(const Mat<T>&) = 0;

	virtual void operator*=(T) = 0;

	Mat<T> solve(const Mat<T>&);
	Mat<T> solve(const SpMat<T>&);
	Mat<T> solve(Mat<T>&&);
	Mat<T> solve(SpMat<T>&&);

	virtual int solve(Mat<T>&, const Mat<T>&) = 0;
	virtual int solve(Mat<T>&, const SpMat<T>&) = 0;
	virtual int solve(Mat<T>&, Mat<T>&&) = 0;
	virtual int solve(Mat<T>&, SpMat<T>&&) = 0;

	[[nodiscard]] virtual int sign_det() const = 0;

	void save(const char*);

	virtual void csc_condense();
	virtual void csr_condense();
};

template<typename T> Mat<T> to_mat(const MetaMat<T>& in_mat) {
	Mat<T> out_mat(in_mat.n_rows, in_mat.n_cols);
	for(uword I = 0; I < in_mat.n_rows; ++I) for(uword J = 0; J < in_mat.n_cols; ++J) out_mat(I, J) = in_mat(I, J);
	return out_mat;
}

template<typename T> Mat<T> to_mat(const shared_ptr<MetaMat<T>>& in_mat) { return to_mat(*in_mat); }

template<typename T> MetaMat<T>::MetaMat(const uword in_rows, const uword in_cols, const uword in_elem)
	: triplet_mat(in_rows, in_cols)
	, n_rows(in_rows)
	, n_cols(in_cols)
	, n_elem(in_elem) {}

template<typename T> void MetaMat<T>::set_tolerance(const double TOL) { tolerance = TOL; }

template<typename T> void MetaMat<T>::set_precision(const Precision P) { precision = P; }

template<typename T> void MetaMat<T>::set_refinement(const unsigned R) { refinement = R; }

template<typename T> Mat<T> MetaMat<T>::solve(const Mat<T>& B) {
	Mat<T> X;
	if(0 != this->solve(X, B)) X.reset();
	return X;
}

template<typename T> Mat<T> MetaMat<T>::solve(const SpMat<T>& B) {
	Mat<T> X;
	if(0 != this->solve(X, B)) X.reset();
	return X;
}

template<typename T> Mat<T> MetaMat<T>::solve(Mat<T>&& B) {
	Mat<T> X;
	if(0 != this->solve(X, std::forward<Mat<T>>(B))) X.reset();
	return X;
}

template<typename T> Mat<T> MetaMat<T>::solve(SpMat<T>&& B) {
	Mat<T> X;
	if(0 != this->solve(X, std::forward<SpMat<T>>(B))) X.reset();
	return X;
}

template<typename T> int MetaMat<T>::solve(Mat<T>& X, const SpMat<T>& B) { return this->solve(X, Mat<T>(B)); }

template<typename T> int MetaMat<T>::solve(Mat<T>& X, Mat<T>&& B) { return this->solve(X, B); }

template<typename T> void MetaMat<T>::save(const char* name) { if(!to_mat(*this).save(name)) suanpan_error("cannot save matrix to file.\n"); }

template<typename T> void MetaMat<T>::csc_condense() {}

template<typename T> void MetaMat<T>::csr_condense() {}

#endif

//! @}
