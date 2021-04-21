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
 * @class SparseMat
 * @brief A SparseMat class that holds matrices.
 *
 * @author tlc
 * @date 06/05/2018
 * @version 0.1.0
 * @file SparseMat.hpp
 * @addtogroup MetaMat
 * @{
 */

#ifndef SPARSEMAT_HPP
#define SPARSEMAT_HPP

#include "MetaMat.hpp"

template<typename T> class SparseMat : public MetaMat<T> {
public:
	using MetaMat<T>::triplet_mat;
	using MetaMat<T>::solve;

	SparseMat(uword, uword, uword = 0);

	[[nodiscard]] bool is_empty() const override;
	void zeros() override;

	void unify(uword) override;

	T max() const override;

	const T& operator()(uword, uword) const override;
	T& at(uword, uword) override;

	const T* memptr() const override;
	T* memptr() override;

	void operator+=(const shared_ptr<MetaMat<T>>&) override;
	void operator-=(const shared_ptr<MetaMat<T>>&) override;

	Mat<T> operator*(const Mat<T>&) override;

	void operator*=(T) override;

	int solve(Mat<T>&, const SpMat<T>&) override;
	int solve(Mat<T>&, Mat<T>&&) override;
	int solve(Mat<T>&, SpMat<T>&&) override;

	[[nodiscard]] int sign_det() const override;

	void csc_condense() override;
	void csr_condense() override;
};

template<typename T> SparseMat<T>::SparseMat(const uword in_row, const uword in_col, const uword in_elem)
	: MetaMat<T>(in_row, in_col, 0) { triplet_mat.init(in_elem); }

template<typename T> bool SparseMat<T>::is_empty() const { return triplet_mat.is_empty(); }

template<typename T> void SparseMat<T>::zeros() {
	triplet_mat.zeros();
	this->factored = false;
}

template<typename T> void SparseMat<T>::unify(const uword idx) {
	using index_t = typename decltype(triplet_mat)::index_type;

	const auto t_idx = static_cast<index_t>(idx);
#ifdef SUANPAN_MT
	tbb::parallel_for(static_cast<index_t>(0), triplet_mat.n_elem, [&](const index_t I) { if(triplet_mat.row_idx[I] == t_idx || triplet_mat.col_idx[I] == t_idx) triplet_mat.val_idx[I] = 0.; });
#else
	for(index_t I = 0; I < triplet_mat.n_elem; ++I) if(triplet_mat.row_idx[I] == t_idx || triplet_mat.col_idx[I] == t_idx) triplet_mat.val_idx[I] = 0.;
#endif
	triplet_mat.at(t_idx, t_idx) = 1.;
}

template<typename T> T SparseMat<T>::max() const { return triplet_mat.max(); }

template<typename T> const T& SparseMat<T>::operator()(const uword in_row, const uword in_col) const {
	using index_t = typename decltype(triplet_mat)::index_type;
	return triplet_mat(static_cast<index_t>(in_row), static_cast<index_t>(in_col));
}

template<typename T> T& SparseMat<T>::at(const uword in_row, const uword in_col) {
	using index_t = typename decltype(triplet_mat)::index_type;
	return triplet_mat.at(static_cast<index_t>(in_row), static_cast<index_t>(in_col));
}

template<typename T> const T* SparseMat<T>::memptr() const { throw invalid_argument("not supproted"); }

template<typename T> T* SparseMat<T>::memptr() { throw invalid_argument("not supproted"); }

template<typename T> void SparseMat<T>::operator+=(const shared_ptr<MetaMat<T>>& in_mat) {
	triplet_mat += in_mat->triplet_mat;
	this->factored = false;
}

template<typename T> void SparseMat<T>::operator-=(const shared_ptr<MetaMat<T>>& in_mat) {
	triplet_mat -= in_mat->triplet_mat;
	this->factored = false;
}

template<typename T> Mat<T> SparseMat<T>::operator*(const Mat<T>& in_mat) { return triplet_mat * in_mat; }

template<typename T> void SparseMat<T>::operator*=(const T scalar) { triplet_mat *= scalar; }

template<typename T> int SparseMat<T>::solve(Mat<T>& X, const SpMat<T>& B) { return this->solve(X, Mat<T>(B)); }

template<typename T> int SparseMat<T>::solve(Mat<T>& X, Mat<T>&& B) { return this->solve(X, B); }

template<typename T> int SparseMat<T>::solve(Mat<T>& X, SpMat<T>&& B) { return this->solve(X, B); }

template<typename T> int SparseMat<T>::sign_det() const { throw invalid_argument("not supproted"); }

template<typename T> void SparseMat<T>::csc_condense() { triplet_mat.csc_condense(); }

template<typename T> void SparseMat<T>::csr_condense() { triplet_mat.csr_condense(); }

#endif

//! @}
