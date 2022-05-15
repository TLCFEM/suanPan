/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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

template<sp_d T> class SparseMat : public MetaMat<T> {
public:
    using MetaMat<T>::triplet_mat;
    using MetaMat<T>::solve;

    SparseMat(uword, uword, uword = 0);

    [[nodiscard]] bool is_empty() const override;
    void zeros() override;

    void unify(uword) override;
    void nullify(uword) override;

    [[nodiscard]] T max() const override;

    const T& operator()(uword, uword) const override;
    T& at(uword, uword) override;

    [[nodiscard]] const T* memptr() const override;
    T* memptr() override;

    void operator+=(const shared_ptr<MetaMat<T>>&) override;
    void operator-=(const shared_ptr<MetaMat<T>>&) override;

    void operator+=(const triplet_form<T, uword>&) override;
    void operator-=(const triplet_form<T, uword>&) override;

    Mat<T> operator*(const Mat<T>&) override;

    void operator*=(T) override;

    [[nodiscard]] int sign_det() const override;

    void csc_condense() override;
    void csr_condense() override;
};

template<sp_d T> SparseMat<T>::SparseMat(const uword in_row, const uword in_col, const uword in_elem)
    : MetaMat<T>(in_row, in_col, 0) { triplet_mat.init(in_elem); }

template<sp_d T> bool SparseMat<T>::is_empty() const { return triplet_mat.is_empty(); }

template<sp_d T> void SparseMat<T>::zeros() {
    triplet_mat.zeros();
    this->factored = false;
}

template<sp_d T> void SparseMat<T>::unify(const uword idx) {
    nullify(idx);
    triplet_mat.at(idx, idx) = 1.;
}

template<sp_d T> void SparseMat<T>::nullify(const uword idx) {
    using index_t = typename decltype(triplet_mat)::index_type;

    const auto t_idx = static_cast<index_t>(idx);

    suanpan_for(static_cast<index_t>(0), triplet_mat.n_elem, [&](const index_t I) { if(triplet_mat.row(I) == t_idx || triplet_mat.col(I) == t_idx) triplet_mat.val_mem()[I] = 0.; });

    this->factored = false;
}

template<sp_d T> T SparseMat<T>::max() const { return triplet_mat.max(); }

template<sp_d T> const T& SparseMat<T>::operator()(const uword in_row, const uword in_col) const {
    using index_t = typename decltype(triplet_mat)::index_type;
    return triplet_mat(static_cast<index_t>(in_row), static_cast<index_t>(in_col));
}

template<sp_d T> T& SparseMat<T>::at(const uword in_row, const uword in_col) {
    this->factored = false;
    using index_t = typename decltype(triplet_mat)::index_type;
    return triplet_mat.at(static_cast<index_t>(in_row), static_cast<index_t>(in_col));
}

template<sp_d T> const T* SparseMat<T>::memptr() const { throw invalid_argument("not supported"); }

template<sp_d T> T* SparseMat<T>::memptr() { throw invalid_argument("not supported"); }

template<sp_d T> void SparseMat<T>::operator+=(const shared_ptr<MetaMat<T>>& in_mat) {
    if(nullptr == in_mat) return;

    if(!in_mat->triplet_mat.is_empty()) return this->operator+=(in_mat->triplet_mat);

    for(uword I = 0llu; I < in_mat->n_rows; ++I) for(uword J = 0llu; J < in_mat->n_cols; ++J) if(const auto t_val = in_mat->operator()(I, J); !suanpan::approx_equal(0., t_val)) at(I, J) = t_val;
}

template<sp_d T> void SparseMat<T>::operator-=(const shared_ptr<MetaMat<T>>& in_mat) {
    if(nullptr == in_mat) return;

    if(!in_mat->triplet_mat.is_empty()) return this->operator-=(in_mat->triplet_mat);

    for(uword I = 0llu; I < in_mat->n_rows; ++I) for(uword J = 0llu; J < in_mat->n_cols; ++J) if(const auto t_val = in_mat->operator()(I, J); !suanpan::approx_equal(0., t_val)) at(I, J) = -t_val;
}

template<sp_d T> void SparseMat<T>::operator+=(const triplet_form<T, uword>& in_mat) {
    this->triplet_mat += in_mat;
    this->factored = false;
}

template<sp_d T> void SparseMat<T>::operator-=(const triplet_form<T, uword>& in_mat) {
    this->triplet_mat -= in_mat;
    this->factored = false;
}

template<sp_d T> Mat<T> SparseMat<T>::operator*(const Mat<T>& in_mat) { return triplet_mat * in_mat; }

template<sp_d T> void SparseMat<T>::operator*=(const T scalar) { triplet_mat *= scalar; }

template<sp_d T> int SparseMat<T>::sign_det() const { throw invalid_argument("not supported"); }

template<sp_d T> void SparseMat<T>::csc_condense() { triplet_mat.csc_condense(); }

template<sp_d T> void SparseMat<T>::csr_condense() { triplet_mat.csr_condense(); }

#endif

//! @}
