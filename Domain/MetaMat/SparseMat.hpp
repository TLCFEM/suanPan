/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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
protected:
    using MetaMat<T>::direct_solve;

    int direct_solve(Mat<T>& X, Mat<T>&& B) override { return this->direct_solve(X, B); }

public:
    SparseMat(const uword in_row, const uword in_col, const uword in_elem = 0)
        : MetaMat<T>(in_row, in_col, 0) { this->triplet_mat.init(in_elem); }

    [[nodiscard]] bool is_empty() const override { return this->triplet_mat.is_empty(); }

    void zeros() override {
        this->factored = false;
        this->triplet_mat.zeros();
    }

    void nullify(const uword idx) override {
        this->factored = false;
        suanpan::for_each(static_cast<uword>(0), this->triplet_mat.n_elem, [&](const uword I) { if(this->triplet_mat.row(I) == idx || this->triplet_mat.col(I) == idx) this->triplet_mat.val_mem()[I] = T(0); });
    }

    [[nodiscard]] T max() const override { return this->triplet_mat.max(); }

    [[nodiscard]] Col<T> diag() const override { return this->triplet_mat.diag(); }

    T operator()(const uword in_row, const uword in_col) const override { return this->triplet_mat(in_row, in_col); }

    T& at(const uword in_row, const uword in_col) override {
        this->factored = false;
        return this->triplet_mat.at(in_row, in_col);
    }

    [[nodiscard]] const T* memptr() const override { throw invalid_argument("not supported"); }

    T* memptr() override { throw invalid_argument("not supported"); }

    void operator+=(const shared_ptr<MetaMat<T>>& in_mat) override {
        if(nullptr == in_mat) return;
        if(!in_mat->triplet_mat.is_empty()) return this->operator+=(in_mat->triplet_mat);
        this->factored = false;
        for(uword I = 0llu; I < in_mat->n_rows; ++I) for(uword J = 0llu; J < in_mat->n_cols; ++J) if(const auto t_val = in_mat->operator()(I, J); !suanpan::approx_equal(T(0), t_val)) at(I, J) = t_val;
    }

    void operator-=(const shared_ptr<MetaMat<T>>& in_mat) override {
        if(nullptr == in_mat) return;
        if(!in_mat->triplet_mat.is_empty()) return this->operator-=(in_mat->triplet_mat);
        this->factored = false;
        for(uword I = 0llu; I < in_mat->n_rows; ++I) for(uword J = 0llu; J < in_mat->n_cols; ++J) if(const auto t_val = in_mat->operator()(I, J); !suanpan::approx_equal(T(0), t_val)) at(I, J) = -t_val;
    }

    void operator+=(const triplet_form<T, uword>& in_mat) override {
        this->factored = false;
        this->triplet_mat += in_mat;
    }

    void operator-=(const triplet_form<T, uword>& in_mat) override {
        this->factored = false;
        this->triplet_mat -= in_mat;
    }

    Mat<T> operator*(const Mat<T>& in_mat) const override { return this->triplet_mat * in_mat; }

    void operator*=(const T scalar) override {
        this->factored = false;
        this->triplet_mat *= scalar;
    }

    [[nodiscard]] int sign_det() const override { throw invalid_argument("not supported"); }

    void csc_condense() override { this->triplet_mat.csc_condense(); }

    void csr_condense() override { this->triplet_mat.csr_condense(); }
};

#endif

//! @}
