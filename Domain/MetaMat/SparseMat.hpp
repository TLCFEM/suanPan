/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
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
        this->triplet_mat.nullify(idx);
    }

    [[nodiscard]] T max() const override { return this->triplet_mat.max(); }

    T operator()(const uword in_row, const uword in_col) const override { return this->triplet_mat(in_row, in_col); }

    T& at(const uword in_row, const uword in_col) override {
        this->factored = false;
        return this->triplet_mat.at(in_row, in_col);
    }

    [[nodiscard]] const T* memptr() const override { throw std::invalid_argument("not supported"); }

    T* memptr() override { throw std::invalid_argument("not supported"); }

    void scale_accu(const T scalar, const shared_ptr<MetaMat<T>>& in_mat) override {
        if(nullptr == in_mat) return;
        if(!in_mat->triplet_mat.is_empty()) return this->scale_accu(scalar, in_mat->triplet_mat);
        this->factored = false;
        for(auto I = 0llu; I < in_mat->n_rows; ++I)
            for(auto J = 0llu; J < in_mat->n_cols; ++J)
                if(const auto t_val = in_mat->operator()(I, J); !suanpan::approx_equal(T(0), t_val)) at(I, J) = scalar * t_val;
    }

    void scale_accu(const T scalar, const triplet_form<T, uword>& in_mat) override {
        this->factored = false;
        if(1. == scalar) this->triplet_mat += in_mat;
        else if(-1. == scalar) this->triplet_mat -= in_mat;
        else this->triplet_mat.assemble(in_mat, 0, 0, scalar);
    }

    Mat<T> operator*(const Mat<T>& in_mat) const override { return this->triplet_mat * in_mat; }

    void operator*=(const T scalar) override {
        this->factored = false;
        this->triplet_mat *= scalar;
    }

    void allreduce() override {
#ifdef SUANPAN_DISTRIBUTED
        const auto coo_elem = this->triplet_mat.n_elem;
        std::vector<uword> dist_elem(comm_size);
        comm_world.allgather(coo_elem, dist_elem.data());
        this->triplet_mat.hack_size(std::accumulate(dist_elem.begin(), dist_elem.end(), uword{0}));

        mpl::irequest_pool requests;

        auto accu_elem = coo_elem;
        for(auto comm_n = 0; comm_n < comm_size; ++comm_n)
            if(comm_n == comm_rank) {
                requests.push(comm_world.ibcast(comm_n, this->triplet_mat.row_mem(), mpl::contiguous_layout<uword>{coo_elem}));
                requests.push(comm_world.ibcast(comm_n, this->triplet_mat.col_mem(), mpl::contiguous_layout<uword>{coo_elem}));
                requests.push(comm_world.ibcast(comm_n, this->triplet_mat.val_mem(), mpl::contiguous_layout<T>{coo_elem}));
            }
            else {
                requests.push(comm_world.ibcast(comm_n, this->triplet_mat.row_mem() + accu_elem, mpl::contiguous_layout<uword>{dist_elem[comm_n]}));
                requests.push(comm_world.ibcast(comm_n, this->triplet_mat.col_mem() + accu_elem, mpl::contiguous_layout<uword>{dist_elem[comm_n]}));
                requests.push(comm_world.ibcast(comm_n, this->triplet_mat.val_mem() + accu_elem, mpl::contiguous_layout<T>{dist_elem[comm_n]}));
                accu_elem += dist_elem[comm_n];
            }

        requests.waitall();
#endif
        csc_condense();
    }

    void csc_condense() override { this->triplet_mat.csc_condense(); }

    void csr_condense() override { this->triplet_mat.csr_condense(); }
};

#endif

//! @}
