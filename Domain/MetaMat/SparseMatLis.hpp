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
 * @class SparseMatLis
 * @brief A SparseMatLis class that holds matrices.
 *
 * @author tlc
 * @date 15/06/2023
 * @version 0.1.0
 * @file SparseMatLis.hpp
 * @addtogroup MetaMat
 * @{
 */

#ifndef SPARSEMATLIS_HPP
#define SPARSEMATLIS_HPP

#include <lis/lislib.h>
#include "SparseMat.hpp"
#include "csr_form.hpp"

template<sp_d T> class SparseMatLis final : public SparseMat<T> {
    LIS_SOLVER solver = nullptr;

protected:
    using SparseMat<T>::direct_solve;
    using SparseMat<T>::setting;

    int direct_solve(Mat<T>&, const Mat<T>&) override;

public:
    SparseMatLis(const uword in_row, const uword in_col, const uword in_elem = 0)
        : SparseMat<T>(in_row, in_col, in_elem) {
        lis_solver_create(&solver);
        lis_solver_set_option("-i fgmres", solver);
    }

    SparseMatLis(const SparseMatLis& other)
        : SparseMat<T>(other) {
        lis_solver_create(&solver);
        lis_solver_set_option("-i fgmres", solver);
    }

    SparseMatLis(SparseMatLis&&) noexcept = delete;
    SparseMatLis& operator=(const SparseMatLis&) = delete;
    SparseMatLis& operator=(SparseMatLis&&) noexcept = delete;

    ~SparseMatLis() override { lis_solver_destroy(solver); }

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SparseMatLis>(*this); }
};

template<sp_d T> int SparseMatLis<T>::direct_solve(Mat<T>& X, const Mat<T>& B) {
    X.set_size(B.n_rows, B.n_cols);

    csr_form<double, LIS_INT> csr_mat(this->triplet_mat, SparseBase::ZERO, true);

    const sp_i auto n = csr_mat.n_rows;
    const sp_i auto nnz = csr_mat.n_elem;

    LIS_MATRIX A;
    LIS_VECTOR b, x;

    lis_matrix_create(0, &A);
    lis_matrix_set_size(A, n, 0);
    lis_matrix_set_csr(nnz, csr_mat.row_mem(), csr_mat.col_mem(), csr_mat.val_mem(), A);
    lis_matrix_assemble(A);

    lis_vector_create(0, &b);
    lis_vector_create(0, &x);
    lis_vector_set_size(b, n, 0);
    lis_vector_set_size(x, n, 0);

    lis_solver_set_option(setting.lis_options.c_str(), solver);

    for(uword I = 0; I < B.n_cols; ++I) {
        // ReSharper disable CppCStyleCast
        lis_vector_set(b, (double*)B.colptr(I));
        lis_vector_set(x, (double*)X.colptr(I));
        // ReSharper restore CppCStyleCast

        lis_solve(A, b, x, solver);
    }

    A->ptr = nullptr;
    A->index = nullptr;
    A->value = nullptr;
    b->value = nullptr;
    x->value = nullptr;

    lis_matrix_destroy(A);
    lis_vector_destroy(b);
    lis_vector_destroy(x);

    return SUANPAN_SUCCESS;
}

#endif

//! @}
