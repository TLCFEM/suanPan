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

#include "SparseMat.hpp"
#include <lis/lislib.h>

template<sp_d T> class SparseMatLis final : public SparseMat<T> {
protected:
    using SparseMat<T>::direct_solve;

    int direct_solve(Mat<T>&, const Mat<T>&) override;

public:
    SparseMatLis(const uword in_row, const uword in_col, const uword in_elem = 0)
        : SparseMat<T>(in_row, in_col, in_elem) {}

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SparseMatLis>(*this); }
};

template<sp_d T> int SparseMatLis<T>::direct_solve(Mat<T>& X, const Mat<T>& B) {
    X.set_size(B.n_rows, B.n_cols);

    csr_form<double, LIS_INT> csr_mat(this->triplet_mat, SparseBase::ZERO, true);

    const auto n = csr_mat.n_rows;
    const auto nnz = csr_mat.n_elem;

    LIS_MATRIX A;
    LIS_VECTOR b, x;
    LIS_SOLVER solver;

    lis_matrix_create(0, &A);
    lis_matrix_set_size(A, n, 0);
    lis_matrix_set_csr(nnz, csr_mat.row_mem(), csr_mat.col_mem(), csr_mat.val_mem(), A);
    lis_matrix_assemble(A);

    lis_solver_create(&solver);

    lis_vector_create(0, &b);
    lis_vector_create(0, &x);
    lis_vector_set_size(b, n, 0);
    lis_vector_set_size(x, n, 0);
    lis_vector_set(b, (double*)B.memptr());
    lis_vector_set(x, (double*)X.memptr());

    lis_solve(A, b, x, solver);

    lis_matrix_destroy(A);
    lis_vector_destroy(b);
    lis_vector_destroy(x);
    lis_solver_destroy(solver);

    return SUANPAN_SUCCESS;
}

#endif

//! @}
