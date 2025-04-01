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
 * @class SparseMatClusterLIS
 * @brief A SparseMatClusterLIS class that holds matrices.
 *
 * @author tlc
 * @date 01/04/2025
 * @version 0.1.0
 * @file SparseMatClusterLIS.hpp
 * @addtogroup MetaMat
 * @{
 */

#ifndef SPARSEMATMPILIS_HPP
#define SPARSEMATMPILIS_HPP

#include "../SparseMat.hpp"

#include <ezp/ezp/lis.hpp>

template<sp_d T> class SparseMatClusterLIS final : public SparseMat<T> {
    ezp::lis solver{};

    csr_form<LIS_SCALAR, LIS_INT> csr_mat{};

    int solve_full(Mat<T>&);

protected:
    int direct_solve(Mat<T>& X, Mat<T>&& B) override { return this->solve_full(X = std::move(B)); }

    int direct_solve(Mat<T>& X, const Mat<T>& B) override { return this->solve_full(X = B); }

public:
    using SparseMat<T>::SparseMat;

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SparseMatClusterLIS>(*this); }
};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnarrowing"
template<sp_d T> int SparseMatClusterLIS<T>::solve_full(Mat<T>& X) {
    la_it info{-1};

    solver.set_option(this->setting.lis_options.c_str());

    if(this->factored) info = solver.solve({X.n_rows, X.n_cols, X.memptr()});
    else {
        this->factored = true;

        csr_mat = csr_form<LIS_SCALAR, LIS_INT>(this->triplet_mat, SparseBase::ZERO, false);

        info = solver.solve({csr_mat.n_rows, csr_mat.n_elem, csr_mat.row_mem(), csr_mat.col_mem(), csr_mat.val_mem()}, {X.n_rows, X.n_cols, X.memptr()});
    }

    if(0 == info) bcast_from_root(X);
    else suanpan_error("Error code {} received.\n", info);

    return info;
}
#pragma GCC diagnostic pop

#endif

//! @}
