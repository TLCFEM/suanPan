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
 * @class SparseMatFGMRES
 * @brief A SparseMatFGMRES class that holds matrices.
 *
 * @author tlc
 * @date 27/03/2025
 * @version 0.2.0
 * @file SparseMatFGMRES.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef SPARSEMATFGMRES_HPP
#define SPARSEMATFGMRES_HPP

#include "../SparseMat.hpp"
#include "../csr_form.hpp"

#include <Toolbox/fgmres.hpp>

template<sp_d T> class SparseMatFGMRES final : public SparseMat<T> {
protected:
    using SparseMat<T>::direct_solve;

    int direct_solve(Mat<T>&, const Mat<T>&) override;

public:
    using SparseMat<T>::SparseMat;

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SparseMatFGMRES>(*this); }
};

template<sp_d T> int SparseMatFGMRES<T>::direct_solve(Mat<T>& X, const Mat<T>& B) {
    X.zeros(B.n_rows, B.n_cols);

    csr_form<T, la_it> csr_mat(this->triplet_mat);

    const auto precond = this->triplet_mat.diag();

    Col<int> info(B.n_cols);

    for(auto I = 0llu; I < B.n_cols; ++I) info[I] = fgmres_solve(csr_mat, precond, X.colptr(I), (double*)B.colptr(I), this->setting.tolerance);

    return info.min();
}

#endif

//! @}
