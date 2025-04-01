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
 * @class SparseMatClusterMUMPS
 * @brief A SparseMatClusterMUMPS class that holds matrices.
 *
 * @author tlc
 * @date 01/04/2025
 * @version 0.1.0
 * @file SparseMatClusterMUMPS.hpp
 * @addtogroup MetaMat
 * @{
 */

#ifndef SPARSEMATMPIMUMPS_HPP
#define SPARSEMATMPIMUMPS_HPP

#include "../SparseMat.hpp"

#include <ezp/ezp/mumps.hpp>

template<sp_d T, ezp::symmetric_pattern sym> class SparseMatBaseClusterMUMPS final : public SparseMat<T> {
    ezp::mumps<T, la_it> solver{sym, ezp::parallel_mode::no_host};

    std::vector<la_it> row_mem, col_mem;
    std::vector<T> val_mem;

    int solve_full(Mat<T>&);

protected:
    int direct_solve(Mat<T>& X, Mat<T>&& B) override { return this->solve_full(X = std::move(B)); }

    int direct_solve(Mat<T>& X, const Mat<T>& B) override { return this->solve_full(X = B); }

public:
    using SparseMat<T>::SparseMat;

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SparseMatBaseClusterMUMPS>(*this); }
};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnarrowing"
template<sp_d T, ezp::symmetric_pattern sym> int SparseMatBaseClusterMUMPS<T, sym>::solve_full(Mat<T>& X) {
    la_it info{-1};

    solver.icntl_printing_level(SUANPAN_VERBOSE ? 1 : 0);

    if(this->factored) info = solver.solve({X.n_rows, X.n_cols, X.memptr()});
    else {
        this->factored = true;

        const auto copy_idx = [&](const auto& coo_mat) {
            const auto plus_one = [](const auto v) { return v + 1; };

            row_mem.resize(coo_mat.n_elem);
            col_mem.resize(coo_mat.n_elem);
            std::transform(coo_mat.row_mem(), coo_mat.row_mem() + coo_mat.n_elem, row_mem.begin(), plus_one);
            std::transform(coo_mat.col_mem(), coo_mat.col_mem() + coo_mat.n_elem, col_mem.begin(), plus_one);
        };

        if(sym == ezp::symmetric_pattern::unsymmetric) {
            auto& coo_mat = this->triplet_mat;

            if(0 == comm_rank) copy_idx(coo_mat);

            info = solver.solve({coo_mat.n_rows, coo_mat.n_elem, row_mem.data(), col_mem.data(), coo_mat.val_mem()}, {X.n_rows, X.n_cols, X.memptr()});
        }
        else {
            // for symmetric matrices, MUMPS takes half the matrix
            auto coo_mat = this->triplet_mat.lower();
            coo_mat.csc_condense();

            if(0 == comm_rank) {
                val_mem = std::vector<T>(coo_mat.val_mem(), coo_mat.val_mem() + coo_mat.n_elem);

                copy_idx(coo_mat);
            }

            info = solver.solve({coo_mat.n_rows, coo_mat.n_elem, row_mem.data(), col_mem.data(), val_mem.data()}, {X.n_rows, X.n_cols, X.memptr()});
        }
    }

    if(0 == info) bcast_from_root(X);
    else suanpan_error("Error code {} received.\n", info);

    return info;
}
#pragma GCC diagnostic pop

template<sp_d T> using SparseMatClusterMUMPS = SparseMatBaseClusterMUMPS<T, ezp::symmetric_pattern::unsymmetric>;
template<sp_d T> using SparseSymmMatClusterMUMPS = SparseMatBaseClusterMUMPS<T, ezp::symmetric_pattern::symmetric_indefinite>;
template<sp_d T> using SparseSPDMatClusterMUMPS = SparseMatBaseClusterMUMPS<T, ezp::symmetric_pattern::symmetric_positive_definite>;

#endif

//! @}
