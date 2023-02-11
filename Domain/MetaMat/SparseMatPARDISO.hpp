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
 * @class SparseMatPARDISO
 * @brief A SparseMatPARDISO class that holds matrices.
 *
 * TODO: improve performance by storing factorization and reusing it
 *
 * @author tlc
 * @date 20/01/2021
 * @version 0.1.0
 * @file SparseMatPARDISO.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef SPARSEMATPARDISO_HPP
#define SPARSEMATPARDISO_HPP

#ifdef SUANPAN_MKL

#include <mkl_pardiso.h>
#include "SparseMat.hpp"

template<sp_d T> class SparseMatPARDISO final : public SparseMat<T> {
    const int maxfct = 1;
    const int mnum = 1;
    const int mtype = 11;
#ifdef SUANPAN_DEBUG
    const int msglvl = 1;
#else
    const int msglvl = 0;
#endif

    int iparm[64];
    std::int64_t pt[64];

    int direct_solve(Mat<T>&, const Mat<T>&) override;

public:
    SparseMatPARDISO(uword, uword, uword = 0);

    unique_ptr<MetaMat<T>> make_copy() override;
};

template<sp_d T> SparseMatPARDISO<T>::SparseMatPARDISO(const uword in_row, const uword in_col, const uword in_elem)
    : SparseMat<T>(in_row, in_col, in_elem)
    , iparm{}
    , pt{} {
    pardisoinit(pt, &mtype, iparm);

    iparm[1] = 3;   // nested dissection algorithm
    iparm[23] = 10; // parallel factorization
    iparm[34] = 1;  // zero-based indexing
    if(std::is_same_v<T, float>) iparm[27] = 1;
}

template<sp_d T> unique_ptr<MetaMat<T>> SparseMatPARDISO<T>::make_copy() { return std::make_unique<SparseMatPARDISO>(*this); }

template<sp_d T> int SparseMatPARDISO<T>::direct_solve(Mat<T>& X, const Mat<T>& B) {
    X.set_size(B.n_rows, B.n_cols);

    csr_form<T, int> csr_mat(this->triplet_mat, SparseBase::ZERO, true);

    const auto n = static_cast<int>(B.n_rows);
    const auto nrhs = static_cast<int>(B.n_cols);
    int error;

    auto phase = 12;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, (void*)csr_mat.val_mem(), csr_mat.row_mem(), csr_mat.col_mem(), nullptr, &nrhs, iparm, &msglvl, nullptr, nullptr, &error);

    if(0 != error) {
        suanpan_error("Error code {} received.\n", error);
        return SUANPAN_FAIL;
    }

    phase = 33;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, (void*)csr_mat.val_mem(), csr_mat.row_mem(), csr_mat.col_mem(), nullptr, &nrhs, iparm, &msglvl, (void*)B.memptr(), (void*)X.memptr(), &error);

    if(0 != error) {
        suanpan_error("Error code {} received.\n", error);
        return SUANPAN_FAIL;
    }

    phase = -1;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, nullptr, csr_mat.row_mem(), csr_mat.col_mem(), nullptr, &nrhs, iparm, &msglvl, nullptr, nullptr, &error);

    return 0 == error ? SUANPAN_SUCCESS : SUANPAN_FAIL;
}

#endif

#endif

//! @}
