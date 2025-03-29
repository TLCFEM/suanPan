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
 * @class SparseMatPARDISO
 * @brief A SparseMatPARDISO class that holds matrices.
 *
 * TODO: improve performance by storing and reusing factorization
 *
 * @author tlc
 * @date 21/03/2025
 * @version 0.1.0
 * @file SparseMatPARDISO.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef SPARSEMATPARDISO_HPP
#define SPARSEMATPARDISO_HPP

#ifdef SUANPAN_MKL

#include "../SparseMat.hpp"

#include <mkl_pardiso.h>

template<sp_d T> class SparseMatPARDISO final : public SparseMat<T> {
    const la_it maxfct = 1;
    const la_it mnum = 1;
    const la_it mtype = 11;
#ifdef SUANPAN_DEBUG
    const la_it msglvl = 1;
#else
    const la_it msglvl = 0;
#endif

    la_it iparm[64]{};
    std::int64_t pt[64]{};

    auto init_config() {
        pardisoinit(pt, &mtype, iparm);

        iparm[1] = 3;   // nested dissection algorithm
        iparm[23] = 10; // parallel factorization
        iparm[34] = 1;  // zero-based indexing
        if(std::is_same_v<T, float>) iparm[27] = 1;
    }

protected:
    using SparseMat<T>::direct_solve;

    int direct_solve(Mat<T>&, const Mat<T>&) override;

public:
    SparseMatPARDISO(const uword in_row, const uword in_col, const uword in_elem = 0)
        : SparseMat<T>(in_row, in_col, in_elem) { init_config(); }

    SparseMatPARDISO(const SparseMatPARDISO& other)
        : SparseMat<T>(other) { init_config(); }

    SparseMatPARDISO(SparseMatPARDISO&&) noexcept = delete;
    SparseMatPARDISO& operator=(const SparseMatPARDISO&) = delete;
    SparseMatPARDISO& operator=(SparseMatPARDISO&&) noexcept = delete;

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SparseMatPARDISO>(*this); }
};

template<sp_d T> int SparseMatPARDISO<T>::direct_solve(Mat<T>& X, const Mat<T>& B) {
    X.set_size(B.n_rows, B.n_cols);

    csr_form<T, la_it> csr_mat(this->triplet_mat, SparseBase::ZERO, true);

    const auto n = static_cast<la_it>(B.n_rows);
    const auto nrhs = static_cast<la_it>(B.n_cols);
    la_it info;

    la_it phase = 13;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, (void*)csr_mat.val_mem(), csr_mat.row_mem(), csr_mat.col_mem(), nullptr, &nrhs, iparm, &msglvl, (void*)B.memptr(), (void*)X.memptr(), &info);

    const auto error = info;
    if(0 != error) suanpan_error("Error code {} received.\n", error);

    phase = -1;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, nullptr, csr_mat.row_mem(), csr_mat.col_mem(), nullptr, &nrhs, iparm, &msglvl, nullptr, nullptr, &info);

    return 0 == error ? SUANPAN_SUCCESS : SUANPAN_FAIL;
}

#endif

#endif

//! @}
