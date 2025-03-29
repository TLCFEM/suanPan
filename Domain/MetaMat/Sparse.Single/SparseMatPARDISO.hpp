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
    static constexpr la_it negone{-1}, PARDISO_ANA_FACT{12}, PARDISO_SOLVE{33}, PARDISO_RELEASE{-1};

    const la_it maxfct{1}, mnum{1}, mtype{11};
#ifdef SUANPAN_DEBUG
    const la_it msglvl{1};
#else
    const la_it msglvl{0};
#endif

    la_it iparm[64]{};
    std::int64_t pt[64]{};

    csr_form<T, la_it> csr_mat{};

    bool is_allocated{false};

    auto init_config() {
        if constexpr(sizeof(la_it) == 4) pardisoinit(pt, &mtype, iparm);

        if constexpr(std::is_same_v<T, float>) iparm[27] = 1;
    }

    auto alloc() {
        dealloc();
        is_allocated = true;

        csr_mat = csr_form<T, la_it>(this->triplet_mat, SparseBase::ONE, true);

        la_it info{-1};
        if constexpr(sizeof(la_it) == 8) {
            using E = long long;
            pardiso_64(pt, (E*)&maxfct, (E*)&mnum, (E*)&mtype, (E*)&PARDISO_ANA_FACT, (E*)&csr_mat.n_rows, csr_mat.val_mem(), csr_mat.row_mem(), csr_mat.col_mem(), nullptr, (E*)&negone, iparm, (E*)&msglvl, nullptr, nullptr, (E*)&info);
        }
        else if constexpr(sizeof(la_it) == 4) {
            using E = int;
            pardiso(pt, (E*)&maxfct, (E*)&mnum, (E*)&mtype, (E*)&PARDISO_ANA_FACT, (E*)&csr_mat.n_rows, csr_mat.val_mem(), csr_mat.row_mem(), csr_mat.col_mem(), nullptr, (E*)&negone, iparm, (E*)&msglvl, nullptr, nullptr, (E*)&info);
        }
        return info;
    }

    auto dealloc() {
        if(!is_allocated) return;
        is_allocated = false;

        la_it info{-1};
        if constexpr(sizeof(la_it) == 8) {
            using E = long long;
            pardiso_64(pt, (E*)&maxfct, (E*)&mnum, (E*)&mtype, (E*)&PARDISO_RELEASE, (E*)&negone, nullptr, nullptr, nullptr, nullptr, (E*)&negone, (E*)iparm, (E*)&msglvl, nullptr, nullptr, (E*)&info);
        }
        else if constexpr(sizeof(la_it) == 4) {
            using E = int;
            pardiso(pt, (E*)&maxfct, (E*)&mnum, (E*)&mtype, (E*)&PARDISO_RELEASE, (E*)&negone, nullptr, nullptr, nullptr, nullptr, (E*)&negone, (E*)iparm, (E*)&msglvl, nullptr, nullptr, (E*)&info);
        }

        for(auto& i : pt) i = 0;
    }

protected:
    using SparseMat<T>::direct_solve;

    int direct_solve(Mat<T>&, const Mat<T>&) override;

public:
    SparseMatPARDISO(const uword in_row, const uword in_col, const uword in_elem = 0)
        : SparseMat<T>(in_row, in_col, in_elem) { init_config(); }

    SparseMatPARDISO(const SparseMatPARDISO& other)
        : SparseMat<T>(other) {
        init_config();
        this->factored = false;
    }

    SparseMatPARDISO(SparseMatPARDISO&&) noexcept = delete;
    SparseMatPARDISO& operator=(const SparseMatPARDISO&) = delete;
    SparseMatPARDISO& operator=(SparseMatPARDISO&&) noexcept = delete;

    ~SparseMatPARDISO() override { dealloc(); }

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SparseMatPARDISO>(*this); }
};

template<sp_d T> int SparseMatPARDISO<T>::direct_solve(Mat<T>& X, const Mat<T>& B) {
    if(!this->factored) {
        if(const auto info = alloc(); 0 != info) {
            suanpan_error("Error code {} received.\n", info);
            return SUANPAN_FAIL;
        }
        this->factored = true;
    }

    X.set_size(B.n_rows, B.n_cols);

    const la_it nrhs{static_cast<la_it>(B.n_cols)};

    la_it info{-1};
    if constexpr(sizeof(la_it) == 8) {
        using E = long long;
        pardiso_64(pt, (E*)&maxfct, (E*)&mnum, (E*)&mtype, (E*)&PARDISO_SOLVE, (E*)&csr_mat.n_rows, csr_mat.val_mem(), csr_mat.row_mem(), csr_mat.col_mem(), nullptr, (E*)&nrhs, iparm, (E*)&msglvl, (void*)B.memptr(), X.memptr(), (E*)&info);
    }
    else if constexpr(sizeof(la_it) == 4) {
        using E = int;
        pardiso(pt, (E*)&maxfct, (E*)&mnum, (E*)&mtype, (E*)&PARDISO_SOLVE, (E*)&csr_mat.n_rows, csr_mat.val_mem(), csr_mat.row_mem(), csr_mat.col_mem(), nullptr, (E*)&nrhs, iparm, (E*)&msglvl, (void*)B.memptr(), X.memptr(), (E*)&info);
    }

    if(0 != info) {
        suanpan_error("Error code {} received.\n", info);
        return SUANPAN_FAIL;
    }

    return SUANPAN_SUCCESS;
}

#endif

#endif

//! @}
