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

#include "../SparseMat.hpp"

extern "C" {
void pardisoinit(void* pt, const la_it* mtype, la_it* iparm);
void pardiso(void* pt, const la_it* maxfct, const la_it* mnum, const la_it* mtype, const la_it* phase, const la_it* n, const void* a, const la_it* ia, const la_it* ja, la_it* perm, const la_it* nrhs, la_it* iparm, const la_it* msglvl, void* b, void* x, la_it* error);
}

template<sp_d T, la_it MT> class SparseMatBasePARDISO final : public SparseMat<T> {
    static constexpr la_it negone{-1}, PARDISO_ANA_FACT{12}, PARDISO_SOLVE{33}, PARDISO_RELEASE{-1};

    const la_it maxfct{1}, mnum{1}, mtype{MT}, msglvl{SUANPAN_VERBOSE ? 1 : 0};

    la_it iparm[64]{};
    std::int64_t pt[64]{};

    csr_form<T, la_it> csr_mat{};

    bool is_allocated{false};

    auto init_config() {
        pardisoinit(pt, &mtype, iparm);

        if constexpr(std::is_same_v<T, float>) iparm[27] = 1;
    }

    auto alloc() {
        dealloc();
        is_allocated = true;

        // ReSharper disable CppDFAConstantConditions
        // ReSharper disable once CppDFAUnreachableCode
        csr_mat = csr_form<T, la_it>(1 == mtype || 11 == mtype ? this->triplet_mat : this->triplet_mat.upper(), SparseBase::ONE, true);
        // ReSharper restore CppDFAConstantConditions

        la_it info{-1};
        pardiso(pt, &maxfct, &mnum, &mtype, &PARDISO_ANA_FACT, &csr_mat.n_rows, csr_mat.val_mem(), csr_mat.row_mem(), csr_mat.col_mem(), nullptr, &negone, iparm, &msglvl, nullptr, nullptr, &info);

        return info;
    }

    auto dealloc() {
        if(!is_allocated) return;
        is_allocated = false;

        la_it info{-1};
        pardiso(pt, &maxfct, &mnum, &mtype, &PARDISO_RELEASE, &negone, nullptr, nullptr, nullptr, nullptr, &negone, iparm, &msglvl, nullptr, nullptr, &info);

        for(auto& i : pt) i = 0;
    }

protected:
    using SparseMat<T>::direct_solve;

    int direct_solve(Mat<T>&, const Mat<T>&) override;

public:
    SparseMatBasePARDISO(const uword in_row, const uword in_col, const uword in_elem = 0)
        : SparseMat<T>(in_row, in_col, in_elem) { init_config(); }

    SparseMatBasePARDISO(const SparseMatBasePARDISO& other)
        : SparseMat<T>(other) {
        init_config();
        this->factored = false;
    }

    SparseMatBasePARDISO(SparseMatBasePARDISO&&) noexcept = delete;
    SparseMatBasePARDISO& operator=(const SparseMatBasePARDISO&) = delete;
    SparseMatBasePARDISO& operator=(SparseMatBasePARDISO&&) noexcept = delete;

    ~SparseMatBasePARDISO() override { dealloc(); }

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SparseMatBasePARDISO>(*this); }
};

template<sp_d T, la_it MT> int SparseMatBasePARDISO<T, MT>::direct_solve(Mat<T>& X, const Mat<T>& B) {
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
    pardiso(pt, &maxfct, &mnum, &mtype, &PARDISO_SOLVE, &csr_mat.n_rows, csr_mat.val_mem(), csr_mat.row_mem(), csr_mat.col_mem(), nullptr, &nrhs, iparm, &msglvl, (void*)B.memptr(), X.memptr(), &info);

    if(0 != info) {
        suanpan_error("Error code {} received.\n", info);
        return SUANPAN_FAIL;
    }

    return SUANPAN_SUCCESS;
}

template<sp_d T> using SparseMatPARDISO = SparseMatBasePARDISO<T, 11>;
template<sp_d T> using SparseSymmMatPARDISO = SparseMatBasePARDISO<T, -2>;

#endif

//! @}
