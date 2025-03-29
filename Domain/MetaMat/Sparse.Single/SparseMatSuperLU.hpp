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
 * @class SparseMatSuperLU
 * @brief A SparseMatSuperLU class that holds matrices.
 *
 * @author tlc
 * @date 14/08/2020
 * @version 0.1.0
 * @file SparseMatSuperLU.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef SPARSEMATSUPERLU_HPP
#define SPARSEMATSUPERLU_HPP

#include "../SparseMat.hpp"
#include "../csc_form.hpp"

#include <superlu-mt/superlu-mt.h>

template<sp_d T> class SparseMatSuperLU final : public SparseMat<T> {
    SuperMatrix A{}, L{}, U{}, B{};

#ifndef SUANPAN_SUPERLUMT
    superlu_options_t options{};

    SuperLUStat_t stat{};
#else
    const int ordering_num = 1;

    Gstat_t stat{};
#endif

    void* t_val = nullptr;
    int* t_row = nullptr;
    int* t_col = nullptr;

    int* perm_r = nullptr;
    int* perm_c = nullptr;

    bool allocated = false;

    auto init_config();

    template<sp_d ET> void alloc(csc_form<ET, int>&&);
    void dealloc();

    template<sp_d ET> void wrap_b(const Mat<ET>&);
    template<sp_d ET> void tri_solve(int&);
    template<sp_d ET> void full_solve(int&);

    int solve_trs(Mat<T>&, Mat<T>&&);

protected:
    int direct_solve(Mat<T>& out_mat, const Mat<T>& in_mat) override { return this->direct_solve(out_mat, Mat<T>(in_mat)); }

    int direct_solve(Mat<T>&, Mat<T>&&) override;

public:
    SparseMatSuperLU(uword, uword, uword = 0);
    SparseMatSuperLU(const SparseMatSuperLU&);
    SparseMatSuperLU(SparseMatSuperLU&&) noexcept = delete;
    SparseMatSuperLU& operator=(const SparseMatSuperLU&) = delete;
    SparseMatSuperLU& operator=(SparseMatSuperLU&&) noexcept = delete;
    ~SparseMatSuperLU() override;

    unique_ptr<MetaMat<T>> make_copy() override;
};

template<sp_d T> auto SparseMatSuperLU<T>::init_config() {
#ifndef SUANPAN_SUPERLUMT
    set_default_options(&options);
    options.IterRefine = std::is_same_v<T, float> ? superlu::IterRefine_t::SLU_SINGLE : superlu::IterRefine_t::SLU_DOUBLE;
    options.Equil = superlu::yes_no_t::NO;

    arrayops::fill_zeros(reinterpret_cast<char*>(&stat), sizeof(SuperLUStat_t));

    StatInit(&stat);
#else
    StatAlloc(static_cast<int>(this->n_cols), SUANPAN_NUM_THREADS, sp_ienv(1), sp_ienv(2), &stat);
    StatInit(static_cast<int>(this->n_cols), SUANPAN_NUM_THREADS, &stat);
#endif
}

template<sp_d T> template<sp_d ET> void SparseMatSuperLU<T>::alloc(csc_form<ET, int>&& in) {
    dealloc();

    auto t_size = sizeof(ET) * in.n_elem;
    t_val = superlu_malloc(t_size);
    memcpy(t_val, (void*)in.val_mem(), t_size);

    t_size = sizeof(int) * in.n_elem;
    t_row = (int*)superlu_malloc(t_size);
    memcpy(t_row, (void*)in.row_mem(), t_size);

    t_size = sizeof(int) * (in.n_cols + 1llu);
    t_col = (int*)superlu_malloc(t_size);
    memcpy(t_col, (void*)in.col_mem(), t_size);

    if constexpr(std::is_same_v<ET, double>) {
        using E = double;
        dCreate_CompCol_Matrix(&A, in.n_rows, in.n_cols, in.n_elem, (E*)t_val, t_row, t_col, Stype_t::SLU_NC, Dtype_t::SLU_D, Mtype_t::SLU_GE);
    }
    else {
        using E = float;
        sCreate_CompCol_Matrix(&A, in.n_rows, in.n_cols, in.n_elem, (E*)t_val, t_row, t_col, Stype_t::SLU_NC, Dtype_t::SLU_S, Mtype_t::SLU_GE);
    }

    perm_r = (int*)superlu_malloc(sizeof(int) * (this->n_rows + 1));
    perm_c = (int*)superlu_malloc(sizeof(int) * (this->n_cols + 1));

    allocated = true;
}

template<sp_d T> void SparseMatSuperLU<T>::dealloc() {
    if(!allocated) return;

    Destroy_SuperMatrix_Store(&A);
#ifdef SUANPAN_SUPERLUMT
    Destroy_SuperNode_SCP(&L);
    Destroy_CompCol_NCP(&U);
#else
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
#endif

    if(t_val) superlu_free(t_val);
    if(t_row) superlu_free(t_row);
    if(t_col) superlu_free(t_col);
    if(perm_r) superlu_free(perm_r);
    if(perm_c) superlu_free(perm_c);

    allocated = false;
}

template<sp_d T> template<sp_d ET> void SparseMatSuperLU<T>::wrap_b(const Mat<ET>& in_mat) {
    if constexpr(std::is_same_v<ET, float>) {
        using E = float;
        sCreate_Dense_Matrix(&B, (int)in_mat.n_rows, (int)in_mat.n_cols, (E*)in_mat.memptr(), (int)in_mat.n_rows, Stype_t::SLU_DN, Dtype_t::SLU_S, Mtype_t::SLU_GE);
    }
    else {
        using E = double;
        dCreate_Dense_Matrix(&B, (int)in_mat.n_rows, (int)in_mat.n_cols, (E*)in_mat.memptr(), (int)in_mat.n_rows, Stype_t::SLU_DN, Dtype_t::SLU_D, Mtype_t::SLU_GE);
    }
}

template<sp_d T> template<sp_d ET> void SparseMatSuperLU<T>::tri_solve(int& flag) {
#ifdef SUANPAN_SUPERLUMT
    if(std::is_same_v<ET, float>) sgstrs(NOTRANS, &L, &U, perm_c, perm_r, &B, &stat, &flag);
    else dgstrs(NOTRANS, &L, &U, perm_c, perm_r, &B, &stat, &flag);
#else
    superlu::gstrs<ET>(options.Trans, &L, &U, perm_c, perm_r, &B, &stat, &flag);
#endif

    Destroy_SuperMatrix_Store(&B);
}

template<sp_d T> template<sp_d ET> void SparseMatSuperLU<T>::full_solve(int& flag) {
#ifdef SUANPAN_SUPERLUMT
    get_perm_c(ordering_num, &A, perm_c);
    if(std::is_same_v<ET, float>) psgssv(SUANPAN_NUM_THREADS, &A, perm_c, perm_r, &L, &U, &B, &flag);
    else pdgssv(SUANPAN_NUM_THREADS, &A, perm_c, perm_r, &L, &U, &B, &flag);
#else
    superlu::gssv<ET>(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &flag);
#endif

    Destroy_SuperMatrix_Store(&B);
}

template<sp_d T> SparseMatSuperLU<T>::SparseMatSuperLU(const uword in_row, const uword in_col, const uword in_elem)
    : SparseMat<T>(in_row, in_col, in_elem) { init_config(); }

template<sp_d T> SparseMatSuperLU<T>::SparseMatSuperLU(const SparseMatSuperLU& other)
    : SparseMat<T>(other) {
    init_config();
    this->factored = false;
}

template<sp_d T> SparseMatSuperLU<T>::~SparseMatSuperLU() {
    dealloc();
    StatFree(&stat);
}

template<sp_d T> unique_ptr<MetaMat<T>> SparseMatSuperLU<T>::make_copy() { return std::make_unique<SparseMatSuperLU>(*this); }

template<sp_d T> int SparseMatSuperLU<T>::direct_solve(Mat<T>& out_mat, Mat<T>&& in_mat) {
    if(this->factored) return solve_trs(out_mat, std::forward<Mat<T>>(in_mat));

    this->factored = true;

    alloc(csc_form<T, int>(this->triplet_mat));

    wrap_b(in_mat);

    auto flag = 0;

    full_solve<T>(flag);

    out_mat = std::move(in_mat);

    return flag;
}

template<sp_d T> int SparseMatSuperLU<T>::solve_trs(Mat<T>& out_mat, Mat<T>&& in_mat) {
    wrap_b(in_mat);

    auto flag = 0;

    tri_solve<T>(flag);

    out_mat = std::move(in_mat);

    return flag;
}
#endif

//! @}
