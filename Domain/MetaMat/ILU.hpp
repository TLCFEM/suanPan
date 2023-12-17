/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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
 * @class ILU
 * @brief A ILU class.
 *
 * @author tlc
 * @date 24/07/2022
 * @version 0.1.0
 * @file ILU.hpp
 * @addtogroup Preconditioner
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef ILU_HPP
#define ILU_HPP

#ifndef SUANPAN_SUPERLUMT

#include <superlu-mt/superlu-mt.h>
#include "Preconditioner.hpp"
#include "csc_form.hpp"
#include "triplet_form.hpp"

template<sp_d data_t> class ILU final : public Preconditioner<data_t> {
    inline static char equed = 'N';

    SuperMatrix A{}, L{}, U{};

    SuperLUStat_t stat{};
    GlobalLU_t global_setting{};
    superlu_options_t options{};
    mem_usage_t usage{};

    void* t_val = nullptr;
    int* t_row = nullptr;
    int* t_col = nullptr;

    int* perm_r = nullptr;
    int* perm_c = nullptr;

    void* scale_r = nullptr;
    void* scale_c = nullptr;

    int* etree = nullptr;

    static void wrap_b(SuperMatrix*, const Col<data_t>&);

    int apply_solver(SuperMatrix*, SuperMatrix*);

    void apply_row_scale(Col<data_t>&);
    void apply_column_scale(Col<data_t>&);

    template<sp_i index_t> void allocate(triplet_form<data_t, index_t>&);

public:
    template<sp_i index_t> explicit ILU(triplet_form<data_t, index_t>&&);
    template<sp_i index_t> explicit ILU(triplet_form<data_t, index_t>&);

    ~ILU() override;

    int init() override;

    [[nodiscard]] Col<data_t> apply(const Col<data_t>&) override;
};

template<sp_d data_t> void ILU<data_t>::wrap_b(SuperMatrix* out, const Col<data_t>& in) {
    if(std::is_same_v<data_t, float>) {
        using E = float;
        sCreate_Dense_Matrix(out, (int)in.n_rows, (int)in.n_cols, (E*)in.memptr(), (int)in.n_rows, Stype_t::SLU_DN, Dtype_t::SLU_S, Mtype_t::SLU_GE); // NOLINT(clang-diagnostic-cast-qual)
    }
    else {
        using E = double;
        dCreate_Dense_Matrix(out, (int)in.n_rows, (int)in.n_cols, (E*)in.memptr(), (int)in.n_rows, Stype_t::SLU_DN, Dtype_t::SLU_D, Mtype_t::SLU_GE); // NOLINT(clang-diagnostic-cast-qual)
    }
}

template<sp_d data_t> int ILU<data_t>::apply_solver(SuperMatrix* X, SuperMatrix* B) {
    int info;

    if(std::is_same_v<data_t, float>) {
        using E = float;
        sgsisx(&options, &A, perm_c, perm_r, etree, &equed, (E*)scale_r, (E*)scale_c, &L, &U, nullptr, 0, B, X, nullptr, nullptr, &global_setting, &usage, &stat, &info);
    }
    else {
        using E = double;
        dgsisx(&options, &A, perm_c, perm_r, etree, &equed, (E*)scale_r, (E*)scale_c, &L, &U, nullptr, 0, B, X, nullptr, nullptr, &global_setting, &usage, &stat, &info);
    }

    return info;
}

template<sp_d data_t> void ILU<data_t>::apply_row_scale(Col<data_t>& in_mat) { in_mat %= Col<data_t>((data_t*)scale_r, in_mat.n_elem); }

template<sp_d data_t> void ILU<data_t>::apply_column_scale(Col<data_t>& in_mat) { in_mat %= Col<data_t>((data_t*)scale_c, in_mat.n_elem); }

template<sp_d data_t> template<sp_i index_t> void ILU<data_t>::allocate(triplet_form<data_t, index_t>& triplet_mat) {
    csc_form<data_t, int> csc_mat(triplet_mat);

    auto t_size = sizeof(data_t) * csc_mat.n_elem;
    t_val = superlu_malloc(t_size);
    std::memcpy(t_val, (void*)csc_mat.val_mem(), t_size);

    t_size = sizeof(int) * csc_mat.n_elem;
    t_row = (int*)superlu_malloc(t_size);
    std::memcpy(t_row, (void*)csc_mat.row_mem(), t_size);

    t_size = sizeof(int) * (csc_mat.n_cols + 1llu);
    t_col = (int*)superlu_malloc(t_size);
    std::memcpy(t_col, (void*)csc_mat.col_mem(), t_size);

    if(std::is_same_v<data_t, float>) {
        using E = float;
        sCreate_CompCol_Matrix(&A, csc_mat.n_rows, csc_mat.n_cols, csc_mat.n_elem, (E*)t_val, t_row, t_col, Stype_t::SLU_NC, Dtype_t::SLU_S, Mtype_t::SLU_GE);
    }
    else {
        using E = double;
        dCreate_CompCol_Matrix(&A, csc_mat.n_rows, csc_mat.n_cols, csc_mat.n_elem, (E*)t_val, t_row, t_col, Stype_t::SLU_NC, Dtype_t::SLU_D, Mtype_t::SLU_GE);
    }

    perm_r = (int*)superlu_malloc(sizeof(int) * (csc_mat.n_rows + 1llu));
    perm_c = (int*)superlu_malloc(sizeof(int) * (csc_mat.n_cols + 1llu));
    scale_r = superlu_malloc(sizeof(data_t) * (csc_mat.n_rows + 1llu));
    scale_c = superlu_malloc(sizeof(data_t) * (csc_mat.n_cols + 1llu));

    etree = (int*)superlu_malloc(sizeof(int) * (csc_mat.n_cols + 1llu));

    ilu_set_default_options(&options);

    StatInit(&stat);
}

template<sp_d data_t> template<sp_i index_t> ILU<data_t>::ILU(triplet_form<data_t, index_t>&& triplet_mat)
    : Preconditioner<data_t>() { this->allocate(triplet_mat); }

template<sp_d data_t> template<sp_i index_t> ILU<data_t>::ILU(triplet_form<data_t, index_t>& triplet_mat)
    : Preconditioner<data_t>() { this->allocate(triplet_mat); }

template<sp_d data_t> ILU<data_t>::~ILU() {
    Destroy_SuperMatrix_Store(&A);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);

    if(t_val) superlu_free(t_val);
    if(t_row) superlu_free(t_row);
    if(t_col) superlu_free(t_col);
    if(etree) superlu_free(etree);
    if(perm_r) superlu_free(perm_r);
    if(perm_c) superlu_free(perm_c);
    if(scale_r) superlu_free(scale_r);
    if(scale_c) superlu_free(scale_c);
}

template<sp_d data_t> int ILU<data_t>::init() {
    Col<data_t> in(A.ncol, 1, fill::none);
    Col<data_t> out(A.nrow, 1, fill::none);

    SuperMatrix B, X;
    this->wrap_b(&X, out);
    this->wrap_b(&B, in);

    B.ncol = 0;

    if(const auto code = this->apply_solver(&X, &B); 0 != code) {
        suanpan_error("Error code {} received.\n", code);
        return SUANPAN_FAIL;
    }

    Destroy_SuperMatrix_Store(&X);
    Destroy_SuperMatrix_Store(&B);

    options.Fact = superlu::FACTORED;

    return SUANPAN_SUCCESS;
}

template<sp_d data_t> Col<data_t> ILU<data_t>::apply(const Col<data_t>& in) {
    Col<data_t> out(arma::size(in), fill::none);

    SuperMatrix X, B;
    this->wrap_b(&X, out);
    this->wrap_b(&B, in);

    this->apply_solver(&X, &B);

    // this->apply_column_scale(out);

    Destroy_SuperMatrix_Store(&X);
    Destroy_SuperMatrix_Store(&B);

    return out;
}

#endif

#endif

//! @}
