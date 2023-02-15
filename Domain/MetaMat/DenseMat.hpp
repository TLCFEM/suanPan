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
 * @class DenseMat
 * @brief A DenseMat class that holds matrices.
 *
 * @author tlc
 * @date 19/04/2021
 * @version 0.1.0
 * @file DenseMat.hpp
 * @addtogroup MetaMat
 * @{
 */

#ifndef DENSEMAT_HPP
#define DENSEMAT_HPP

#include "MetaMat.hpp"

template<sp_d T> class DenseMat : public MetaMat<T> {
protected:
    using MetaMat<T>::direct_solve;

    int direct_solve(Mat<T>& X, const Mat<T>& B) override { return this->direct_solve(X, Mat<T>(B)); }

    podarray<int> pivot;
    podarray<float> s_memory; // float storage used in mixed precision algorithm

    std::unique_ptr<T[]> memory = nullptr;

    podarray<float> to_float() {
        podarray<float> f_memory(this->n_elem);
        suanpan_for(0llu, this->n_elem, [&](const uword I) { f_memory(I) = static_cast<float>(memory[I]); });
        return f_memory;
    }

public:
    DenseMat(const uword in_rows, const uword in_cols, const uword in_elem)
        : MetaMat<T>(in_rows, in_cols, in_elem)
        , memory(std::unique_ptr<T[]>(new T[this->n_elem])) { DenseMat::zeros(); }

    DenseMat(const DenseMat& old_mat)
        : MetaMat<T>(old_mat)
        , pivot(old_mat.pivot)
        , s_memory(old_mat.s_memory)
        , memory(std::unique_ptr<T[]>(new T[this->n_elem])) { suanpan_for(0llu, this->n_elem, [&](const uword I) { memory[I] = old_mat.memory[I]; }); }

    DenseMat(DenseMat&&) noexcept = delete;
    DenseMat& operator=(const DenseMat&) = delete;
    DenseMat& operator=(DenseMat&&) noexcept = delete;
    ~DenseMat() override = default;

    [[nodiscard]] bool is_empty() const override { return 0 == this->n_elem; }

    void zeros() override {
        this->factored = false;
        arrayops::fill_zeros(memptr(), this->n_elem);
    }

    [[nodiscard]] T max() const override {
        T max_value = T(1);
        for(uword I = 0; I < std::min(this->n_rows, this->n_cols); ++I) if(const auto t_val = this->operator()(I, I); t_val > max_value) max_value = t_val;
        return max_value;
    }

    [[nodiscard]] Col<T> diag() const override {
        Col<T> diag_vec(std::min(this->n_rows, this->n_cols), fill::none);
        suanpan_for(0llu, diag_vec.n_elem, [&](const uword I) { diag_vec(I) = this->operator()(I, I); });
        return diag_vec;
    }

    [[nodiscard]] const T* memptr() const override { return memory.get(); }

    T* memptr() override { return memory.get(); }

    void operator+=(const shared_ptr<MetaMat<T>>& M) override {
        if(nullptr == M) return;
        if(!M->triplet_mat.is_empty()) return this->operator+=(M->triplet_mat);
        if(this->n_rows != M->n_rows || this->n_cols != M->n_cols || this->n_elem != M->n_elem) throw invalid_argument("size mismatch");
        if(nullptr == M->memptr()) return;
        this->factored = false;
        arrayops::inplace_plus(memptr(), M->memptr(), this->n_elem);
    }

    void operator-=(const shared_ptr<MetaMat<T>>& M) override {
        if(nullptr == M) return;
        if(!M->triplet_mat.is_empty()) return this->operator-=(M->triplet_mat);
        if(this->n_rows != M->n_rows || this->n_cols != M->n_cols || this->n_elem != M->n_elem) throw invalid_argument("size mismatch");
        if(nullptr == M->memptr()) return;
        this->factored = false;
        arrayops::inplace_minus(memptr(), M->memptr(), this->n_elem);
    }

    void operator+=(const triplet_form<T, uword>& M) override {
        if(this->n_rows != M.n_rows || this->n_cols != M.n_cols) throw invalid_argument("size mismatch");
        this->factored = false;
        const auto row = M.row_mem();
        const auto col = M.col_mem();
        const auto val = M.val_mem();
        for(uword I = 0llu; I < M.n_elem; ++I) this->at(row[I], col[I]) += val[I];
    }

    void operator-=(const triplet_form<T, uword>& M) override {
        if(this->n_rows != M.n_rows || this->n_cols != M.n_cols) throw invalid_argument("size mismatch");
        this->factored = false;
        const auto row = M.row_mem();
        const auto col = M.col_mem();
        const auto val = M.val_mem();
        for(uword I = 0llu; I < M.n_elem; ++I) this->at(row[I], col[I]) -= val[I];
    }

    void operator*=(const T value) override {
        this->factored = false;
        arrayops::inplace_mul(memptr(), value, this->n_elem);
    }

    [[nodiscard]] int sign_det() const override {
        if(IterativeSolver::NONE != this->setting.iterative_solver) throw invalid_argument("analysis requires the sign of determinant but iterative solver does not support it");
        auto det_sign = 1;
        for(unsigned I = 0; I < pivot.n_elem; ++I) if((this->operator()(I, I) < T(0)) ^ (static_cast<int>(I) + 1 != pivot(I))) det_sign = -det_sign;
        return det_sign;
    }
};

#endif

//! @}
