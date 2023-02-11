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
    podarray<int> pivot;
    podarray<float> s_memory; // float storage used in mixed precision algorithm

    std::unique_ptr<T[]> memory = nullptr;

    podarray<float> to_float();

public:
    DenseMat(uword, uword, uword);
    DenseMat(const DenseMat&);
    DenseMat(DenseMat&&) noexcept;
    DenseMat& operator=(const DenseMat&);
    DenseMat& operator=(DenseMat&&) noexcept;
    ~DenseMat() override = default;

    void unify(uword) override;

    [[nodiscard]] bool is_empty() const override;
    void zeros() override;

    [[nodiscard]] T max() const override;
    [[nodiscard]] Col<T> diag() const override;

    const T* memptr() const override;
    T* memptr() override;

    void operator+=(const shared_ptr<MetaMat<T>>&) override;
    void operator-=(const shared_ptr<MetaMat<T>>&) override;

    void operator+=(const triplet_form<T, uword>&) override;
    void operator-=(const triplet_form<T, uword>&) override;

    void operator*=(T) override;

    [[nodiscard]] int sign_det() const override;
};

template<sp_d T> podarray<float> DenseMat<T>::to_float() {
    podarray<float> f_memory(this->n_elem);

    suanpan_for(0llu, this->n_elem, [&](const uword I) { f_memory(I) = static_cast<float>(memory[I]); });

    return f_memory;
}

template<sp_d T> DenseMat<T>::DenseMat(const uword in_rows, const uword in_cols, const uword in_elem)
    : MetaMat<T>(in_rows, in_cols, in_elem)
    , memory(std::unique_ptr<T[]>(new T[this->n_elem])) { DenseMat::zeros(); }

template<sp_d T> DenseMat<T>::DenseMat(const DenseMat& old_mat)
    : MetaMat<T>(old_mat)
    , pivot(old_mat.pivot)
    , s_memory(old_mat.s_memory)
    , memory(std::unique_ptr<T[]>(new T[this->n_elem])) { std::copy(old_mat.memory.get(), old_mat.memory.get() + old_mat.n_elem, DenseMat::memptr()); }

template<sp_d T> DenseMat<T>::DenseMat(DenseMat&& old_mat) noexcept
    : MetaMat<T>(std::move(old_mat))
    , pivot(std::move(old_mat.pivot))
    , s_memory(std::move(old_mat.s_memory))
    , memory(std::move(old_mat.memory)) {}

template<sp_d T> DenseMat<T>& DenseMat<T>::operator=(const DenseMat& old_mat) {
    if(this == &old_mat) return *this;
    MetaMat<T>::operator=(old_mat);
    pivot = old_mat.pivot;
    s_memory = old_mat.s_memory;
    memory = std::unique_ptr<T[]>(new T[this->n_elem]);
    std::copy(old_mat.memory.get(), old_mat.memory.get() + old_mat.n_elem, memptr());
    return *this;
}

template<sp_d T> DenseMat<T>& DenseMat<T>::operator=(DenseMat&& old_mat) noexcept {
    if(this == &old_mat) return *this;
    MetaMat<T>::operator=(std::move(old_mat));
    pivot = std::move(old_mat.pivot);
    s_memory = std::move(old_mat.s_memory);
    memory = std::move(old_mat.memory);
    return *this;
}

template<sp_d T> void DenseMat<T>::unify(const uword K) {
    this->nullify(K);
    this->at(K, K) = 1.;
}

template<sp_d T> bool DenseMat<T>::is_empty() const { return 0 == this->n_elem; }

template<sp_d T> void DenseMat<T>::zeros() {
    arrayops::fill_zeros(memptr(), this->n_elem);
    this->factored = false;
}

template<sp_d T> T DenseMat<T>::max() const {
    T max_value = T(1);
    for(uword I = 0; I < std::min(this->n_rows, this->n_cols); ++I) if(const auto t_val = this->operator()(I, I); t_val > max_value) max_value = t_val;
    return max_value;
}

template<sp_d T> Col<T> DenseMat<T>::diag() const {
    Col<T> diag_vec(std::min(this->n_rows, this->n_cols), fill::none);

    suanpan_for(0llu, diag_vec.n_elem, [&](const uword I) { diag_vec(I) = this->operator()(I, I); });

    return diag_vec;
}

template<sp_d T> const T* DenseMat<T>::memptr() const { return memory.get(); }

template<sp_d T> T* DenseMat<T>::memptr() { return memory.get(); }

template<sp_d T> void DenseMat<T>::operator+=(const shared_ptr<MetaMat<T>>& M) {
    if(nullptr == M) return;
    if(!M->triplet_mat.is_empty()) return this->operator+=(M->triplet_mat);
    if(this->n_rows != M->n_rows || this->n_cols != M->n_cols || this->n_elem != M->n_elem) return;
    if(nullptr == M->memptr()) return;
    arrayops::inplace_plus(memptr(), M->memptr(), this->n_elem);
    this->factored = false;
}

template<sp_d T> void DenseMat<T>::operator-=(const shared_ptr<MetaMat<T>>& M) {
    if(nullptr == M) return;
    if(!M->triplet_mat.is_empty()) return this->operator-=(M->triplet_mat);
    if(this->n_rows != M->n_rows || this->n_cols != M->n_cols || this->n_elem != M->n_elem) return;
    if(nullptr == M->memptr()) return;
    arrayops::inplace_minus(memptr(), M->memptr(), this->n_elem);
    this->factored = false;
}

template<sp_d T> void DenseMat<T>::operator+=(const triplet_form<T, uword>& M) {
    if(this->n_rows != M.n_rows || this->n_cols != M.n_cols) return;

    const auto row = M.row_mem();
    const auto col = M.col_mem();
    const auto val = M.val_mem();
    for(uword I = 0llu; I < M.n_elem; ++I) this->at(row[I], col[I]) += val[I];
}

template<sp_d T> void DenseMat<T>::operator-=(const triplet_form<T, uword>& M) {
    if(this->n_rows != M.n_rows || this->n_cols != M.n_cols) return;

    const auto row = M.row_mem();
    const auto col = M.col_mem();
    const auto val = M.val_mem();
    for(uword I = 0llu; I < M.n_elem; ++I) this->at(row[I], col[I]) -= val[I];
}

template<sp_d T> void DenseMat<T>::operator*=(const T value) { arrayops::inplace_mul(memptr(), value, this->n_elem); }

template<sp_d T> int DenseMat<T>::sign_det() const {
    if(IterativeSolver::NONE != this->setting.iterative_solver) throw invalid_argument("analysis requires the sign of determinant but iterative solver does not support it");
    auto det_sign = 1;
    for(unsigned I = 0; I < pivot.n_elem; ++I) if((this->operator()(I, I) < 0.) ^ (static_cast<int>(I) + 1 != pivot(I))) det_sign = -det_sign;
    return det_sign;
}

#endif

//! @}
