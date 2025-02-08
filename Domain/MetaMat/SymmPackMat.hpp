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
 * @class SymmPackMat
 * @brief A SymmPackMat class that holds matrices.
 *
 * @author tlc
 * @date 13/11/2022
 * @version 0.2.0
 * @file SymmPackMat.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef SYMMPACKMAT_HPP
#define SYMMPACKMAT_HPP

#include "DenseMat.hpp"

template<sp_d T> class SymmPackMat final : public DenseMat<T> {
    static constexpr char UPLO = 'L';

    static T bin;

    const uword length; // 2n-1

    int solve_trs(Mat<T>&, Mat<T>&&);

protected:
    using DenseMat<T>::direct_solve;

    int direct_solve(Mat<T>&, Mat<T>&&) override;

public:
    explicit SymmPackMat(const uword in_size)
        : DenseMat<T>(in_size, in_size, (in_size + 1) * in_size / 2)
        , length(2 * in_size - 1) {}

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SymmPackMat>(*this); }

    void nullify(const uword K) override {
        this->factored = false;
        suanpan::for_each(K, [&](const uword I) { this->memory[K + (length - I) * I / 2] = T(0); });
        const auto t_factor = (length - K) * K / 2;
        suanpan::for_each(K, this->n_rows, [&](const uword I) { this->memory[I + t_factor] = T(0); });
    }

    T operator()(const uword in_row, const uword in_col) const override { return this->memory[in_row >= in_col ? in_row + (length - in_col) * in_col / 2 : in_col + (length - in_row) * in_row / 2]; }

    T& unsafe_at(const uword in_row, const uword in_col) override {
        this->factored = false;
        return this->memory[in_row + (length - in_col) * in_col / 2];
    }

    T& at(const uword in_row, const uword in_col) override {
        if(in_row < in_col) [[unlikely]] return bin = T(0);
        return this->unsafe_at(in_row, in_col);
    }

    Mat<T> operator*(const Mat<T>&) const override;
};

template<sp_d T> T SymmPackMat<T>::bin = T(0);

template<sp_d T> Mat<T> SymmPackMat<T>::operator*(const Mat<T>& X) const {
    auto Y = Mat<T>(arma::size(X), fill::none);

    const auto N = static_cast<int>(this->n_rows);
    constexpr auto INC = 1;
    T ALPHA = T(1);
    T BETA = T(0);

    if constexpr(std::is_same_v<T, float>) {
        using E = float;
        suanpan::for_each(X.n_cols, [&](const uword I) { arma_fortran(arma_sspmv)(&UPLO, &N, (E*)&ALPHA, (E*)this->memptr(), (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
    }
    else {
        using E = double;
        suanpan::for_each(X.n_cols, [&](const uword I) { arma_fortran(arma_dspmv)(&UPLO, &N, (E*)&ALPHA, (E*)this->memptr(), (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
    }

    return Y;
}

template<sp_d T> int SymmPackMat<T>::direct_solve(Mat<T>& X, Mat<T>&& B) {
    if(this->factored) return this->solve_trs(X, std::forward<Mat<T>>(B));

    suanpan_assert([&] { if(this->n_rows != this->n_cols) throw invalid_argument("requires a square matrix"); });

    auto INFO = 0;

    const auto N = static_cast<int>(this->n_rows);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDB = static_cast<int>(B.n_rows);
    this->factored = true;

    if constexpr(std::is_same_v<T, float>) {
        using E = float;
        arma_fortran(arma_sppsv)(&UPLO, &N, &NRHS, (E*)this->memptr(), (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else if(Precision::FULL == this->setting.precision) {
        using E = double;
        arma_fortran(arma_dppsv)(&UPLO, &N, &NRHS, (E*)this->memptr(), (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else {
        this->s_memory = this->to_float();
        arma_fortran(arma_spptrf)(&UPLO, &N, this->s_memory.memptr(), &INFO);
        if(0 == INFO) INFO = this->solve_trs(X, std::forward<Mat<T>>(B));
    }

    if(0 != INFO)
        suanpan_error("Error code {} received, the matrix is probably singular.\n", INFO);

    return INFO;
}

template<sp_d T> int SymmPackMat<T>::solve_trs(Mat<T>& X, Mat<T>&& B) {
    auto INFO = 0;

    const auto N = static_cast<int>(this->n_rows);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDB = static_cast<int>(B.n_rows);

    if constexpr(std::is_same_v<T, float>) {
        using E = float;
        arma_fortran(arma_spptrs)(&UPLO, &N, &NRHS, (E*)this->memptr(), (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else if(Precision::FULL == this->setting.precision) {
        using E = double;
        arma_fortran(arma_dpptrs)(&UPLO, &N, &NRHS, (E*)this->memptr(), (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else
        this->mixed_trs(X, std::forward<Mat<T>>(B), [&](fmat& residual) {
            arma_fortran(arma_spptrs)(&UPLO, &N, &NRHS, this->s_memory.memptr(), residual.memptr(), &LDB, &INFO);
            return INFO;
        });

    if(0 != INFO)
        suanpan_error("Error code {} received, the matrix is probably singular.\n", INFO);

    return INFO;
}

#endif

//! @}
