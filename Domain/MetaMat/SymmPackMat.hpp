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
    int solve_trs(Mat<T>&, const Mat<T>&);

public:
    explicit SymmPackMat(uword);

    unique_ptr<MetaMat<T>> make_copy() override;

    void unify(uword) override;
    void nullify(uword) override;

    const T& operator()(uword, uword) const override;
    T& unsafe_at(uword, uword) override;
    T& at(uword, uword) override;

    Mat<T> operator*(const Mat<T>&) const override;

    int direct_solve(Mat<T>&, Mat<T>&&) override;
    int direct_solve(Mat<T>&, const Mat<T>&) override;
};

template<sp_d T> T SymmPackMat<T>::bin = 0.;

template<sp_d T> SymmPackMat<T>::SymmPackMat(const uword in_size)
    : DenseMat<T>(in_size, in_size, (in_size + 1) * in_size / 2)
    , length(2 * in_size - 1) {}

template<sp_d T> unique_ptr<MetaMat<T>> SymmPackMat<T>::make_copy() { return std::make_unique<SymmPackMat<T>>(*this); }

template<sp_d T> void SymmPackMat<T>::unify(const uword K) {
    nullify(K);
    access::rw(this->memory[(length - K + 2) * K / 2]) = 1.;
}

template<sp_d T> void SymmPackMat<T>::nullify(const uword K) {
    suanpan_for(0llu, K, [&](const uword I) { access::rw(this->memory[K + (length - I) * I / 2]) = 0.; });
    const auto t_factor = (length - K) * K / 2;
    suanpan_for(K, this->n_rows, [&](const uword I) { access::rw(this->memory[I + t_factor]) = 0.; });

    this->factored = false;
}

template<sp_d T> const T& SymmPackMat<T>::operator()(const uword in_row, const uword in_col) const { return this->memory[in_row >= in_col ? in_row + (length - in_col) * in_col / 2 : in_col + (length - in_row) * in_row / 2]; }

template<sp_d T> T& SymmPackMat<T>::unsafe_at(const uword in_row, const uword in_col) {
    this->factored = false;
    return access::rw(this->memory[in_row + (length - in_col) * in_col / 2]);
}

template<sp_d T> T& SymmPackMat<T>::at(const uword in_row, const uword in_col) {
    if(in_row < in_col) [[unlikely]] return bin;
    this->factored = false;
    return access::rw(this->memory[in_row + (length - in_col) * in_col / 2]);
}

template<sp_d T> Mat<T> SymmPackMat<T>::operator*(const Mat<T>& X) const {
    auto Y = Mat<T>(arma::size(X), fill::none);

    const auto N = static_cast<int>(this->n_rows);
    const auto INC = 1;
    T ALPHA = 1.;
    T BETA = 0.;

    if(std::is_same_v<T, float>) {
        using E = float;
        suanpan_for(0llu, X.n_cols, [&](const uword I) { arma_fortran(arma_sspmv)(&UPLO, &N, (E*)&ALPHA, (E*)this->memptr(), (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
    }
    else if(std::is_same_v<T, double>) {
        using E = double;
        suanpan_for(0llu, X.n_cols, [&](const uword I) { arma_fortran(arma_dspmv)(&UPLO, &N, (E*)&ALPHA, (E*)this->memptr(), (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
    }

    return Y;
}

template<sp_d T> int SymmPackMat<T>::direct_solve(Mat<T>& X, const Mat<T>& B) {
    if(this->factored) return this->solve_trs(X, B);

    const auto N = static_cast<int>(this->n_rows);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDB = static_cast<int>(B.n_rows);
    auto INFO = 0;

    this->factored = true;

    if(std::is_same_v<T, float>) {
        using E = float;
        X = B;
        arma_fortran(arma_sppsv)(&UPLO, &N, &NRHS, (E*)this->memptr(), (E*)X.memptr(), &LDB, &INFO);
    }
    else if(Precision::FULL == this->setting.precision) {
        using E = double;
        X = B;
        arma_fortran(arma_dppsv)(&UPLO, &N, &NRHS, (E*)this->memptr(), (E*)X.memptr(), &LDB, &INFO);
    }
    else {
        this->s_memory = this->to_float();
        arma_fortran(arma_spptrf)(&UPLO, &N, this->s_memory.memptr(), &INFO);
        if(0 == INFO) INFO = this->solve_trs(X, B);
    }

    if(0 != INFO) suanpan_error("solve() receives error code %u from the base driver, the matrix is probably singular.\n", INFO);

    return INFO;
}

template<sp_d T> int SymmPackMat<T>::solve_trs(Mat<T>& X, const Mat<T>& B) {
    const auto N = static_cast<int>(this->n_rows);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDB = static_cast<int>(B.n_rows);
    auto INFO = 0;

    if(std::is_same_v<T, float>) {
        using E = float;
        X = B;
        arma_fortran(arma_spptrs)(&UPLO, &N, &NRHS, (E*)this->memptr(), (E*)X.memptr(), &LDB, &INFO);
    }
    else if(Precision::FULL == this->setting.precision) {
        using E = double;
        X = B;
        arma_fortran(arma_dpptrs)(&UPLO, &N, &NRHS, (E*)this->memptr(), (E*)X.memptr(), &LDB, &INFO);
    }
    else {
        X = arma::zeros(B.n_rows, B.n_cols);

        mat full_residual = B;

        auto multiplier = norm(full_residual);

        auto counter = 0u;
        while(counter++ < this->setting.iterative_refinement) {
            if(multiplier < this->setting.tolerance) break;

            auto residual = conv_to<fmat>::from(full_residual / multiplier);

            arma_fortran(arma_spptrs)(&UPLO, &N, &NRHS, this->s_memory.memptr(), residual.memptr(), &LDB, &INFO);
            if(0 != INFO) break;

            const mat incre = multiplier * conv_to<mat>::from(residual);

            X += incre;

            suanpan_debug("mixed precision algorithm multiplier: %.5E.\n", multiplier = arma::norm(full_residual -= this->operator*(incre)));
        }
    }

    if(INFO != 0) suanpan_error("solve() receives error code %u from base driver, the matrix is probably singular.\n", INFO);

    return INFO;
}

template<sp_d T> int SymmPackMat<T>::direct_solve(Mat<T>& X, Mat<T>&& B) {
    if(this->factored) return this->solve_trs(X, std::forward<Mat<T>>(B));

    const auto N = static_cast<int>(this->n_rows);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDB = static_cast<int>(B.n_rows);
    auto INFO = 0;

    this->factored = true;

    if(std::is_same_v<T, float>) {
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

    if(0 != INFO) suanpan_error("solve() receives error code %u from the base driver, the matrix is probably singular.\n", INFO);

    return INFO;
}

template<sp_d T> int SymmPackMat<T>::solve_trs(Mat<T>& X, Mat<T>&& B) {
    const auto N = static_cast<int>(this->n_rows);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDB = static_cast<int>(B.n_rows);
    auto INFO = 0;

    if(std::is_same_v<T, float>) {
        using E = float;
        arma_fortran(arma_spptrs)(&UPLO, &N, &NRHS, (E*)this->memptr(), (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else if(Precision::FULL == this->setting.precision) {
        using E = double;
        arma_fortran(arma_dpptrs)(&UPLO, &N, &NRHS, (E*)this->memptr(), (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else {
        X = arma::zeros(B.n_rows, B.n_cols);

        auto multiplier = arma::norm(B);

        auto counter = 0u;
        while(counter++ < this->setting.iterative_refinement) {
            if(multiplier < this->setting.tolerance) break;

            auto residual = conv_to<fmat>::from(B / multiplier);

            arma_fortran(arma_spptrs)(&UPLO, &N, &NRHS, this->s_memory.memptr(), residual.memptr(), &LDB, &INFO);
            if(0 != INFO) break;

            const mat incre = multiplier * conv_to<mat>::from(residual);

            X += incre;

            suanpan_debug("mixed precision algorithm multiplier: %.5E.\n", multiplier = arma::norm(B -= this->operator*(incre)));
        }
    }

    if(INFO != 0) suanpan_error("solve() receives error code %u from base driver, the matrix is probably singular.\n", INFO);

    return INFO;
}

#endif

//! @}
