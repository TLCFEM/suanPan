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
 * @class FullMat
 * @brief A FullMat class that holds matrices.
 *
 * @author tlc
 * @date 06/09/2017
 * @version 0.1.0
 * @file FullMat.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef FULLMAT_HPP
#define FULLMAT_HPP

#include "DenseMat.hpp"

template<sp_d T> class FullMat : public DenseMat<T> {
    static constexpr char TRAN = 'N';

    int solve_trs(Mat<T>&, Mat<T>&&);

protected:
    using MetaMat<T>::direct_solve;

    int direct_solve(Mat<T>&, Mat<T>&&) override;

public:
    FullMat(const uword in_rows, const uword in_cols)
        : DenseMat<T>(in_rows, in_cols, in_rows * in_cols) {}

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<FullMat>(*this); }

    void nullify(const uword K) override {
        this->factored = false;
        suanpan_for(0llu, this->n_rows, [&](const uword I) { at(I, K) = T(0); });
        suanpan_for(0llu, this->n_cols, [&](const uword I) { at(K, I) = T(0); });
    }

    T operator()(const uword in_row, const uword in_col) const override { return this->memory[in_row + in_col * this->n_rows]; }

    T& at(const uword in_row, const uword in_col) override {
        this->factored = false;
        return this->memory[in_row + in_col * this->n_rows];
    }

    Mat<T> operator*(const Mat<T>&) const override;
};

template<sp_d T> Mat<T> FullMat<T>::operator*(const Mat<T>& B) const {
    Mat<T> C(arma::size(B));

    const auto M = static_cast<int>(this->n_rows);
    const auto N = static_cast<int>(this->n_cols);

    T ALPHA = T(1), BETA = T(0);

    if(1 == B.n_cols) {
        constexpr auto INC = 1;

        if constexpr(std::is_same_v<T, float>) {
            using E = float;
            arma_fortran(arma_sgemv)(&TRAN, &M, &N, (E*)&ALPHA, (E*)this->memptr(), &M, (E*)B.memptr(), &INC, (E*)&BETA, (E*)C.memptr(), &INC);
        }
        else {
            using E = double;
            arma_fortran(arma_dgemv)(&TRAN, &M, &N, (E*)&ALPHA, (E*)this->memptr(), &M, (E*)B.memptr(), &INC, (E*)&BETA, (E*)C.memptr(), &INC);
        }
    }
    else {
        const auto K = static_cast<int>(B.n_cols);

        if constexpr(std::is_same_v<T, float>) {
            using E = float;
            arma_fortran(arma_sgemm)(&TRAN, &TRAN, &M, &K, &N, (E*)&ALPHA, (E*)this->memptr(), &M, (E*)B.memptr(), &N, (E*)&BETA, (E*)C.memptr(), &M);
        }
        else {
            using E = double;
            arma_fortran(arma_dgemm)(&TRAN, &TRAN, &M, &K, &N, (E*)&ALPHA, (E*)this->memptr(), &M, (E*)B.memptr(), &N, (E*)&BETA, (E*)C.memptr(), &M);
        }
    }

    return C;
}

template<sp_d T> int FullMat<T>::direct_solve(Mat<T>& X, Mat<T>&& B) {
    if(this->factored) return this->solve_trs(X, std::forward<Mat<T>>(B));

    suanpan_assert([&] { if(this->n_rows != this->n_cols) throw invalid_argument("requires a square matrix"); });

    auto INFO = 0;

    auto N = static_cast<int>(this->n_rows);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDB = static_cast<int>(B.n_rows);
    this->pivot.zeros(N);
    this->factored = true;

    if constexpr(std::is_same_v<T, float>) {
        using E = float;
        arma_fortran(arma_sgesv)(&N, &NRHS, (E*)this->memptr(), &N, this->pivot.memptr(), (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else if(Precision::FULL == this->setting.precision) {
        using E = double;
        arma_fortran(arma_dgesv)(&N, &NRHS, (E*)this->memptr(), &N, this->pivot.memptr(), (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else {
        this->s_memory = this->to_float();
        arma_fortran(arma_sgetrf)(&N, &N, this->s_memory.memptr(), &N, this->pivot.memptr(), &INFO);
        if(0 == INFO) INFO = this->solve_trs(X, std::forward<Mat<T>>(B));
    }

    if(0 != INFO)
        suanpan_error("Error code {} received, the matrix is probably singular.\n", INFO);

    return INFO;
}

template<sp_d T> int FullMat<T>::solve_trs(Mat<T>& X, Mat<T>&& B) {
    auto INFO = 0;

    const auto N = static_cast<int>(this->n_rows);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDB = static_cast<int>(B.n_rows);

    if constexpr(std::is_same_v<T, float>) {
        using E = float;
        arma_fortran(arma_sgetrs)(&TRAN, &N, &NRHS, (E*)this->memptr(), &N, this->pivot.memptr(), (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else if(Precision::FULL == this->setting.precision) {
        using E = double;
        arma_fortran(arma_dgetrs)(&TRAN, &N, &NRHS, (E*)this->memptr(), &N, this->pivot.memptr(), (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else
        this->mixed_trs(X, std::forward<Mat<T>>(B), [&](fmat& residual) {
            arma_fortran(arma_sgetrs)(&TRAN, &N, &NRHS, this->s_memory.memptr(), &N, this->pivot.memptr(), residual.memptr(), &LDB, &INFO);
            return INFO;
        });

    return INFO;
}

#endif

//! @}
