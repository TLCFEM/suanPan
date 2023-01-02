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
    int solve_trs(Mat<T>&, const Mat<T>&);

public:
    FullMat(uword, uword);

    unique_ptr<MetaMat<T>> make_copy() override;

    void unify(uword) override;
    void nullify(uword) override;

    const T& operator()(uword, uword) const override;

    T& at(uword, uword) override;

    Mat<T> operator*(const Mat<T>&) const override;

    int direct_solve(Mat<T>&, Mat<T>&&) override;
    int direct_solve(Mat<T>&, const Mat<T>&) override;
};

template<sp_d T> FullMat<T>::FullMat(const uword in_rows, const uword in_cols)
    : DenseMat<T>(in_rows, in_cols, in_rows * in_cols) {}

template<sp_d T> unique_ptr<MetaMat<T>> FullMat<T>::make_copy() { return std::make_unique<FullMat<T>>(*this); }

template<sp_d T> void FullMat<T>::unify(const uword K) {
    nullify(K);
    at(K, K) = 1.;
}

template<sp_d T> void FullMat<T>::nullify(const uword K) {
    suanpan_for(0llu, this->n_rows, [&](const uword I) { at(I, K) = 0.; });
    suanpan_for(0llu, this->n_cols, [&](const uword I) { at(K, I) = 0.; });

    this->factored = false;
}

template<sp_d T> const T& FullMat<T>::operator()(const uword in_row, const uword in_col) const { return this->memory[in_row + in_col * this->n_rows]; }

template<sp_d T> T& FullMat<T>::at(const uword in_row, const uword in_col) {
    this->factored = false;
    return access::rw(this->memory[in_row + in_col * this->n_rows]);
}

template<sp_d T> Mat<T> FullMat<T>::operator*(const Mat<T>& B) const {
    Mat<T> C(arma::size(B));

    const auto M = static_cast<int>(this->n_rows);
    const auto N = static_cast<int>(this->n_cols);

    T ALPHA = 1., BETA = 0.;

    if(1 == B.n_cols) {
        constexpr auto INCX = 1, INCY = 1;

        if(std::is_same_v<T, float>) {
            using E = float;
            arma_fortran(arma_sgemv)(&TRAN, &M, &N, (E*)&ALPHA, (E*)this->memptr(), &M, (E*)B.memptr(), &INCX, (E*)&BETA, (E*)C.memptr(), &INCY);
        }
        else if(std::is_same_v<T, double>) {
            using E = double;
            arma_fortran(arma_dgemv)(&TRAN, &M, &N, (E*)&ALPHA, (E*)this->memptr(), &M, (E*)B.memptr(), &INCX, (E*)&BETA, (E*)C.memptr(), &INCY);
        }
    }
    else {
        const auto K = static_cast<int>(B.n_cols);

        if(std::is_same_v<T, float>) {
            using E = float;
            arma_fortran(arma_sgemm)(&TRAN, &TRAN, &M, &K, &N, (E*)&ALPHA, (E*)this->memptr(), &M, (E*)B.memptr(), &N, (E*)&BETA, (E*)C.memptr(), &M);
        }
        else if(std::is_same_v<T, double>) {
            using E = double;
            arma_fortran(arma_dgemm)(&TRAN, &TRAN, &M, &K, &N, (E*)&ALPHA, (E*)this->memptr(), &M, (E*)B.memptr(), &N, (E*)&BETA, (E*)C.memptr(), &M);
        }
    }

    return C;
}

template<sp_d T> int FullMat<T>::direct_solve(Mat<T>& X, const Mat<T>& B) {
    if(this->factored) return this->solve_trs(X, B);

    auto N = static_cast<int>(this->n_rows);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDB = static_cast<int>(B.n_rows);
    auto INFO = 0;
    this->pivot.zeros(N);
    this->factored = true;

    if(std::is_same_v<T, float>) {
        using E = float;
        X = B;
        arma_fortran(arma_sgesv)(&N, &NRHS, (E*)this->memptr(), &N, this->pivot.memptr(), (E*)X.memptr(), &LDB, &INFO);
    }
    else if(Precision::FULL == this->setting.precision) {
        using E = double;
        X = B;
        arma_fortran(arma_dgesv)(&N, &NRHS, (E*)this->memptr(), &N, this->pivot.memptr(), (E*)X.memptr(), &LDB, &INFO);
    }
    else {
        this->s_memory = this->to_float();
        arma_fortran(arma_sgetrf)(&N, &N, this->s_memory.memptr(), &N, this->pivot.memptr(), &INFO);
        if(0 == INFO) INFO = this->solve_trs(X, B);
    }

    if(0 != INFO) SP_E("Error code {} received, the matrix is probably singular.\n", INFO);

    return INFO;
}

template<sp_d T> int FullMat<T>::solve_trs(Mat<T>& X, const Mat<T>& B) {
    const auto N = static_cast<int>(this->n_rows);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDB = static_cast<int>(B.n_rows);
    auto INFO = 0;

    if(std::is_same_v<T, float>) {
        using E = float;
        X = B;
        arma_fortran(arma_sgetrs)(&TRAN, &N, &NRHS, (E*)this->memptr(), &N, this->pivot.memptr(), (E*)X.memptr(), &LDB, &INFO);
    }
    else if(Precision::FULL == this->setting.precision) {
        using E = double;
        X = B;
        arma_fortran(arma_dgetrs)(&TRAN, &N, &NRHS, (E*)this->memptr(), &N, this->pivot.memptr(), (E*)X.memptr(), &LDB, &INFO);
    }
    else {
        X = arma::zeros(B.n_rows, B.n_cols);

        mat full_residual = B;

        auto multiplier = norm(full_residual);

        auto counter = 0u;
        while(counter++ < this->setting.iterative_refinement) {
            if(multiplier < this->setting.tolerance) break;

            auto residual = conv_to<fmat>::from(full_residual / multiplier);

            arma_fortran(arma_sgetrs)(&TRAN, &N, &NRHS, this->s_memory.memptr(), &N, this->pivot.memptr(), residual.memptr(), &LDB, &INFO);
            if(0 != INFO) break;

            const mat incre = multiplier * conv_to<mat>::from(residual);

            X += incre;

            suanpan_debug("mixed precision algorithm multiplier: %.5E.\n", multiplier = arma::norm(full_residual -= this->operator*(incre)));
        }
    }

    return INFO;
}

template<sp_d T> int FullMat<T>::direct_solve(Mat<T>& X, Mat<T>&& B) {
    if(this->factored) return this->solve_trs(X, std::forward<Mat<T>>(B));

    auto N = static_cast<int>(this->n_rows);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDB = static_cast<int>(B.n_rows);
    auto INFO = 0;

    this->pivot.zeros(N);

    this->factored = true;

    if(std::is_same_v<T, float>) {
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

    if(0 != INFO) SP_E("Error code {} received, the matrix is probably singular.\n", INFO);

    return INFO;
}

template<sp_d T> int FullMat<T>::solve_trs(Mat<T>& X, Mat<T>&& B) {
    const auto N = static_cast<int>(this->n_rows);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDB = static_cast<int>(B.n_rows);
    auto INFO = 0;

    if(std::is_same_v<T, float>) {
        using E = float;
        arma_fortran(arma_sgetrs)(&TRAN, &N, &NRHS, (E*)this->memptr(), &N, this->pivot.memptr(), (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else if(Precision::FULL == this->setting.precision) {
        using E = double;
        arma_fortran(arma_dgetrs)(&TRAN, &N, &NRHS, (E*)this->memptr(), &N, this->pivot.memptr(), (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else {
        X = arma::zeros(B.n_rows, B.n_cols);

        auto multiplier = arma::norm(B);

        auto counter = 0u;
        while(counter++ < this->setting.iterative_refinement) {
            if(multiplier < this->setting.tolerance) break;

            auto residual = conv_to<fmat>::from(B / multiplier);

            arma_fortran(arma_sgetrs)(&TRAN, &N, &NRHS, this->s_memory.memptr(), &N, this->pivot.memptr(), residual.memptr(), &LDB, &INFO);
            if(0 != INFO) break;

            const mat incre = multiplier * conv_to<mat>::from(residual);

            X += incre;

            suanpan_debug("mixed precision algorithm multiplier: %.5E.\n", multiplier = arma::norm(B -= this->operator*(incre)));
        }
    }

    return INFO;
}

#endif

//! @}
