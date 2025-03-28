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
 * @class BandMat
 * @brief A BandMat class that holds matrices.
 *
 * @author tlc
 * @date 06/09/2017
 * @version 0.1.0
 * @file BandMat.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef BANDMAT_HPP
#define BANDMAT_HPP

#include "DenseMat.hpp"

template<sp_d T> class BandMat : public DenseMat<T> {
    static constexpr char TRAN = 'N';

    static T bin;

    const uword s_band;

    int solve_trs(Mat<T>&, Mat<T>&&);

protected:
    const uword m_rows; // memory block layout

    const uword l_band;
    const uword u_band;

    using DenseMat<T>::direct_solve;

    int direct_solve(Mat<T>&, Mat<T>&&) override;

public:
    BandMat(const uword in_size, const uword in_l, const uword in_u)
        : DenseMat<T>(in_size, in_size, (2 * in_l + in_u + 1) * in_size)
        , s_band(in_l + in_u)
        , m_rows(2 * in_l + in_u + 1)
        , l_band(in_l)
        , u_band(in_u) {
        if(m_rows >= in_size)
            suanpan_warning("The storage requirement for the banded matrix is larger than that of a full matrix, consider using a full/sparse matrix instead.\n");
    }

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<BandMat>(*this); }

    void nullify(const uword K) override {
        this->factored = false;
        suanpan::for_each(std::max(K, u_band) - u_band, std::min(this->n_rows, K + l_band + 1), [&](const uword I) { this->memory[I + s_band + K * (m_rows - 1)] = T(0); });
        suanpan::for_each(std::max(K, l_band) - l_band, std::min(this->n_cols, K + u_band + 1), [&](const uword I) { this->memory[K + s_band + I * (m_rows - 1)] = T(0); });
    }

    T operator()(const uword in_row, const uword in_col) const override {
        if(in_row > in_col + l_band || in_row + u_band < in_col) [[unlikely]] return bin = T(0);
        return this->memory[in_row + s_band + in_col * (m_rows - 1)];
    }

    T& unsafe_at(const uword in_row, const uword in_col) override {
        this->factored = false;
        return this->memory[in_row + s_band + in_col * (m_rows - 1)];
    }

    T& at(const uword in_row, const uword in_col) override {
        if(in_row > in_col + l_band || in_row + u_band < in_col) [[unlikely]] return bin = T(0);
        return this->unsafe_at(in_row, in_col);
    }

    Mat<T> operator*(const Mat<T>&) const override;
};

template<sp_d T> T BandMat<T>::bin = T(0);

template<sp_d T> Mat<T> BandMat<T>::operator*(const Mat<T>& X) const {
    Mat<T> Y(arma::size(X));

    const auto M = static_cast<blas_int>(this->n_rows);
    const auto N = static_cast<blas_int>(this->n_cols);
    const auto KL = static_cast<blas_int>(l_band);
    const auto KU = static_cast<blas_int>(u_band);
    const auto LDA = static_cast<blas_int>(m_rows);
    constexpr blas_int INC = 1;
    T ALPHA = T(1);
    T BETA = T(0);

    if constexpr(std::is_same_v<T, float>) {
        using E = float;
        suanpan::for_each(X.n_cols, [&](const uword I) { arma_fortran(arma_sgbmv)(&TRAN, &M, &N, &KL, &KU, (E*)&ALPHA, (E*)(this->memptr() + l_band), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
    }
    else {
        using E = double;
        suanpan::for_each(X.n_cols, [&](const uword I) { arma_fortran(arma_dgbmv)(&TRAN, &M, &N, &KL, &KU, (E*)&ALPHA, (E*)(this->memptr() + l_band), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
    }

    return Y;
}

template<sp_d T> int BandMat<T>::direct_solve(Mat<T>& X, Mat<T>&& B) {
    if(this->factored) return this->solve_trs(X, std::forward<Mat<T>>(B));

    suanpan_assert([&] { if(this->n_rows != this->n_cols) throw invalid_argument("requires a square matrix"); });

    blas_int INFO = 0;

    auto N = static_cast<blas_int>(this->n_rows);
    const auto KL = static_cast<blas_int>(l_band);
    const auto KU = static_cast<blas_int>(u_band);
    const auto NRHS = static_cast<blas_int>(B.n_cols);
    const auto LDAB = static_cast<blas_int>(m_rows);
    const auto LDB = static_cast<blas_int>(B.n_rows);
    this->pivot.zeros(N);
    this->factored = true;

    if constexpr(std::is_same_v<T, float>) {
        using E = float;
        arma_fortran(arma_sgbsv)(&N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, this->pivot.memptr(), (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else if(Precision::FULL == this->setting.precision) {
        using E = double;
        arma_fortran(arma_dgbsv)(&N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, this->pivot.memptr(), (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else {
        this->s_memory = this->to_float();
        arma_fortran(arma_sgbtrf)(&N, &N, &KL, &KU, this->s_memory.memptr(), &LDAB, this->pivot.memptr(), &INFO);
        if(0 == INFO) INFO = this->solve_trs(X, std::forward<Mat<T>>(B));
    }

    if(0 != INFO)
        suanpan_error("Error code {} received, the matrix is probably singular.\n", INFO);

    return INFO;
}

template<sp_d T> int BandMat<T>::solve_trs(Mat<T>& X, Mat<T>&& B) {
    blas_int INFO = 0;

    const auto N = static_cast<blas_int>(this->n_rows);
    const auto KL = static_cast<blas_int>(l_band);
    const auto KU = static_cast<blas_int>(u_band);
    const auto NRHS = static_cast<blas_int>(B.n_cols);
    const auto LDAB = static_cast<blas_int>(m_rows);
    const auto LDB = static_cast<blas_int>(B.n_rows);

    if constexpr(std::is_same_v<T, float>) {
        using E = float;
        arma_fortran(arma_sgbtrs)(&TRAN, &N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, this->pivot.memptr(), (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else if(Precision::FULL == this->setting.precision) {
        using E = double;
        arma_fortran(arma_dgbtrs)(&TRAN, &N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, this->pivot.memptr(), (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else
        this->mixed_trs(X, std::forward<Mat<T>>(B), [&](fmat& residual) {
            arma_fortran(arma_sgbtrs)(&TRAN, &N, &KL, &KU, &NRHS, this->s_memory.memptr(), &LDAB, this->pivot.memptr(), residual.memptr(), &LDB, &INFO);
            return INFO;
        });

    if(0 != INFO)
        suanpan_error("Error code {} received, the matrix is probably singular.\n", INFO);

    return INFO;
}

#endif

//! @}
