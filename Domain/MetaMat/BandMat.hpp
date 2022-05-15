/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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

template<sp_d T> class BandMat final : public DenseMat<T> {
    static constexpr char TRAN = 'N';

    static T bin;

    const uword l_band;
    const uword u_band;
    const uword s_band;
    const uword m_rows; // memory block layout

    int solve_trs(Mat<T>&, Mat<T>&&);
    int solve_trs(Mat<T>&, const Mat<T>&);

public:
    BandMat(uword, uword, uword);

    unique_ptr<MetaMat<T>> make_copy() override;

    void unify(uword) override;
    void nullify(uword) override;

    const T& operator()(uword, uword) const override;
    T& at(uword, uword) override;

    Mat<T> operator*(const Mat<T>&) override;

    int solve(Mat<T>&, Mat<T>&&) override;
    int solve(Mat<T>&, const Mat<T>&) override;
};

template<sp_d T> T BandMat<T>::bin = 0.;

template<sp_d T> BandMat<T>::BandMat(const uword in_size, const uword in_l, const uword in_u)
    : DenseMat<T>(in_size, in_size, (2 * in_l + in_u + 1) * in_size)
    , l_band(in_l)
    , u_band(in_u)
    , s_band(in_l + in_u)
    , m_rows(2 * in_l + in_u + 1) {}

template<sp_d T> unique_ptr<MetaMat<T>> BandMat<T>::make_copy() { return std::make_unique<BandMat<T>>(*this); }

template<sp_d T> void BandMat<T>::unify(const uword K) {
    nullify(K);
    access::rw(this->memory[s_band + K * m_rows]) = 1.;
}

template<sp_d T> void BandMat<T>::nullify(const uword K) {
    suanpan_for(std::max(K, u_band) - u_band, std::min(this->n_rows, K + l_band + 1), [&](const uword I) { access::rw(this->memory[I + s_band + K * (m_rows - 1)]) = 0.; });
    suanpan_for(std::max(K, l_band) - l_band, std::min(this->n_cols, K + u_band + 1), [&](const uword I) { access::rw(this->memory[K + s_band + I * (m_rows - 1)]) = 0.; });

    this->factored = false;
}

template<sp_d T> const T& BandMat<T>::operator()(const uword in_row, const uword in_col) const {
    if(in_row > in_col + l_band || in_row + u_band < in_col) return bin = 0.;
    return this->memory[in_row + s_band + in_col * (m_rows - 1)];
}

template<sp_d T> T& BandMat<T>::at(const uword in_row, const uword in_col) {
    if(in_row > in_col + l_band || in_row + u_band < in_col) return bin = 0.;
    this->factored = false;
    return access::rw(this->memory[in_row + s_band + in_col * (m_rows - 1)]);
}

template<sp_d T> Mat<T> BandMat<T>::operator*(const Mat<T>& X) {
    Mat<T> Y(arma::size(X));

    const auto M = static_cast<int>(this->n_rows);
    const auto N = static_cast<int>(this->n_cols);
    const auto KL = static_cast<int>(l_band);
    const auto KU = static_cast<int>(u_band);
    const auto LDA = static_cast<int>(m_rows);
    const auto INC = 1;
    T ALPHA = 1.;
    T BETA = 0.;

    if(std::is_same_v<T, float>) {
        using E = float;
        suanpan_for(0llu, X.n_cols, [&](const uword I) { arma_fortran(arma_sgbmv)(&TRAN, &M, &N, &KL, &KU, (E*)&ALPHA, (E*)(this->memptr() + l_band), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
    }
    else if(std::is_same_v<T, double>) {
        using E = double;
        suanpan_for(0llu, X.n_cols, [&](const uword I) { arma_fortran(arma_dgbmv)(&TRAN, &M, &N, &KL, &KU, (E*)&ALPHA, (E*)(this->memptr() + l_band), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
    }

    return Y;
}

template<sp_d T> int BandMat<T>::solve(Mat<T>& X, const Mat<T>& B) {
    if(this->factored) return this->solve_trs(X, B);

    suanpan_debug([&] { if(this->n_rows != this->n_cols) throw invalid_argument("requires a square matrix"); });

    auto INFO = 0;

    auto N = static_cast<int>(this->n_rows);
    const auto KL = static_cast<int>(l_band);
    const auto KU = static_cast<int>(u_band);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDAB = static_cast<int>(m_rows);
    const auto LDB = static_cast<int>(B.n_rows);
    this->pivot.zeros(N);
    this->factored = true;

    if(std::is_same_v<T, float>) {
        using E = float;
        X = B;
        arma_fortran(arma_sgbsv)(&N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, this->pivot.memptr(), (E*)X.memptr(), &LDB, &INFO);
    }
    else if(Precision::FULL == this->precision) {
        using E = double;
        X = B;
        arma_fortran(arma_dgbsv)(&N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, this->pivot.memptr(), (E*)X.memptr(), &LDB, &INFO);
    }
    else {
        this->s_memory = this->to_float();
        arma_fortran(arma_sgbtrf)(&N, &N, &KL, &KU, this->s_memory.memptr(), &LDAB, this->pivot.memptr(), &INFO);
        if(0 == INFO) INFO = this->solve_trs(X, B);
    }

    if(0 != INFO) suanpan_error("solve() receives error code %u from base driver, the matrix is probably singular.\n", INFO);

    return INFO;
}

template<sp_d T> int BandMat<T>::solve_trs(Mat<T>& X, const Mat<T>& B) {
    auto INFO = 0;

    const auto N = static_cast<int>(this->n_rows);
    const auto KL = static_cast<int>(l_band);
    const auto KU = static_cast<int>(u_band);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDAB = static_cast<int>(m_rows);
    const auto LDB = static_cast<int>(B.n_rows);

    if(std::is_same_v<T, float>) {
        using E = float;
        X = B;
        arma_fortran(arma_sgbtrs)(&TRAN, &N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, this->pivot.memptr(), (E*)X.memptr(), &LDB, &INFO);
    }
    else if(Precision::FULL == this->precision) {
        using E = double;
        X = B;
        arma_fortran(arma_dgbtrs)(&TRAN, &N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, this->pivot.memptr(), (E*)X.memptr(), &LDB, &INFO);
    }
    else {
        X = arma::zeros(B.n_rows, B.n_cols);

        mat full_residual = B;

        auto multiplier = norm(full_residual);

        auto counter = 0u;
        while(counter++ < this->refinement) {
            if(multiplier < this->tolerance) break;

            auto residual = conv_to<fmat>::from(full_residual / multiplier);

            arma_fortran(arma_sgbtrs)(&TRAN, &N, &KL, &KU, &NRHS, this->s_memory.memptr(), &LDAB, this->pivot.memptr(), residual.memptr(), &LDB, &INFO);
            if(0 != INFO) break;

            const mat incre = multiplier * conv_to<mat>::from(residual);

            X += incre;

            suanpan_debug("mixed precision algorithm multiplier: %.5E.\n", multiplier = arma::norm(full_residual -= this->operator*(incre)));
        }
    }

    if(INFO != 0) suanpan_error("solve() receives error code %u from base driver, the matrix is probably singular.\n", INFO);

    return INFO;
}

template<sp_d T> int BandMat<T>::solve(Mat<T>& X, Mat<T>&& B) {
    if(this->factored) return this->solve_trs(X, std::forward<Mat<T>>(B));

    suanpan_debug([&] { if(this->n_rows != this->n_cols) throw invalid_argument("requires a square matrix"); });

    auto INFO = 0;

    auto N = static_cast<int>(this->n_rows);
    const auto KL = static_cast<int>(l_band);
    const auto KU = static_cast<int>(u_band);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDAB = static_cast<int>(m_rows);
    const auto LDB = static_cast<int>(B.n_rows);
    this->pivot.zeros(N);
    this->factored = true;

    if(std::is_same_v<T, float>) {
        using E = float;
        arma_fortran(arma_sgbsv)(&N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, this->pivot.memptr(), (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else if(Precision::FULL == this->precision) {
        using E = double;
        arma_fortran(arma_dgbsv)(&N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, this->pivot.memptr(), (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else {
        this->s_memory = this->to_float();
        arma_fortran(arma_sgbtrf)(&N, &N, &KL, &KU, this->s_memory.memptr(), &LDAB, this->pivot.memptr(), &INFO);
        if(0 == INFO) INFO = this->solve_trs(X, std::forward<Mat<T>>(B));
    }

    if(0 != INFO) suanpan_error("solve() receives error code %u from base driver, the matrix is probably singular.\n", INFO);

    return INFO;
}

template<sp_d T> int BandMat<T>::solve_trs(Mat<T>& X, Mat<T>&& B) {
    auto INFO = 0;

    const auto N = static_cast<int>(this->n_rows);
    const auto KL = static_cast<int>(l_band);
    const auto KU = static_cast<int>(u_band);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDAB = static_cast<int>(m_rows);
    const auto LDB = static_cast<int>(B.n_rows);

    if(std::is_same_v<T, float>) {
        using E = float;
        arma_fortran(arma_sgbtrs)(&TRAN, &N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, this->pivot.memptr(), (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else if(Precision::FULL == this->precision) {
        using E = double;
        arma_fortran(arma_dgbtrs)(&TRAN, &N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, this->pivot.memptr(), (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else {
        X = arma::zeros(B.n_rows, B.n_cols);

        auto multiplier = norm(B);

        auto counter = 0u;
        while(counter++ < this->refinement) {
            if(multiplier < this->tolerance) break;

            auto residual = conv_to<fmat>::from(B / multiplier);

            arma_fortran(arma_sgbtrs)(&TRAN, &N, &KL, &KU, &NRHS, this->s_memory.memptr(), &LDAB, this->pivot.memptr(), residual.memptr(), &LDB, &INFO);
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
