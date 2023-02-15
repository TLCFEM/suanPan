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
 * @class BandSymmMat
 * @brief A BandSymmMat class that holds matrices.
 *
 * @author tlc
 * @date 06/09/2017
 * @version 0.1.0
 * @file BandSymmMat.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef BANDSYMMMAT_HPP
#define BANDSYMMMAT_HPP

#include "DenseMat.hpp"

template<sp_d T> class BandSymmMat final : public DenseMat<T> {
    static constexpr char UPLO = 'L';

    static T bin;

    const uword band;
    const uword m_rows; // memory block layout

    int solve_trs(Mat<T>&, Mat<T>&&);

protected:
    using MetaMat<T>::direct_solve;

    int direct_solve(Mat<T>&, Mat<T>&&) override;

public:
    BandSymmMat(const uword in_size, const uword in_bandwidth)
        : DenseMat<T>(in_size, in_size, (in_bandwidth + 1) * in_size)
        , band(in_bandwidth)
        , m_rows(in_bandwidth + 1) {}

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<BandSymmMat>(*this); }

    void nullify(const uword K) override {
        this->factored = false;
        suanpan_for(std::max(band, K) - band, K, [&](const uword I) { this->memory[K - I + I * m_rows] = T(0); });
        const auto t_factor = K * m_rows - K;
        suanpan_for(K, std::min(this->n_rows, K + band + 1), [&](const uword I) { this->memory[I + t_factor] = T(0); });
    }

    T operator()(const uword in_row, const uword in_col) const override {
        if(in_row > band + in_col || in_col > in_row + band) [[unlikely]] return bin = T(0);
        return this->memory[in_row > in_col ? in_row - in_col + in_col * m_rows : in_col - in_row + in_row * m_rows];
    }

    T& unsafe_at(const uword in_row, const uword in_col) override {
        this->factored = false;
        return this->memory[in_row - in_col + in_col * m_rows];
    }

    T& at(const uword in_row, const uword in_col) override {
        if(in_row > band + in_col || in_row < in_col) [[unlikely]] return bin = T(0);
        return this->unsafe_at(in_row, in_col);
    }

    Mat<T> operator*(const Mat<T>&) const override;
};

template<sp_d T> T BandSymmMat<T>::bin = T(0);

template<sp_d T> Mat<T> BandSymmMat<T>::operator*(const Mat<T>& X) const {
    Mat<T> Y(arma::size(X));

    const auto N = static_cast<int>(this->n_cols);
    const auto K = static_cast<int>(band);
    const auto LDA = static_cast<int>(m_rows);
    constexpr auto INC = 1;
    T ALPHA = T(1);
    T BETA = T(0);

    if constexpr(std::is_same_v<T, float>) {
        using E = float;
        suanpan_for(0llu, X.n_cols, [&](const uword I) { arma_fortran(arma_ssbmv)(&UPLO, &N, &K, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
    }
    else {
        using E = double;
        suanpan_for(0llu, X.n_cols, [&](const uword I) { arma_fortran(arma_dsbmv)(&UPLO, &N, &K, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
    }

    return Y;
}

template<sp_d T> int BandSymmMat<T>::direct_solve(Mat<T>& X, Mat<T>&& B) {
    if(this->factored) return this->solve_trs(X, std::forward<Mat<T>>(B));

    suanpan_assert([&] { if(this->n_rows != this->n_cols) throw invalid_argument("requires a square matrix"); });

    auto INFO = 0;

    const auto N = static_cast<int>(this->n_rows);
    const auto KD = static_cast<int>(band);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDAB = static_cast<int>(m_rows);
    const auto LDB = static_cast<int>(B.n_rows);
    this->factored = true;

    if constexpr(std::is_same_v<T, float>) {
        using E = float;
        arma_fortran(arma_spbsv)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else if(Precision::FULL == this->setting.precision) {
        using E = double;
        arma_fortran(arma_dpbsv)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else {
        this->s_memory = this->to_float();
        arma_fortran(arma_spbtrf)(&UPLO, &N, &KD, this->s_memory.memptr(), &LDAB, &INFO);
        if(0 == INFO) INFO = this->solve_trs(X, std::forward<Mat<T>>(B));
    }

    if(0 != INFO)
        suanpan_error("Error code {} received, the matrix is probably singular.\n", INFO);

    return INFO;
}

template<sp_d T> int BandSymmMat<T>::solve_trs(Mat<T>& X, Mat<T>&& B) {
    auto INFO = 0;

    const auto N = static_cast<int>(this->n_rows);
    const auto KD = static_cast<int>(band);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDAB = static_cast<int>(m_rows);
    const auto LDB = static_cast<int>(B.n_rows);

    if constexpr(std::is_same_v<T, float>) {
        using E = float;
        arma_fortran(arma_spbtrs)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else if(Precision::FULL == this->setting.precision) {
        using E = double;
        arma_fortran(arma_dpbtrs)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else
        this->mixed_trs(X, std::forward<Mat<T>>(B), [&](fmat& residual) {
            arma_fortran(arma_spbtrs)(&UPLO, &N, &KD, &NRHS, this->s_memory.memptr(), &LDAB, residual.memptr(), &LDB, &INFO);
            return INFO;
        });

    if(0 != INFO)
        suanpan_error("Error code {} received, the matrix is probably singular.\n", INFO);

    return INFO;
}

#endif

//! @}
