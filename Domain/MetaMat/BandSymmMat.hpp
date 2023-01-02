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
    int solve_trs(Mat<T>&, const Mat<T>&);

public:
    BandSymmMat(uword, uword);

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

template<sp_d T> T BandSymmMat<T>::bin = 0.;

template<sp_d T> BandSymmMat<T>::BandSymmMat(const uword in_size, const uword in_bandwidth)
    : DenseMat<T>(in_size, in_size, (in_bandwidth + 1) * in_size)
    , band(in_bandwidth)
    , m_rows(in_bandwidth + 1) {}

template<sp_d T> unique_ptr<MetaMat<T>> BandSymmMat<T>::make_copy() { return std::make_unique<BandSymmMat<T>>(*this); }

template<sp_d T> void BandSymmMat<T>::unify(const uword K) {
    nullify(K);
    access::rw(this->memory[K * m_rows]) = 1.;
}

template<sp_d T> void BandSymmMat<T>::nullify(const uword K) {
    suanpan_for(std::max(band, K) - band, K, [&](const uword I) { access::rw(this->memory[K - I + I * m_rows]) = 0.; });
    const auto t_factor = K * m_rows - K;
    suanpan_for(K, std::min(this->n_rows, K + band + 1), [&](const uword I) { access::rw(this->memory[I + t_factor]) = 0.; });

    this->factored = false;
}

template<sp_d T> const T& BandSymmMat<T>::operator()(const uword in_row, const uword in_col) const {
    if(in_row > band + in_col) return bin = 0.;
    return this->memory[in_row > in_col ? in_row - in_col + in_col * m_rows : in_col - in_row + in_row * m_rows];
}

template<sp_d T> T& BandSymmMat<T>::unsafe_at(const uword in_row, const uword in_col) {
    this->factored = false;
    return access::rw(this->memory[in_row - in_col + in_col * m_rows]);
}

template<sp_d T> T& BandSymmMat<T>::at(const uword in_row, const uword in_col) {
    if(in_row > band + in_col || in_row < in_col) [[unlikely]] return bin = 0.;
    this->factored = false;
    return access::rw(this->memory[in_row - in_col + in_col * m_rows]);
}

template<sp_d T> Mat<T> BandSymmMat<T>::operator*(const Mat<T>& X) const {
    Mat<T> Y(arma::size(X));

    const auto N = static_cast<int>(this->n_cols);
    const auto K = static_cast<int>(band);
    const auto LDA = static_cast<int>(m_rows);
    const auto INC = 1;
    T ALPHA = 1.;
    T BETA = 0.;

    if(std::is_same_v<T, float>) {
        using E = float;
        suanpan_for(0llu, X.n_cols, [&](const uword I) { arma_fortran(arma_ssbmv)(&UPLO, &N, &K, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
    }
    else if(std::is_same_v<T, double>) {
        using E = double;
        suanpan_for(0llu, X.n_cols, [&](const uword I) { arma_fortran(arma_dsbmv)(&UPLO, &N, &K, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
    }

    return Y;
}

template<sp_d T> int BandSymmMat<T>::direct_solve(Mat<T>& X, const Mat<T>& B) {
    if(this->factored) return this->solve_trs(X, B);

    const auto N = static_cast<int>(this->n_rows);
    const auto KD = static_cast<int>(band);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDAB = static_cast<int>(m_rows);
    const auto LDB = static_cast<int>(B.n_rows);
    auto INFO = 0;

    this->factored = true;

    if(std::is_same_v<T, float>) {
        using E = float;
        X = B;
        arma_fortran(arma_spbsv)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)X.memptr(), &LDB, &INFO);
    }
    else if(Precision::FULL == this->setting.precision) {
        using E = double;
        X = B;
        arma_fortran(arma_dpbsv)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)X.memptr(), &LDB, &INFO);
    }
    else {
        this->s_memory = this->to_float();
        arma_fortran(arma_spbtrf)(&UPLO, &N, &KD, this->s_memory.memptr(), &LDAB, &INFO);
        if(0 == INFO) INFO = this->solve_trs(X, B);
    }

    if(0 != INFO)
        SP_E("Error code {} received, the matrix is probably singular.\n", INFO);

    return INFO;
}

template<sp_d T> int BandSymmMat<T>::solve_trs(Mat<T>& X, const Mat<T>& B) {
    const auto N = static_cast<int>(this->n_rows);
    const auto KD = static_cast<int>(band);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDAB = static_cast<int>(m_rows);
    const auto LDB = static_cast<int>(B.n_rows);
    auto INFO = 0;

    if(std::is_same_v<T, float>) {
        using E = float;
        X = B;
        arma_fortran(arma_spbtrs)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)X.memptr(), &LDB, &INFO);
    }
    else if(Precision::FULL == this->setting.precision) {
        using E = double;
        X = B;
        arma_fortran(arma_dpbtrs)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)X.memptr(), &LDB, &INFO);
    }
    else {
        X = arma::zeros(B.n_rows, B.n_cols);

        mat full_residual = B;

        auto multiplier = norm(full_residual);

        auto counter = 0u;
        while(counter++ < this->setting.iterative_refinement) {
            if(multiplier < this->setting.tolerance) break;

            auto residual = conv_to<fmat>::from(full_residual / multiplier);

            arma_fortran(arma_spbtrs)(&UPLO, &N, &KD, &NRHS, this->s_memory.memptr(), &LDAB, residual.memptr(), &LDB, &INFO);
            if(0 != INFO) break;

            const mat incre = multiplier * conv_to<mat>::from(residual);

            X += incre;

            suanpan_debug("mixed precision algorithm multiplier: %.5E.\n", multiplier = arma::norm(full_residual -= this->operator*(incre)));
        }
    }

    if(0 != INFO)
        SP_E("Error code {} received, the matrix is probably singular.\n", INFO);

    return INFO;
}

template<sp_d T> int BandSymmMat<T>::direct_solve(Mat<T>& X, Mat<T>&& B) {
    if(this->factored) return this->solve_trs(X, std::forward<Mat<T>>(B));

    const auto N = static_cast<int>(this->n_rows);
    const auto KD = static_cast<int>(band);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDAB = static_cast<int>(m_rows);
    const auto LDB = static_cast<int>(B.n_rows);
    auto INFO = 0;

    this->factored = true;

    if(std::is_same_v<T, float>) {
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
        SP_E("Error code {} received, the matrix is probably singular.\n", INFO);

    return INFO;
}

template<sp_d T> int BandSymmMat<T>::solve_trs(Mat<T>& X, Mat<T>&& B) {
    const auto N = static_cast<int>(this->n_rows);
    const auto KD = static_cast<int>(band);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDAB = static_cast<int>(m_rows);
    const auto LDB = static_cast<int>(B.n_rows);
    auto INFO = 0;

    if(std::is_same_v<T, float>) {
        using E = float;
        arma_fortran(arma_spbtrs)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else if(Precision::FULL == this->setting.precision) {
        using E = double;
        arma_fortran(arma_dpbtrs)(&UPLO, &N, &KD, &NRHS, (E*)this->memptr(), &LDAB, (E*)B.memptr(), &LDB, &INFO);
        X = std::move(B);
    }
    else {
        X = arma::zeros(B.n_rows, B.n_cols);

        auto multiplier = norm(B);

        auto counter = 0u;
        while(counter++ < this->setting.iterative_refinement) {
            if(multiplier < this->setting.tolerance) break;

            auto residual = conv_to<fmat>::from(B / multiplier);

            arma_fortran(arma_spbtrs)(&UPLO, &N, &KD, &NRHS, this->s_memory.memptr(), &LDAB, residual.memptr(), &LDB, &INFO);
            if(0 != INFO) break;

            const mat incre = multiplier * conv_to<mat>::from(residual);

            X += incre;

            suanpan_debug("mixed precision algorithm multiplier: %.5E.\n", multiplier = arma::norm(B -= this->operator*(incre)));
        }
    }

    if(0 != INFO)
        SP_E("Error code {} received, the matrix is probably singular.\n", INFO);

    return INFO;
}

#endif

//! @}
