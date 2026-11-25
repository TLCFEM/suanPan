/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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
 * @class BandMatSpikeBase
 * @brief A BandMatSpikeBase class that holds matrices.
 *
 * The `BandMatSpikeBase` class stores a matrix in the banded layout shared by
 * the CPU and GPU based SPIKE solvers.
 *
 * @author tlc
 * @date 25/11/2026
 * @version 0.1.0
 * @file BandMatSpikeBase.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef BANDMATSPIKEBASE_HPP
#define BANDMATSPIKEBASE_HPP

#include "../DenseMat.hpp"

template<sp_d T> class BandMatSpikeBase : public DenseMat<T> {
    static constexpr auto TRAN = 'N';

    static T bin;

protected:
    const uword l_band;
    const uword u_band;
    const uword m_rows; // memory block layout

    using DenseMat<T>::direct_solve;

public:
    BandMatSpikeBase(const uword in_size, const uword in_l, const uword in_u)
        : DenseMat<T>(in_size, in_size, (in_l + in_u + 1) * in_size)
        , l_band(in_l)
        , u_band(in_u)
        , m_rows(in_l + in_u + 1) {}

    BandMatSpikeBase(const BandMatSpikeBase& other)
        : DenseMat<T>(other)
        , l_band(other.l_band)
        , u_band(other.u_band)
        , m_rows(other.m_rows) {}

    BandMatSpikeBase(BandMatSpikeBase&&) = delete;
    BandMatSpikeBase& operator=(const BandMatSpikeBase&) = delete;
    BandMatSpikeBase& operator=(BandMatSpikeBase&&) = delete;

    void nullify(const uword K) override {
        this->factored = false;
        suanpan::for_each(std::max(K, u_band) - u_band, std::min(this->n_rows, K + l_band + 1), [&](const uword I) { this->memory[I + u_band + K * (m_rows - 1)] = T(0); });
        suanpan::for_each(std::max(K, l_band) - l_band, std::min(this->n_cols, K + u_band + 1), [&](const uword I) { this->memory[K + u_band + I * (m_rows - 1)] = T(0); });
    }

    T operator()(const uword in_row, const uword in_col) const override {
        if(in_row > in_col + l_band || in_row + u_band < in_col) [[unlikely]]
            return bin = T(0);
        return this->memory[in_row + u_band + in_col * (m_rows - 1)];
    }

    T& unsafe_at(const uword in_row, const uword in_col) override {
        this->factored = false;
        return this->memory[in_row + u_band + in_col * (m_rows - 1)];
    }

    T& at(const uword in_row, const uword in_col) override {
        if(in_row > in_col + l_band || in_row + u_band < in_col) [[unlikely]]
            return bin = T(0);
        return this->unsafe_at(in_row, in_col);
    }

    Mat<T> operator*(const Mat<T>&) const override;
};

template<sp_d T> T BandMatSpikeBase<T>::bin = T(0);

template<sp_d T> Mat<T> BandMatSpikeBase<T>::operator*(const Mat<T>& X) const {
    Mat<T> Y(arma::size(X));

    const auto M = static_cast<blas_int>(this->n_rows);
    const auto N = static_cast<blas_int>(this->n_cols);
    const auto KL = static_cast<blas_int>(l_band);
    const auto KU = static_cast<blas_int>(u_band);
    const auto LDA = static_cast<blas_int>(m_rows);
    static constexpr blas_int INC = 1;
    static constexpr T ALPHA{1}, BETA{0};

    if constexpr(std::is_same_v<T, float>) suanpan::for_each(X.n_cols, [&](const uword I) { arma_fortran(arma_sgbmv)(&TRAN, &M, &N, &KL, &KU, &ALPHA, this->memptr(), &LDA, X.colptr(I), &INC, &BETA, Y.colptr(I), &INC); });
    else suanpan::for_each(X.n_cols, [&](const uword I) { arma_fortran(arma_dgbmv)(&TRAN, &M, &N, &KL, &KU, &ALPHA, this->memptr(), &LDA, X.colptr(I), &INC, &BETA, Y.colptr(I), &INC); });
    return Y;
}

#endif

//! @}
