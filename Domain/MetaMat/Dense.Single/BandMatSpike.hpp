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
 * @class BandMatSpike
 * @brief A BandMatSpike class that holds matrices.
 *
 * @author tlc
 * @date 08/01/2020
 * @version 0.1.0
 * @file BandMatSpike.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef BANDMATSPIKE_HPP
#define BANDMATSPIKE_HPP

#include "../DenseMat.hpp"

extern "C" {
void spikeinit_(la_it*, la_it*, la_it*);
void dspike_tune_(la_it*);
void dspike_gbsv_(la_it*, la_it*, la_it*, la_it*, la_it*, double*, la_it*, double*, la_it*, la_it*);
void dspike_gbtrf_(la_it*, la_it*, la_it*, la_it*, double*, la_it*, double*, la_it*);
void dspike_gbtrs_(la_it*, const char*, la_it*, la_it*, la_it*, la_it*, double*, la_it*, double*, double*, la_it*);
void sspike_tune_(la_it*);
void sspike_gbsv_(la_it*, la_it*, la_it*, la_it*, la_it*, float*, la_it*, float*, la_it*, la_it*);
void sspike_gbtrf_(la_it*, la_it*, la_it*, la_it*, float*, la_it*, float*, la_it*);
void sspike_gbtrs_(la_it*, const char*, la_it*, la_it*, la_it*, la_it*, float*, la_it*, float*, float*, la_it*);
}

template<sp_d T> class BandMatSpike final : public DenseMat<T> {
    static constexpr char TRAN = 'N';

    static la_it SPROTO, DPROTO;

    static T bin;

    const uword l_band;
    const uword u_band;
    const uword m_rows; // memory block layout

    la_it SPIKE[64]{};

    podarray<T> WORK;
    podarray<float> SWORK;

    void init_spike() {
        auto N = static_cast<la_it>(this->n_rows);
        auto KLU = static_cast<la_it>(std::max(l_band, u_band));

        spikeinit_(SPIKE, &N, &KLU);

        SPIKE[6] = std::is_same_v<T, float> ? SPROTO : DPROTO;
        SPIKE[4] = SPIKE[6] + SPIKE[6] / 2 + 10;
        SPIKE[3] = SPIKE[4] / 2;
    }

    int solve_trs(Mat<T>&, Mat<T>&&);

protected:
    using DenseMat<T>::direct_solve;

    int direct_solve(Mat<T>&, Mat<T>&&) override;

public:
    BandMatSpike(const uword in_size, const uword in_l, const uword in_u)
        : DenseMat<T>(in_size, in_size, (in_l + in_u + 1) * in_size)
        , l_band(in_l)
        , u_band(in_u)
        , m_rows(in_l + in_u + 1) { init_spike(); }

    BandMatSpike(const BandMatSpike& other)
        : DenseMat<T>(other)
        , l_band(other.l_band)
        , u_band(other.u_band)
        , m_rows(other.m_rows) { init_spike(); }

    BandMatSpike(BandMatSpike&&) noexcept = delete;
    BandMatSpike& operator=(const BandMatSpike&) = delete;
    BandMatSpike& operator=(BandMatSpike&&) noexcept = delete;

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<BandMatSpike>(*this); }

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

    [[nodiscard]] int sign_det() const override { throw invalid_argument("not supported"); }
};

template<sp_d T> la_it BandMatSpike<T>::SPROTO = [] {
    la_it PROTO[64]{};
    sspike_tune_(PROTO);
    return PROTO[6];
}();

template<sp_d T> la_it BandMatSpike<T>::DPROTO = [] {
    la_it PROTO[64]{};
    dspike_tune_(PROTO);
    return PROTO[6];
}();

template<sp_d T> T BandMatSpike<T>::bin = T(0);

template<sp_d T> Mat<T> BandMatSpike<T>::operator*(const Mat<T>& X) const {
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
        suanpan::for_each(X.n_cols, [&](const uword I) { arma_fortran(arma_sgbmv)(&TRAN, &M, &N, &KL, &KU, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
    }
    else {
        using E = double;
        suanpan::for_each(X.n_cols, [&](const uword I) { arma_fortran(arma_dgbmv)(&TRAN, &M, &N, &KL, &KU, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
    }

    return Y;
}

template<sp_d T> int BandMatSpike<T>::direct_solve(Mat<T>& X, Mat<T>&& B) {
    if(!this->factored) {
        suanpan_assert([&] { if(this->n_rows != this->n_cols) throw invalid_argument("requires a square matrix"); });

        la_it INFO = 0;

        auto N = static_cast<la_it>(this->n_rows);
        auto KL = static_cast<la_it>(l_band);
        auto KU = static_cast<la_it>(u_band);
        auto LDAB = static_cast<la_it>(m_rows);
        const auto KLU = std::max(l_band, u_band);
        this->factored = true;

        if constexpr(std::is_same_v<T, float>) {
            using E = float;
            WORK.zeros(KLU * KLU * SPIKE[9]);
            sspike_gbtrf_(SPIKE, &N, &KL, &KU, (E*)this->memptr(), &LDAB, (E*)WORK.memptr(), &INFO);
        }
        else if(Precision::FULL == this->setting.precision) {
            using E = double;
            WORK.zeros(KLU * KLU * SPIKE[9]);
            dspike_gbtrf_(SPIKE, &N, &KL, &KU, (E*)this->memptr(), &LDAB, (E*)WORK.memptr(), &INFO);
        }
        else {
            this->s_memory = this->to_float();
            SWORK.zeros(KLU * KLU * SPIKE[9]);
            sspike_gbtrf_(SPIKE, &N, &KL, &KU, this->s_memory.memptr(), &LDAB, SWORK.memptr(), &INFO);
        }

        if(0 != INFO) {
            suanpan_error("Error code {} received, the matrix is probably singular.\n", INFO);
            return INFO;
        }
    }

    return this->solve_trs(X, std::move(B));
}

template<sp_d T> int BandMatSpike<T>::solve_trs(Mat<T>& X, Mat<T>&& B) {
    auto N = static_cast<la_it>(this->n_rows);
    auto KL = static_cast<la_it>(l_band);
    auto KU = static_cast<la_it>(u_band);
    auto NRHS = static_cast<la_it>(B.n_cols);
    auto LDAB = static_cast<la_it>(m_rows);
    auto LDB = static_cast<la_it>(B.n_rows);

    if constexpr(std::is_same_v<T, float>) {
        using E = float;
        sspike_gbtrs_(SPIKE, &TRAN, &N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, (E*)WORK.memptr(), (E*)B.memptr(), &LDB);
        X = std::move(B);
    }
    else if(Precision::FULL == this->setting.precision) {
        using E = double;
        dspike_gbtrs_(SPIKE, &TRAN, &N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, (E*)WORK.memptr(), (E*)B.memptr(), &LDB);
        X = std::move(B);
    }
    else
        this->mixed_trs(X, std::move(B), [&](fmat& residual) {
            sspike_gbtrs_(SPIKE, &TRAN, &N, &KL, &KU, &NRHS, this->s_memory.memptr(), &LDAB, SWORK.memptr(), residual.memptr(), &LDB);
            return 0;
        });

    return SUANPAN_SUCCESS;
}

#endif

//! @}
