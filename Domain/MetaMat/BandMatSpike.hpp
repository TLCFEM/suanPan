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

#include <feast/spike.h>
#include "DenseMat.hpp"

template<sp_d T> class BandMatSpike final : public DenseMat<T> {
    static constexpr char TRAN = 'N';

    static T bin;

    const uword l_band;
    const uword u_band;
    const uword m_rows; // memory block layout

    podarray<int> SPIKE = podarray<int>(64);
    podarray<T> WORK;
    podarray<float> SWORK;

    void init_spike();

    int solve_trs(Mat<T>&, Mat<T>&&);
    int solve_trs(Mat<T>&, const Mat<T>&);

public:
    BandMatSpike(uword, uword, uword);

    unique_ptr<MetaMat<T>> make_copy() override;

    void unify(uword) override;
    void nullify(uword) override;

    const T& operator()(uword, uword) const override;
    T& at(uword, uword) override;

    Mat<T> operator*(const Mat<T>&) const override;

    int solve(Mat<T>&, Mat<T>&&) override;
    int solve(Mat<T>&, const Mat<T>&) override;

    [[nodiscard]] int sign_det() const override;
};

template<sp_d T> T BandMatSpike<T>::bin = 0.;

template<sp_d T> void BandMatSpike<T>::init_spike() {
    auto N = static_cast<int>(this->n_rows);
    auto KLU = static_cast<int>(std::max(l_band, u_band));

    spikeinit_(SPIKE.memptr(), &N, &KLU);

    std::is_same_v<T, float> ? sspike_tune_(SPIKE.memptr()) : dspike_tune_(SPIKE.memptr());
}

template<sp_d T> BandMatSpike<T>::BandMatSpike(const uword in_size, const uword in_l, const uword in_u)
    : DenseMat<T>(in_size, in_size, (in_l + in_u + 1) * in_size)
    , l_band(in_l)
    , u_band(in_u)
    , m_rows(in_l + in_u + 1) { init_spike(); }

template<sp_d T> unique_ptr<MetaMat<T>> BandMatSpike<T>::make_copy() { return std::make_unique<BandMatSpike<T>>(*this); }

template<sp_d T> void BandMatSpike<T>::unify(const uword K) {
    nullify(K);
    access::rw(this->memory[u_band + K * m_rows]) = 1.;
}

template<sp_d T> void BandMatSpike<T>::nullify(const uword K) {
    suanpan_for(std::max(K, u_band) - u_band, std::min(this->n_rows, K + l_band + 1), [&](const uword I) { access::rw(this->memory[I + u_band + K * (m_rows - 1)]) = 0.; });
    suanpan_for(std::max(K, l_band) - l_band, std::min(this->n_cols, K + u_band + 1), [&](const uword I) { access::rw(this->memory[K + u_band + I * (m_rows - 1)]) = 0.; });

    this->factored = false;
}

template<sp_d T> const T& BandMatSpike<T>::operator()(const uword in_row, const uword in_col) const {
    if(in_row > in_col + l_band || in_row + u_band < in_col) return bin = 0.;
    return this->memory[in_row + u_band + in_col * (m_rows - 1)];
}

template<sp_d T> T& BandMatSpike<T>::at(const uword in_row, const uword in_col) {
    if(in_row > in_col + l_band || in_row + u_band < in_col) return bin = 0.;
    this->factored = false;
    return access::rw(this->memory[in_row + u_band + in_col * (m_rows - 1)]);
}

template<sp_d T> Mat<T> BandMatSpike<T>::operator*(const Mat<T>& X) const {
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
        suanpan_for(0llu, X.n_cols, [&](const uword I) { arma_fortran(arma_sgbmv)(&TRAN, &M, &N, &KL, &KU, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
    }
    else if(std::is_same_v<T, double>) {
        using E = double;
        suanpan_for(0llu, X.n_cols, [&](const uword I) { arma_fortran(arma_dgbmv)(&TRAN, &M, &N, &KL, &KU, (E*)&ALPHA, (E*)this->memptr(), &LDA, (E*)X.colptr(I), &INC, (E*)&BETA, (E*)Y.colptr(I), &INC); });
    }

    return Y;
}

template<sp_d T> int BandMatSpike<T>::solve(Mat<T>& X, const Mat<T>& B) {
    if(!this->factored) {
        auto N = static_cast<int>(this->n_rows);
        auto KL = static_cast<int>(l_band);
        auto KU = static_cast<int>(u_band);
        auto LDAB = static_cast<int>(m_rows);
        const auto KLU = std::max(l_band, u_band);
        auto INFO = 0;

        if(std::is_same_v<T, float>) {
            using E = float;
            WORK.zeros(KLU * KLU * SPIKE(9));
            sspike_gbtrf_(SPIKE.memptr(), &N, &KL, &KU, (E*)this->memptr(), &LDAB, (E*)WORK.memptr(), &INFO);
        }
        else if(Precision::FULL == this->setting.precision) {
            using E = double;
            WORK.zeros(KLU * KLU * SPIKE(9));
            dspike_gbtrf_(SPIKE.memptr(), &N, &KL, &KU, (E*)this->memptr(), &LDAB, (E*)WORK.memptr(), &INFO);
        }
        else {
            this->s_memory = this->to_float();
            SWORK.zeros(KLU * KLU * SPIKE(9));
            sspike_gbtrf_(SPIKE.memptr(), &N, &KL, &KU, this->s_memory.mem, &LDAB, SWORK.memptr(), &INFO);
        }

        if(INFO != 0) {
            suanpan_error("solve() receives error code %u from the base driver, the matrix is probably singular.\n", INFO);
            return INFO;
        }

        this->factored = true;
    }

    return this->solve_trs(X, B);
}

template<sp_d T> int BandMatSpike<T>::solve_trs(Mat<T>& X, const Mat<T>& B) {
    auto N = static_cast<int>(this->n_rows);
    auto KL = static_cast<int>(l_band);
    auto KU = static_cast<int>(u_band);
    auto NRHS = static_cast<int>(B.n_cols);
    auto LDAB = static_cast<int>(m_rows);
    auto LDB = static_cast<int>(B.n_rows);

    if(std::is_same_v<T, float>) {
        using E = float;
        X = B;
        sspike_gbtrs_(SPIKE.memptr(), &TRAN, &N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, (E*)WORK.memptr(), (E*)X.memptr(), &LDB);
    }
    else if(Precision::FULL == this->setting.precision) {
        using E = double;
        X = B;
        dspike_gbtrs_(SPIKE.memptr(), &TRAN, &N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, (E*)WORK.memptr(), (E*)X.memptr(), &LDB);
    }
    else {
        X = arma::zeros(B.n_rows, B.n_cols);

        mat full_residual = B;

        auto multiplier = norm(full_residual);

        auto counter = 0u;
        while(counter++ < this->setting.iterative_refinement) {
            if(multiplier < this->setting.tolerance) break;

            auto residual = conv_to<fmat>::from(full_residual / multiplier);

            sspike_gbtrs_(SPIKE.memptr(), &TRAN, &N, &KL, &KU, &NRHS, this->s_memory.memptr(), &LDAB, SWORK.memptr(), residual.memptr(), &LDB);

            const mat incre = multiplier * conv_to<mat>::from(residual);

            X += incre;

            suanpan_debug("mixed precision algorithm multiplier: %.5E.\n", multiplier = arma::norm(full_residual -= this->operator*(incre)));
        }
    }

    return SUANPAN_SUCCESS;
}

template<sp_d T> int BandMatSpike<T>::solve(Mat<T>& X, Mat<T>&& B) {
    if(!this->factored) {
        auto N = static_cast<int>(this->n_rows);
        auto KL = static_cast<int>(l_band);
        auto KU = static_cast<int>(u_band);
        auto LDAB = static_cast<int>(m_rows);
        const auto KLU = std::max(l_band, u_band);
        auto INFO = 0;

        if(std::is_same_v<T, float>) {
            using E = float;
            WORK.zeros(KLU * KLU * SPIKE(9));
            sspike_gbtrf_(SPIKE.memptr(), &N, &KL, &KU, (E*)this->memptr(), &LDAB, (E*)WORK.memptr(), &INFO);
        }
        else if(Precision::FULL == this->setting.precision) {
            using E = double;
            WORK.zeros(KLU * KLU * SPIKE(9));
            dspike_gbtrf_(SPIKE.memptr(), &N, &KL, &KU, (E*)this->memptr(), &LDAB, (E*)WORK.memptr(), &INFO);
        }
        else {
            this->s_memory = this->to_float();
            SWORK.zeros(KLU * KLU * SPIKE(9));
            sspike_gbtrf_(SPIKE.memptr(), &N, &KL, &KU, this->s_memory.mem, &LDAB, SWORK.memptr(), &INFO);
        }

        if(INFO != 0) {
            suanpan_error("solve() receives error code %u from the base driver, the matrix is probably singular.\n", INFO);
            return INFO;
        }

        this->factored = true;
    }

    return this->solve_trs(X, std::forward<Mat<T>>(B));
}

template<sp_d T> int BandMatSpike<T>::solve_trs(Mat<T>& X, Mat<T>&& B) {
    auto N = static_cast<int>(this->n_rows);
    auto KL = static_cast<int>(l_band);
    auto KU = static_cast<int>(u_band);
    auto NRHS = static_cast<int>(B.n_cols);
    auto LDAB = static_cast<int>(m_rows);
    auto LDB = static_cast<int>(B.n_rows);

    if(std::is_same_v<T, float>) {
        using E = float;
        sspike_gbtrs_(SPIKE.memptr(), &TRAN, &N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, (E*)WORK.memptr(), (E*)B.memptr(), &LDB);
        X = std::move(B);
    }
    else if(Precision::FULL == this->setting.precision) {
        using E = double;
        dspike_gbtrs_(SPIKE.memptr(), &TRAN, &N, &KL, &KU, &NRHS, (E*)this->memptr(), &LDAB, (E*)WORK.memptr(), (E*)B.memptr(), &LDB);
        X = std::move(B);
    }
    else {
        X = arma::zeros(B.n_rows, B.n_cols);

        auto multiplier = norm(B);

        auto counter = 0u;
        while(counter++ < this->setting.iterative_refinement) {
            if(multiplier < this->setting.tolerance) break;

            auto residual = conv_to<fmat>::from(B / multiplier);

            sspike_gbtrs_(SPIKE.memptr(), &TRAN, &N, &KL, &KU, &NRHS, this->s_memory.memptr(), &LDAB, SWORK.memptr(), residual.memptr(), &LDB);

            const mat incre = multiplier * conv_to<mat>::from(residual);

            X += incre;

            suanpan_debug("mixed precision algorithm multiplier: %.5E.\n", multiplier = arma::norm(B -= this->operator*(incre)));
        }
    }

    return SUANPAN_SUCCESS;
}

template<sp_d T> int BandMatSpike<T>::sign_det() const { throw invalid_argument("not supported"); }

#endif

//! @}
