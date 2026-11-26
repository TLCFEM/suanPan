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

#include "BandMatSpikeBase.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

extern "C" {
void spikeinit_(la_it*, la_it*, la_it*);
void dspike_tune_(la_it*);
void dspike_gbtrf_(la_it*, la_it*, la_it*, la_it*, double*, la_it*, double*, la_it*);
void dspike_gbtrs_(la_it*, const char*, la_it*, la_it*, la_it*, la_it*, double*, la_it*, double*, double*, la_it*);
void sspike_tune_(la_it*);
void sspike_gbtrf_(la_it*, la_it*, la_it*, la_it*, float*, la_it*, float*, la_it*);
void sspike_gbtrs_(la_it*, const char*, la_it*, la_it*, la_it*, la_it*, float*, la_it*, float*, float*, la_it*);
void spikerl_calc_size_partitions_(la_it*, la_it*, la_it*);
}

namespace spike_detail {
    // Accumulate only the sign so overflow or underflow cannot invalidate it.
    template<sp_d T> struct det_sign {
        int sign = 1;
        bool zero = false;

        void mul(const T V) {
            if(V == T(0)) zero = true;
            else if(V < T(0)) sign = -sign;
        }

        void mul(const det_sign& other) {
            sign *= other.sign;
            zero = zero || other.zero;
        }
    };

    // Convert the one-based Fortran partition offsets to zero-based starts.
    inline std::vector<la_it> partition_starts(const la_it* SPIKE, const la_it N) {
        la_it spm[64]{};
        std::copy_n(SPIKE, 64, spm);

        la_it n = N;
        std::vector<la_it> ajmin(static_cast<size_t>(spm[19]) + 1, 1);
        spikerl_calc_size_partitions_(spm, ajmin.data(), &n);
        for(auto& I : ajmin) --I;

        return ajmin;
    }

    // Dense LU of the interface system.
    template<sp_d T> det_sign<T> interface_sign(const T* M0, const la_it M) {
        det_sign<T> result;

        if(0 == M) return result;

        std::vector<T> matrix(M0, M0 + static_cast<size_t>(M) * M);
        const auto at = [&](const la_it I, const la_it J) -> T& { return matrix[static_cast<size_t>(J) * M + I]; };

        for(la_it J = 0; J < M; ++J) {
            la_it pivot = J;
            for(la_it I = J + 1; I < M; ++I)
                if(std::fabs(at(I, J)) > std::fabs(at(pivot, J))) pivot = I;
            if(pivot != J) {
                for(la_it K = 0; K < M; ++K) std::swap(at(J, K), at(pivot, K));
                result.sign = -result.sign;
            }
            const T D = at(J, J);
            if(D == T(0)) {
                result.zero = true;
                return result;
            }
            result.mul(D);
            for(la_it I = J + 1; I < M; ++I) {
                at(I, J) /= D;
                for(la_it K = J + 1; K < M; ++K) at(I, K) -= at(I, J) * at(J, K);
            }
        }

        return result;
    }

    // det(A) = prod_p det(A_pp) * det(R), with R the reduced interface system.
    // Two-thread partitions keep part of their determinant outside WORK.
    template<sp_d T> int determinant_sign(const la_it* SPIKE, const la_it N, const la_it KL, const la_it KU, const T* A, const la_it LDA, const T* WORK) {
        if(nullptr == SPIKE || nullptr == A || N < 1 || KL < 0 || KU < 0 || LDA < KL + KU + 1) return 0;

        const auto klu = std::max(KL, KU);
        const auto nbpart = SPIKE[19];

        if(nbpart < 1 || nbpart > 63) return 0;
        if(0 < SPIKE[22]) return 0;
        if(nbpart > 1 && (klu < 1 || SPIKE[20] < 1 || nullptr == WORK)) return 0;

        det_sign<T> total;

        const auto ajmin = partition_starts(SPIKE, N);

        // spike_gbtrf does not pivot, the factored diagonal is in band row ku.
        for(la_it p = 0; p < nbpart; ++p) {
            const auto start = ajmin[p];
            const auto size = ajmin[p + 1] - start;
            for(la_it j = 0; j < size; ++j) total.mul(A[static_cast<size_t>(start + j) * LDA + KU]);
        }

        // Rebuild the level zero interface from the saved spikes.
        // rV/rW: (2*klu, klu, nblevel, nbpart), Fortran column major.
        if(nbpart > 1 && klu > 0) {
            const auto nblevel = SPIKE[20];
            const size_t plane = static_cast<size_t>(2) * klu * klu * nblevel;
            const size_t rw0 = plane * nbpart;
            const la_it m = 2 * klu * (nbpart - 1);
            std::vector<T> R(static_cast<size_t>(m) * m, T(0));
            const auto r = [&](const la_it I, const la_it J) -> T& { return R[static_cast<size_t>(J) * m + I]; };
            const auto v = [&](const la_it I, const la_it J, const la_it P) -> T { return WORK[static_cast<size_t>(P) * plane + static_cast<size_t>(J) * 2 * klu + I]; };
            const auto w = [&](const la_it I, const la_it J, const la_it P) -> T { return WORK[rw0 + static_cast<size_t>(P) * plane + static_cast<size_t>(J) * 2 * klu + I]; };

            for(la_it i = 0; i < nbpart - 1; ++i) {
                const la_it I0 = 2 * klu * i;
                for(la_it a = 0; a < klu; ++a) {
                    r(I0 + a, I0 + a) = T(1);
                    r(I0 + klu + a, I0 + klu + a) = T(1);
                    for(la_it b = 0; b < klu; ++b) {
                        r(I0 + a, I0 + klu + b) = v(klu + a, b, i); // Vp_i^bot
                        r(I0 + klu + a, I0 + b) = w(a, b, i + 1);   // Vm_{i+1}^top
                        if(i > 0)
                            r(I0 + a, I0 - 2 * klu + b) = w(klu + a, b, i); // Vm_i^bot
                        if(i + 1 < nbpart - 1)
                            r(I0 + klu + a, I0 + 3 * klu + b) = v(a, b, i + 1); // Vp_{i+1}^top
                    }
                }
            }

            total.mul(interface_sign(R.data(), m));
        }

        return total.zero ? 0 : total.sign;
    }
} // namespace spike_detail

template<sp_d T> class BandMatSpike final : public BandMatSpikeBase<T> {
    static constexpr auto TRAN = 'N';

    static la_it SPROTO, DPROTO;

    la_it SPIKE[64]{};

    podarray<T> WORK;
    podarray<float> SWORK;

    void init_spike() {
        auto N = static_cast<la_it>(this->n_rows);
        auto KLU = static_cast<la_it>(std::max(this->l_band, this->u_band));

        spikeinit_(SPIKE, &N, &KLU);

        SPIKE[6] = std::is_same_v<T, float> ? SPROTO : DPROTO;
        SPIKE[4] = SPIKE[6] + SPIKE[6] / 2 + 10;
        SPIKE[3] = SPIKE[4] / 2;
    }

    int solve_trs(Mat<T>&, Mat<T>&&);

protected:
    using BandMatSpikeBase<T>::direct_solve;

    int direct_solve(Mat<T>&, Mat<T>&&) override;

public:
    BandMatSpike(const uword in_size, const uword in_l, const uword in_u)
        : BandMatSpikeBase<T>(in_size, in_l, in_u) { init_spike(); }

    BandMatSpike(const BandMatSpike& other)
        : BandMatSpikeBase<T>(other) { init_spike(); }

    BandMatSpike(BandMatSpike&&) = delete;
    BandMatSpike& operator=(const BandMatSpike&) = delete;
    BandMatSpike& operator=(BandMatSpike&&) = delete;

    unique_ptr<MetaMat<T>> unique_copy() override { return std::make_unique<BandMatSpike>(*this); }

    [[nodiscard]] int sign_det() const override {
        if(!this->factored) return 0;

        const auto N = static_cast<la_it>(this->n_rows);
        const auto KL = static_cast<la_it>(this->l_band);
        const auto KU = static_cast<la_it>(this->u_band);
        const auto LDA = static_cast<la_it>(this->m_rows);

        if constexpr(std::is_same_v<T, float>) return spike_detail::determinant_sign(SPIKE, N, KL, KU, this->memptr(), LDA, WORK.memptr());
        else if(Precision::FULL == this->setting.precision) return spike_detail::determinant_sign(SPIKE, N, KL, KU, this->memptr(), LDA, WORK.memptr());
        else return spike_detail::determinant_sign(SPIKE, N, KL, KU, this->s_memory.memptr(), LDA, SWORK.memptr());
    }
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

template<sp_d T> int BandMatSpike<T>::direct_solve(Mat<T>& X, Mat<T>&& B) {
    if(!this->factored) {
        suanpan_assert([&] { if(this->n_rows != this->n_cols) throw std::invalid_argument("requires a square matrix"); });

        la_it INFO = 0;

        auto N = static_cast<la_it>(this->n_rows);
        auto KL = static_cast<la_it>(this->l_band);
        auto KU = static_cast<la_it>(this->u_band);
        auto LDAB = static_cast<la_it>(this->m_rows);
        const auto KLU = std::max(this->l_band, this->u_band);
        this->factored = true;

        if constexpr(std::is_same_v<T, float>) {
            WORK.zeros(KLU * KLU * SPIKE[9]);
            sspike_gbtrf_(SPIKE, &N, &KL, &KU, this->memptr(), &LDAB, WORK.memptr(), &INFO);
        }
        else if(Precision::FULL == this->setting.precision) {
            WORK.zeros(KLU * KLU * SPIKE[9]);
            dspike_gbtrf_(SPIKE, &N, &KL, &KU, this->memptr(), &LDAB, WORK.memptr(), &INFO);
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
    auto KL = static_cast<la_it>(this->l_band);
    auto KU = static_cast<la_it>(this->u_band);
    auto NRHS = static_cast<la_it>(B.n_cols);
    auto LDAB = static_cast<la_it>(this->m_rows);
    auto LDB = static_cast<la_it>(B.n_rows);

    if constexpr(std::is_same_v<T, float>) {
        sspike_gbtrs_(SPIKE, &TRAN, &N, &KL, &KU, &NRHS, this->memptr(), &LDAB, WORK.memptr(), B.memptr(), &LDB);
        X = std::move(B);
    }
    else if(Precision::FULL == this->setting.precision) {
        dspike_gbtrs_(SPIKE, &TRAN, &N, &KL, &KU, &NRHS, this->memptr(), &LDAB, WORK.memptr(), B.memptr(), &LDB);
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
