/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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
 * @class BandMatMAGMA
 * @brief A BandMatMAGMA class that holds matrices.
 *
 * @author tlc
 * @date 28/04/2024
 * @version 0.1.0
 * @file BandMatMAGMA.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef BANDMATMAGMA_HPP
#define BANDMATMAGMA_HPP

#ifdef SUANPAN_MAGMA

#include "BandMat.hpp"
#include "magma_v2.h"

template<sp_d T> class BandMatMAGMA final : public BandMat<T> {
    magma_queue_t queue{};

    void *APTR{}, *BPTR{};
    magma_int_t* IPIV{};

    void release() const {
        magma_free(APTR);
        magma_free(BPTR);
        magma_free(IPIV);
    }

protected:
    using BandMat<T>::direct_solve;

    int direct_solve(Mat<T>&, Mat<T>&&) override;

public:
    BandMatMAGMA(const uword in_size, const uword in_l, const uword in_u)
        : BandMat<T>(in_size, in_l, in_u) { magma_queue_create(0, &queue); }

    BandMatMAGMA(const BandMatMAGMA& other)
        : BandMat<T>(other) { magma_queue_create(0, &queue); }

    BandMatMAGMA(BandMatMAGMA&&) noexcept = delete;
    BandMatMAGMA& operator=(const BandMatMAGMA&) = delete;
    BandMatMAGMA& operator=(BandMatMAGMA&&) noexcept = delete;

    ~BandMatMAGMA() override {
        release();
        magma_queue_destroy(queue);
    }

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<BandMatMAGMA>(*this); }
};

template<sp_d T> int BandMatMAGMA<T>::direct_solve(Mat<T>& X, Mat<T>&& B) {
    suanpan_assert([&] { if(this->n_rows != this->n_cols) throw invalid_argument("requires a square matrix"); });

    auto INFO = 0;

    const auto N = static_cast<magma_int_t>(this->n_rows);
    const auto KL = static_cast<magma_int_t>(this->l_band);
    const auto KU = static_cast<magma_int_t>(this->u_band);
    const auto NRHS = static_cast<magma_int_t>(B.n_cols);
    const auto LDAB = static_cast<magma_int_t>(this->m_rows);
    const auto LDB = static_cast<magma_int_t>(B.n_rows);

    const auto LDDAB = magma_roundup(LDAB, 32);
    const auto LDDB = magma_roundup(LDB, 32);
    const auto SIZEA = LDDAB * N;
    const auto SIZEB = LDDB * NRHS;

    release();

    magma_imalloc(&IPIV, N);

    if constexpr(std::is_same_v<T, float>) {
        using E = float;

        magma_smalloc((magmaFloat_ptr*)&APTR, SIZEA);
        magma_smalloc((magmaFloat_ptr*)&BPTR, SIZEB);

        magma_ssetmatrix(LDAB, N, (E*)this->memptr(), LDAB, (magmaFloat_ptr)APTR, LDDAB, queue);
        magma_ssetmatrix(N, NRHS, (E*)B.memptr(), LDB, (magmaFloat_ptr)BPTR, LDDB, queue);
        magma_sgbsv_native(N, KL, KU, NRHS, (magmaFloat_ptr)APTR, LDDAB, IPIV, (magmaFloat_ptr)BPTR, LDDB, &INFO);
        magma_sgetmatrix(N, NRHS, (magmaFloat_ptr)BPTR, LDDB, (E*)B.memptr(), LDB, queue);

        X = std::move(B);
    }
    else if(Precision::FULL == this->setting.precision) {
        using E = double;

        magma_dmalloc((magmaDouble_ptr*)&APTR, SIZEA);
        magma_dmalloc((magmaDouble_ptr*)&BPTR, SIZEB);

        magma_dsetmatrix(LDAB, N, (E*)this->memptr(), LDAB, (magmaDouble_ptr)APTR, LDDAB, queue);
        magma_dsetmatrix(N, NRHS, (E*)B.memptr(), LDB, (magmaDouble_ptr)BPTR, LDDB, queue);
        magma_dgbsv_native(N, KL, KU, NRHS, (magmaDouble_ptr)APTR, LDDAB, IPIV, (magmaDouble_ptr)BPTR, LDDB, &INFO);
        magma_dgetmatrix(N, NRHS, (magmaDouble_ptr)BPTR, LDDB, (E*)B.memptr(), LDB, queue);

        X = std::move(B);
    }
    else {
        magma_smalloc((magmaFloat_ptr*)&APTR, SIZEA);
        magma_smalloc((magmaFloat_ptr*)&BPTR, SIZEB);

        this->s_memory = this->to_float();

        this->mixed_trs(X, std::forward<Mat<T>>(B), [&](fmat& residual) {
            magma_ssetmatrix(LDAB, N, this->s_memory.memptr(), LDAB, (magmaFloat_ptr)APTR, LDDAB, queue);
            magma_ssetmatrix(N, NRHS, residual.memptr(), LDB, (magmaFloat_ptr)BPTR, LDDB, queue);
            magma_sgbsv_native(N, KL, KU, NRHS, (magmaFloat_ptr)APTR, LDDAB, IPIV, (magmaFloat_ptr)BPTR, LDDB, &INFO);
            magma_sgetmatrix(N, NRHS, (magmaFloat_ptr)BPTR, LDDB, residual.memptr(), LDB, queue);

            return INFO;
        });
    }

    if(0 != INFO)
        suanpan_error("Error code {} received, the matrix is probably singular.\n", INFO);

    this->pivot.zeros(N);
    magma_igetmatrix(N, 1, IPIV, N, this->pivot.memptr(), N, queue);

    release();

    return INFO;
}

#endif

#endif

//! @}
