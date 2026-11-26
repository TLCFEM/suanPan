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
 * @class BandMatSpikeGPU
 * @brief A BandMatSpikeGPU class that holds matrices.
 *
 * The `BandMatSpikeGPU` class uses GSPIKE, a GPU based SPIKE solver. Single
 * precision matrices are solved directly, while double precision matrices
 * always use the mixed precision solver with iterative refinement.
 *
 * @author tlc
 * @date 25/11/2026
 * @version 0.1.0
 * @file BandMatSpikeGPU.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef BANDMATSPIKEGPU_HPP
#define BANDMATSPIKEGPU_HPP

#include "../cuda_ptr.hpp"
#include "BandMatSpikeBase.hpp"

#include <gspike/spike.h>

template<sp_d T> class BandMatSpikeGPU final : public BandMatSpikeBase<T> {
    cuda_ptr d_A{}, d_B{};

    std::unique_ptr<gspike::solver<float>> solver;
    std::unique_ptr<gspike::mixed_solver> mixer;

    void init_config() {
        this->factored = false;
        solver.reset();
        mixer.reset();
    }

protected:
    using BandMatSpikeBase<T>::direct_solve;

    int direct_solve(Mat<T>&, Mat<T>&&) override;

public:
    BandMatSpikeGPU(const uword in_size, const uword in_l, const uword in_u)
        : BandMatSpikeBase<T>(in_size, in_l, in_u) { init_config(); }

    BandMatSpikeGPU(const BandMatSpikeGPU& other)
        : BandMatSpikeBase<T>(other) { init_config(); }

    BandMatSpikeGPU(BandMatSpikeGPU&&) = delete;
    BandMatSpikeGPU& operator=(const BandMatSpikeGPU&) = delete;
    BandMatSpikeGPU& operator=(BandMatSpikeGPU&&) = delete;

    unique_ptr<MetaMat<T>> unique_copy() override { return std::make_unique<BandMatSpikeGPU>(*this); }

    [[nodiscard]] int sign_det() const override {
        if(!this->factored) return 0;

        if constexpr(std::is_same_v<T, float>) return nullptr == solver ? 0 : solver->det_sign();
        else return nullptr == mixer ? 0 : mixer->det_sign();
    }
};

template<sp_d T> int BandMatSpikeGPU<T>::direct_solve(Mat<T>& X, Mat<T>&& B) {
    suanpan_assert([&] { if(this->n_rows != this->n_cols) throw std::invalid_argument("requires a square matrix"); });

    const auto N = static_cast<int>(this->n_rows);
    const auto KL = static_cast<int>(this->l_band);
    const auto KU = static_cast<int>(this->u_band);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDAB = static_cast<int>(this->m_rows);
    const auto LDB = static_cast<int>(B.n_rows);

    if(!this->factored) {
        solver.reset();
        mixer.reset();

        d_A = cuda_ptr(sizeof(T), static_cast<int>(this->n_elem));
        d_A.copy_from(this->memptr());

        if constexpr(std::is_same_v<T, float>) solver = std::make_unique<gspike::solver<float>>(N, KL, KU, d_A.get<float>(), LDAB);
        else mixer = std::make_unique<gspike::mixed_solver>(N, KL, KU, d_A.get<double>(), LDAB);

        this->factored = true;
    }

    // sync copies transfer the whole buffer, so it must match the RHS exactly
    if(d_B.size != static_cast<int>(B.n_elem)) d_B = cuda_ptr(sizeof(T), static_cast<int>(B.n_elem));

    d_B.copy_from(B.memptr());

    auto INFO = 0;
    if constexpr(std::is_same_v<T, float>) INFO = solver->solve(NRHS, d_B.get<float>(), LDB);
    else INFO = mixer->solve(NRHS, d_B.get<double>(), LDB);

    if(0 != INFO) {
        suanpan_error("Error code {} received, the matrix is probably singular.\n", INFO);
        return INFO;
    }

    d_B.copy_to(B.memptr());

    X = std::move(B);

    return SUANPAN_SUCCESS;
}

#endif

//! @}
