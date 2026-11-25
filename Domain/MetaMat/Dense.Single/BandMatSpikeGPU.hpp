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
 * The `BandMatSpikeGPU` class uses GSPIKE, a GPU based SPIKE solver, and
 * supports single, double and mixed precision.
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

#include "../DenseMat.hpp"
#include "../cuda_ptr.hpp"

#include <gspike/spike.h>

template<sp_d T> class BandMatSpikeGPU final : public DenseMat<T> {
    static constexpr auto TRAN = 'N';

    static T bin;

    const uword l_band;
    const uword u_band;
    const uword m_rows; // memory block layout

    cuda_ptr d_A{}, d_B{};

    std::unique_ptr<gspike::solver<T>> solver;
    std::unique_ptr<gspike::mixed_solver> mixer;

    void init_config() {
        this->factored = false;
        solver.reset();
        mixer.reset();
    }

protected:
    using DenseMat<T>::direct_solve;

    int direct_solve(Mat<T>&, Mat<T>&&) override;

public:
    BandMatSpikeGPU(const uword in_size, const uword in_l, const uword in_u)
        : DenseMat<T>(in_size, in_size, (in_l + in_u + 1) * in_size)
        , l_band(in_l)
        , u_band(in_u)
        , m_rows(in_l + in_u + 1) { init_config(); }

    BandMatSpikeGPU(const BandMatSpikeGPU& other)
        : DenseMat<T>(other)
        , l_band(other.l_band)
        , u_band(other.u_band)
        , m_rows(other.m_rows) { init_config(); }

    BandMatSpikeGPU(BandMatSpikeGPU&&) = delete;
    BandMatSpikeGPU& operator=(const BandMatSpikeGPU&) = delete;
    BandMatSpikeGPU& operator=(BandMatSpikeGPU&&) = delete;

    unique_ptr<MetaMat<T>> unique_copy() override { return std::make_unique<BandMatSpikeGPU>(*this); }

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

template<sp_d T> T BandMatSpikeGPU<T>::bin = T(0);

template<sp_d T> Mat<T> BandMatSpikeGPU<T>::operator*(const Mat<T>& X) const {
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

template<sp_d T> int BandMatSpikeGPU<T>::direct_solve(Mat<T>& X, Mat<T>&& B) {
    suanpan_assert([&] { if(this->n_rows != this->n_cols) throw std::invalid_argument("requires a square matrix"); });

    const auto N = static_cast<int>(this->n_rows);
    const auto KL = static_cast<int>(l_band);
    const auto KU = static_cast<int>(u_band);
    const auto NRHS = static_cast<int>(B.n_cols);
    const auto LDAB = static_cast<int>(m_rows);
    const auto LDB = static_cast<int>(B.n_rows);

    if(!this->factored) {
        solver.reset();
        mixer.reset();

        d_A = cuda_ptr(sizeof(T), static_cast<int>(this->n_elem));
        cudaMemcpy(d_A.get<T>(), this->memptr(), sizeof(T) * this->n_elem, cudaMemcpyHostToDevice);

        if constexpr(std::is_same_v<T, float>) solver = std::make_unique<gspike::solver<float>>(N, KL, KU, d_A.get<float>(), LDAB);
        else if(Precision::FULL == this->setting.precision) solver = std::make_unique<gspike::solver<double>>(N, KL, KU, d_A.get<double>(), LDAB);
        else mixer = std::make_unique<gspike::mixed_solver>(N, KL, KU, d_A.get<double>(), LDAB);

        this->factored = true;
    }

    if(d_B.size < static_cast<int>(B.n_elem)) d_B = cuda_ptr(sizeof(T), static_cast<int>(B.n_elem));

    cudaMemcpy(d_B.get<T>(), B.memptr(), sizeof(T) * B.n_elem, cudaMemcpyHostToDevice);

    auto INFO = 0;
    if constexpr(std::is_same_v<T, float>) INFO = solver->solve(NRHS, d_B.get<float>(), LDB);
    else if(Precision::FULL == this->setting.precision) INFO = solver->solve(NRHS, d_B.get<double>(), LDB);
    else INFO = mixer->solve(NRHS, d_B.get<double>(), LDB);

    if(0 != INFO) {
        suanpan_error("Error code {} received, the matrix is probably singular.\n", INFO);
        return INFO;
    }

    cudaMemcpy(B.memptr(), d_B.get<T>(), sizeof(T) * B.n_elem, cudaMemcpyDeviceToHost);

    X = std::move(B);

    return SUANPAN_SUCCESS;
}

#endif

//! @}
