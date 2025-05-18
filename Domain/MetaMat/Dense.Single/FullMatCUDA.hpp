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
 * @class FullMatCUDA
 * @brief A FullMatCUDA class that holds matrices.
 *
 * @author tlc
 * @date 17/04/2021
 * @version 0.1.0
 * @file FullMatCUDA.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef FULLMATCUDA_HPP
#define FULLMATCUDA_HPP

#include "../cuda_ptr.hpp"
#include "FullMat.hpp"

#include <cusolverDn.h>

template<sp_d T> class FullMatCUDA final : public FullMat<T> {
    cusolverDnHandle_t handle = nullptr;
    cudaStream_t stream = nullptr;

    cuda_ptr info{sizeof(int), 1}, d_ipiv{sizeof(int), static_cast<int>(this->n_rows)}, d_A{}, d_work{};

    void acquire() {
        cusolverDnCreate(&handle);
        cudaStreamCreate(&stream);
        cusolverDnSetStream(handle, stream);

        int work_size = 0;
        if(std::is_same_v<T, float> || Precision::MIXED == this->setting.precision) {
            d_A = cuda_ptr(sizeof(float), static_cast<int>(this->n_elem));
            cusolverDnSgetrf_bufferSize(handle, static_cast<int>(this->n_rows), static_cast<int>(this->n_cols), d_A.get<float>(), d_A.size, &work_size);
            d_work = cuda_ptr(sizeof(float), work_size);
        }
        else {
            d_A = cuda_ptr(sizeof(double), static_cast<int>(this->n_elem));
            cusolverDnDgetrf_bufferSize(handle, static_cast<int>(this->n_rows), static_cast<int>(this->n_cols), d_A.get<double>(), d_A.size, &work_size);
            d_work = cuda_ptr(sizeof(double), work_size);
        }

        this->factored = false;
    }

    void release() const {
        if(handle) cusolverDnDestroy(handle);
        if(stream) cudaStreamDestroy(stream);
    }

protected:
    int direct_solve(Mat<T>& X, Mat<T>&& B) override { return this->direct_solve(X, B); }

    int direct_solve(Mat<T>&, const Mat<T>&) override;

public:
    FullMatCUDA(const uword in_rows, const uword in_cols)
        : FullMat<T>(in_rows, in_cols) { acquire(); }

    FullMatCUDA(const FullMatCUDA& other)
        : FullMat<T>(other) { acquire(); }

    FullMatCUDA(FullMatCUDA&&) = delete;
    FullMatCUDA& operator=(const FullMatCUDA&) = delete;
    FullMatCUDA& operator=(FullMatCUDA&&) = delete;

    ~FullMatCUDA() override { release(); }

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<FullMatCUDA>(*this); }
};

template<sp_d T> int FullMatCUDA<T>::direct_solve(Mat<T>& X, const Mat<T>& B) {
    const auto NROW = static_cast<int>(this->n_rows), NCOL = static_cast<int>(this->n_cols);

    int flag;

    if constexpr(std::is_same_v<T, float>) {
        // pure float
        if(!this->factored) {
            this->factored = true;
            d_A.copy_from(this->memptr(), stream);
            cusolverDnSgetrf(handle, NROW, NCOL, d_A.get<float>(), NROW, d_work.get<float>(), d_ipiv.get(), info.get());
        }

        const cuda_ptr d_x{sizeof(float), static_cast<int>(B.n_elem)};
        d_x.copy_from(B.memptr(), stream);

        cusolverDnSgetrs(handle, CUBLAS_OP_N, NROW, static_cast<int>(B.n_cols), d_A.get<float>(), NROW, d_ipiv.get(), d_x.get<float>(), NROW, info.get());

        X.set_size(arma::size(B));
        d_x.copy_to(X.memptr(), stream);
    }
    else if(Precision::MIXED == this->setting.precision) {
        // mixed precision
        if(!this->factored) {
            this->factored = true;
            this->s_memory = this->to_float();
            d_A.copy_from(this->s_memory.memptr(), stream);
            cusolverDnSgetrf(handle, NROW, NCOL, d_A.get<float>(), NROW, d_work.get<float>(), d_ipiv.get(), info.get());
        }

        const cuda_ptr d_x{sizeof(float), static_cast<int>(B.n_elem)};

        X = arma::zeros(B.n_rows, B.n_cols);

        mat full_residual = B;

        std::uint8_t counter{0};
        while(counter++ < this->setting.iterative_refinement) {
            const auto multiplier = norm(full_residual);
            if(multiplier < this->setting.tolerance) break;
            suanpan_debug("Mixed precision algorithm multiplier: {:.5E}.\n", multiplier);

            auto residual = conv_to<fmat>::from(full_residual / multiplier);
            d_x.copy_from(residual.memptr(), stream);

            cusolverDnSgetrs(handle, CUBLAS_OP_N, NROW, static_cast<int>(B.n_cols), d_A.get<float>(), NROW, d_ipiv.get(), d_x.get<float>(), NROW, info.get());

            d_x.copy_to(residual.memptr(), stream);
            full_residual = B - this->operator*(X += multiplier * conv_to<mat>::from(residual));
        }
    }
    else {
        // pure double
        if(!this->factored) {
            this->factored = true;
            d_A.copy_from(this->memptr(), stream);
            cusolverDnDgetrf(handle, NROW, NCOL, d_A.get<double>(), NROW, d_work.get<double>(), d_ipiv.get(), info.get());
        }

        const cuda_ptr d_x{sizeof(float), static_cast<int>(B.n_elem)};
        d_x.copy_from(B.memptr(), stream);

        cusolverDnDgetrs(handle, CUBLAS_OP_N, NROW, static_cast<int>(B.n_cols), d_A.get<double>(), NROW, d_ipiv.get(), d_x.get<double>(), NROW, info.get());

        X.set_size(arma::size(B));
        d_x.copy_to(X.memptr(), stream);
    }

    info.copy_to(&flag, stream);

    return flag;
}

#endif

//! @}
