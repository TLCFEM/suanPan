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

#ifdef SUANPAN_CUDA

#include <cuda_runtime.h>
#include <cusolverDn.h>
#include "FullMat.hpp"

template<sp_d T> class FullMatCUDA final : public FullMat<T> {
    cusolverDnHandle_t handle = nullptr;
    cudaStream_t stream = nullptr;

    int* info = nullptr;
    int* ipiv = nullptr;
    void* d_A = nullptr;
    void* buffer = nullptr;

    void acquire();
    void release() const;

public:
    FullMatCUDA(uword, uword);
    FullMatCUDA(const FullMatCUDA&);
    FullMatCUDA(FullMatCUDA&&) noexcept = delete;
    FullMatCUDA& operator=(const FullMatCUDA&) = delete;
    FullMatCUDA& operator=(FullMatCUDA&&) noexcept = delete;
    ~FullMatCUDA() override;

    unique_ptr<MetaMat<T>> make_copy() override;

    int solve(Mat<T>&, Mat<T>&&) override;
    int solve(Mat<T>&, const Mat<T>&) override;
};

template<sp_d T> void FullMatCUDA<T>::acquire() {
    cusolverDnCreate(&handle);
    cudaStreamCreate(&stream);
    cusolverDnSetStream(handle, stream);

    cudaMalloc(&info, sizeof(int));
    cudaMemset(info, 0, sizeof(int));
    cudaMalloc(&ipiv, sizeof(int) * this->n_rows);

    if(int bufferSize = 0; std::is_same_v<T, float> || Precision::MIXED == this->precision) {
        cudaMalloc(&d_A, sizeof(float) * this->n_elem);
        cusolverDnSgetrf_bufferSize(handle, this->n_rows, this->n_cols, (float*)d_A, this->n_elem, &bufferSize);
        cudaMalloc(&buffer, sizeof(float) * bufferSize);
    }
    else {
        cudaMalloc(&d_A, sizeof(double) * this->n_elem);
        cusolverDnDgetrf_bufferSize(handle, this->n_rows, this->n_cols, (double*)d_A, this->n_elem, &bufferSize);
        cudaMalloc(&buffer, sizeof(double) * bufferSize);
    }
}

template<sp_d T> void FullMatCUDA<T>::release() const {
    if(handle) cusolverDnDestroy(handle);
    if(stream) cudaStreamDestroy(stream);

    if(info) cudaFree(info);
    if(d_A) cudaFree(d_A);
    if(buffer) cudaFree(buffer);
    if(ipiv) cudaFree(ipiv);
}

template<sp_d T> FullMatCUDA<T>::FullMatCUDA(const uword in_rows, const uword in_cols)
    : FullMat<T>(in_rows, in_cols) { acquire(); }

template<sp_d T> FullMatCUDA<T>::FullMatCUDA(const FullMatCUDA& other)
    : FullMat<T>(other) { acquire(); }

template<sp_d T> FullMatCUDA<T>::~FullMatCUDA() { release(); }

template<sp_d T> unique_ptr<MetaMat<T>> FullMatCUDA<T>::make_copy() { return make_unique<FullMatCUDA<T>>(*this); }

template<sp_d T> int FullMatCUDA<T>::solve(Mat<T>& X, Mat<T>&& B) { return solve(X, B); }

template<sp_d T> int FullMatCUDA<T>::solve(Mat<T>& X, const Mat<T>& B) {
    if(std::is_same_v<T, float>) {
        // pure float
        if(!this->factored) {
            cudaMemcpyAsync(d_A, this->memptr(), sizeof(float) * this->n_elem, cudaMemcpyHostToDevice, stream);
            cusolverDnSgetrf(handle, this->n_rows, this->n_cols, (float*)d_A, this->n_rows, (float*)buffer, ipiv, info);

            this->factored = true;
        }

        const size_t byte_size = sizeof(float) * B.n_elem;

        void* d_x = nullptr;
        cudaMalloc(&d_x, byte_size);
        cudaMemcpyAsync(d_x, B.memptr(), byte_size, cudaMemcpyHostToDevice, stream);
        cusolverDnSgetrs(handle, CUBLAS_OP_N, this->n_rows, B.n_cols, (float*)d_A, this->n_rows, ipiv, (float*)d_x, this->n_rows, info);

        X.set_size(arma::size(B));

        cudaMemcpyAsync(X.memptr(), d_x, byte_size, cudaMemcpyDeviceToHost, stream);

        cudaDeviceSynchronize();

        if(d_x) cudaFree(d_x);
    }
    else if(Precision::MIXED == this->precision) {
        // mixed precision
        if(!this->factored) {
            this->s_memory = this->to_float();

            cudaMemcpyAsync(d_A, this->s_memory.memptr(), sizeof(float) * this->s_memory.n_elem, cudaMemcpyHostToDevice, stream);
            cusolverDnSgetrf(handle, this->n_rows, this->n_cols, (float*)d_A, this->n_rows, (float*)buffer, ipiv, info);

            this->factored = true;
        }

        const size_t byte_size = sizeof(float) * B.n_elem;

        void* d_x = nullptr;
        cudaMalloc(&d_x, byte_size);

        X = arma::zeros(B.n_rows, B.n_cols);

        mat full_residual = B;

        auto multiplier = arma::norm(full_residual);

        auto counter = 0u;
        while(counter++ < this->refinement) {
            if(multiplier < this->tolerance) break;

            auto residual = conv_to<fmat>::from(full_residual / multiplier);

            cudaMemcpyAsync(d_x, residual.memptr(), byte_size, cudaMemcpyHostToDevice, stream);
            cusolverDnSgetrs(handle, CUBLAS_OP_N, this->n_rows, B.n_cols, (float*)d_A, this->n_rows, ipiv, (float*)d_x, this->n_rows, info);
            cudaMemcpyAsync(residual.memptr(), d_x, byte_size, cudaMemcpyDeviceToHost, stream);

            cudaDeviceSynchronize();

            const mat incre = multiplier * conv_to<mat>::from(residual);

            X += incre;

            suanpan_debug("mixed precision algorithm multiplier: %.5E\n", multiplier = arma::norm(full_residual -= this->operator*(incre)));
        }

        if(d_x) cudaFree(d_x);
    }
    else {
        // pure double
        if(!this->factored) {
            cudaMemcpyAsync(d_A, this->memptr(), sizeof(double) * this->n_elem, cudaMemcpyHostToDevice, stream);
            cusolverDnDgetrf(handle, this->n_rows, this->n_cols, (double*)d_A, this->n_rows, (double*)buffer, ipiv, info);

            this->factored = true;
        }

        const size_t byte_size = sizeof(double) * B.n_elem;

        void* d_x = nullptr;
        cudaMalloc(&d_x, byte_size);
        cudaMemcpyAsync(d_x, B.memptr(), byte_size, cudaMemcpyHostToDevice, stream);
        cusolverDnDgetrs(handle, CUBLAS_OP_N, this->n_rows, B.n_cols, (double*)d_A, this->n_rows, ipiv, (double*)d_x, this->n_rows, info);

        X.set_size(arma::size(B));

        cudaMemcpyAsync(X.memptr(), d_x, byte_size, cudaMemcpyDeviceToHost, stream);

        cudaDeviceSynchronize();

        if(d_x) cudaFree(d_x);
    }

    return 0;
}

#endif

#endif

//! @}
