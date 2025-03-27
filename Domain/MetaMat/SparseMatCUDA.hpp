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
 * @class SparseMatCUDA
 * @brief A SparseMatCUDA class that holds matrices.
 *
 * The `SparseMatCUDA` class uses CUDA and supports single, double, mixed precision.
 *
 * @author tlc
 * @date 21/04/2021
 * @version 0.1.0
 * @file SparseMatCUDA.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef SPARSEMATCUDA_HPP
#define SPARSEMATCUDA_HPP

#ifdef SUANPAN_CUDA

#include <cusolverSp.h>
#include <cusparse.h>
#include "SparseMat.hpp"
#include "csr_form.hpp"

template<sp_d T> class SparseMatCUDA final : public SparseMat<T> {
    cusolverSpHandle_t handle = nullptr;
    cudaStream_t stream = nullptr;
    cusparseMatDescr_t descr = nullptr;

    void* d_val_idx = nullptr;
    void* d_col_idx = nullptr;
    void* d_row_ptr = nullptr;

    void acquire() {
        cusolverSpCreate(&handle);
        cudaStreamCreate(&stream);
        cusolverSpSetStream(handle, stream);
        cusparseCreateMatDescr(&descr);
        cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
        cudaMalloc(&d_row_ptr, sizeof(int) * (this->n_rows + 1));
    }

    void release() const {
        if(handle) cusolverSpDestroy(handle);
        if(stream) cudaStreamDestroy(stream);
        if(descr) cusparseDestroyMatDescr(descr);
        if(d_row_ptr) cudaFree(d_row_ptr);
    }

    template<sp_d ET> void device_alloc(csr_form<ET, int>&& csr_mat) {
        const size_t n_val = sizeof(ET) * csr_mat.n_elem;
        const size_t n_col = sizeof(int) * csr_mat.n_elem;

        cudaMalloc(&d_val_idx, n_val);
        cudaMalloc(&d_col_idx, n_col);

        cudaMemcpyAsync(d_val_idx, csr_mat.val_mem(), n_val, cudaMemcpyHostToDevice, stream);
        cudaMemcpyAsync(d_col_idx, csr_mat.col_mem(), n_col, cudaMemcpyHostToDevice, stream);
        cudaMemcpyAsync(d_row_ptr, csr_mat.row_mem(), sizeof(int) * (csr_mat.n_rows + 1llu), cudaMemcpyHostToDevice, stream);
    }

    void device_dealloc() const {
        if(d_val_idx) cudaFree(d_val_idx);
        if(d_col_idx) cudaFree(d_col_idx);
    }

protected:
    using SparseMat<T>::direct_solve;

    int direct_solve(Mat<T>&, const Mat<T>&) override;

public:
    SparseMatCUDA(const uword in_row, const uword in_col, const uword in_elem = 0)
        : SparseMat<T>(in_row, in_col, in_elem) { acquire(); }

    SparseMatCUDA(const SparseMatCUDA& other)
        : SparseMat<T>(other) { acquire(); }

    SparseMatCUDA(SparseMatCUDA&&) noexcept = delete;
    SparseMatCUDA& operator=(const SparseMatCUDA&) = delete;
    SparseMatCUDA& operator=(SparseMatCUDA&&) noexcept = delete;

    ~SparseMatCUDA() override {
        release();
        device_dealloc();
    }

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SparseMatCUDA>(*this); }
};

template<sp_d T> int SparseMatCUDA<T>::direct_solve(Mat<T>& X, const Mat<T>& B) {
    if(!this->factored) {
        // deallocate memory previously allocated for csr matrix
        device_dealloc();

        std::is_same_v<T, float> || Precision::MIXED == this->setting.precision ? device_alloc(csr_form<float, int>(this->triplet_mat)) : device_alloc(csr_form<double, int>(this->triplet_mat));

        this->factored = true;
    }

    const size_t n_rhs = (std::is_same_v<T, float> || Precision::MIXED == this->setting.precision ? sizeof(float) : sizeof(double)) * B.n_elem;

    void* d_b = nullptr;
    void* d_x = nullptr;

    cudaMalloc(&d_b, n_rhs);
    cudaMalloc(&d_x, n_rhs);

    int singularity;
    auto code = 0;

    if constexpr(std::is_same_v<T, float>) {
        cudaMemcpyAsync(d_b, B.memptr(), n_rhs, cudaMemcpyHostToDevice, stream);

        for(auto I = 0llu; I < B.n_elem; I += B.n_rows) code += cusolverSpScsrlsvqr(handle, int(this->n_rows), int(this->triplet_mat.n_elem), descr, (float*)d_val_idx, (int*)d_row_ptr, (int*)d_col_idx, (float*)d_b + I, float(this->setting.tolerance), 3, (float*)d_x + I, &singularity);

        X.set_size(arma::size(B));

        cudaMemcpyAsync(X.memptr(), d_x, n_rhs, cudaMemcpyDeviceToHost, stream);

        cudaDeviceSynchronize();
    }
    else if(Precision::FULL == this->setting.precision) {
        cudaMemcpyAsync(d_b, B.memptr(), n_rhs, cudaMemcpyHostToDevice, stream);

        for(auto I = 0llu; I < B.n_elem; I += B.n_rows) code += cusolverSpDcsrlsvqr(handle, int(this->n_rows), int(this->triplet_mat.n_elem), descr, (double*)d_val_idx, (int*)d_row_ptr, (int*)d_col_idx, (double*)d_b + I, this->setting.tolerance, 3, (double*)d_x + I, &singularity);

        X.set_size(arma::size(B));

        cudaMemcpyAsync(X.memptr(), d_x, n_rhs, cudaMemcpyDeviceToHost, stream);

        cudaDeviceSynchronize();
    }
    else {
        X = arma::zeros(arma::size(B));

        mat full_residual = B;

        auto multiplier = norm(full_residual);

        auto counter = std::uint8_t{0};
        while(counter++ < this->setting.iterative_refinement) {
            if(multiplier < this->setting.tolerance) break;

            auto residual = conv_to<fmat>::from(full_residual / multiplier);

            cudaMemcpyAsync(d_b, residual.memptr(), n_rhs, cudaMemcpyHostToDevice, stream);

            code = 0;
            for(auto I = 0llu; I < B.n_elem; I += B.n_rows) code += cusolverSpScsrlsvqr(handle, int(this->n_rows), int(this->triplet_mat.n_elem), descr, (float*)d_val_idx, (int*)d_row_ptr, (int*)d_col_idx, (float*)d_b + I, float(this->setting.tolerance), 3, (float*)d_x + I, &singularity);
            if(0 != code) break;

            cudaMemcpyAsync(residual.memptr(), d_x, n_rhs, cudaMemcpyDeviceToHost, stream);

            cudaDeviceSynchronize();

            const mat incre = multiplier * conv_to<mat>::from(residual);

            X += incre;

            suanpan_debug("Mixed precision algorithm multiplier: {:.5E}.\n", multiplier = arma::norm(full_residual -= this->operator*(incre)));
        }
    }

    if(d_b) cudaFree(d_b);
    if(d_x) cudaFree(d_x);

    return 0 == code ? SUANPAN_SUCCESS : SUANPAN_FAIL;
}

#endif

#endif

//! @}
