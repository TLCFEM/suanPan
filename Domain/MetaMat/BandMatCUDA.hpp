/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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
 * @class BandMatCUDA
 * @brief A BandMatCUDA class that holds matrices.
 *
 * @author tlc
 * @date 16/02/2023
 * @version 0.1.0
 * @file BandMatCUDA.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef BANDMATCUDA_HPP
#define BANDMATCUDA_HPP

#ifdef SUANPAN_CUDA

#include "BandMat.hpp"
#include <cusolverSp.h>
#include <cusparse.h>
#include "csr_form.hpp"

template<sp_d T> class BandMatCUDA final : public BandMat<T> {
    cusolverSpHandle_t handle = nullptr;
    cudaStream_t stream = nullptr;
    cusparseMatDescr_t descr = nullptr;

    void* d_val_idx = nullptr;
    void* d_col_idx = nullptr;
    void* d_row_ptr = nullptr;

    triplet_form<float, int> s_mat{static_cast<int>(this->n_rows), static_cast<int>(this->n_cols), static_cast<int>(this->n_elem)};

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

    void device_alloc(csr_form<float, int>&& csr_mat) {
        const size_t n_val = sizeof(float) * csr_mat.n_elem;
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
    using BandMat<T>::direct_solve;

    int direct_solve(Mat<T>&, Mat<T>&&) override;

public:
    BandMatCUDA(const uword in_size, const uword in_l, const uword in_u)
        : BandMat<T>(in_size, in_l, in_u) { acquire(); }

    BandMatCUDA(const BandMatCUDA& other)
        : BandMat<T>(other) { acquire(); }

    BandMatCUDA(BandMatCUDA&&) noexcept = delete;
    BandMatCUDA& operator=(const BandMatCUDA&) = delete;
    BandMatCUDA& operator=(BandMatCUDA&&) noexcept = delete;

    ~BandMatCUDA() override {
        release();
        device_dealloc();
    }

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<BandMatCUDA>(*this); }
};

template<sp_d T> int BandMatCUDA<T>::direct_solve(Mat<T>& X, Mat<T>&& B) {
    suanpan_assert([&] { if(this->n_rows != this->n_cols) throw invalid_argument("requires a square matrix"); });

    if(!this->factored) {
        this->factored = true;

        device_dealloc();

        s_mat.zeros();
        for(auto I = 0; I < static_cast<int>(this->n_rows); ++I) for(auto J = std::max(0, I - static_cast<int>(this->u_band)); J <= std::min(static_cast<int>(this->n_rows) - 1, I + static_cast<int>(this->l_band)); ++J) s_mat.at(J, I) = static_cast<float>(this->at(J, I));

        device_alloc(csr_form<float, int>(s_mat));
    }

    const size_t n_rhs = sizeof(float) * B.n_elem;

    void* d_b = nullptr;
    void* d_x = nullptr;

    cudaMalloc(&d_b, n_rhs);
    cudaMalloc(&d_x, n_rhs);

    auto INFO = this->mixed_trs(X, std::move(B), [&](fmat& residual) {
        cudaMemcpyAsync(d_b, residual.memptr(), n_rhs, cudaMemcpyHostToDevice, stream);

        int singularity;

        auto code = 0;
        for(auto I = 0llu; I < residual.n_elem; I += residual.n_rows) code += cusolverSpScsrlsvqr(handle, static_cast<int>(this->n_rows), static_cast<int>(this->s_mat.n_elem), descr, (float*)d_val_idx, (int*)d_row_ptr, (int*)d_col_idx, (float*)d_b + I, static_cast<float>(this->setting.tolerance), 3, (float*)d_x + I, &singularity);

        cudaMemcpyAsync(residual.memptr(), d_x, n_rhs, cudaMemcpyDeviceToHost, stream);

        cudaDeviceSynchronize();

        return code;
    });

    if(0 != INFO)
        suanpan_error("Error code {} received, the matrix is probably singular.\n", INFO);

    return INFO;
}

#endif

#endif

//! @}
