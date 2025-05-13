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

#include "../SparseMat.hpp"
#include "../csr_form.hpp"

#include <cusolverSp.h>
#include <cusparse.h>

template<sp_d T> class SparseMatCUDA final : public SparseMat<T> {
    cusolverSpHandle_t handle = nullptr;
    cudaStream_t stream = nullptr;
    cusparseMatDescr_t descr = nullptr;

    class cuda_ptr {
        void* ptr{};

    public:
        size_t size{};

        cuda_ptr(const size_t in_size = 0)
            : size(in_size) {
            if(size > 0) cudaMalloc(&ptr, size);
        }
        cuda_ptr(const cuda_ptr& other)
            : cuda_ptr(other.size) {}
        cuda_ptr(cuda_ptr&&) = delete;
        cuda_ptr& operator=(const cuda_ptr&) = delete;
        cuda_ptr& operator=(cuda_ptr&& other) noexcept {
            if(this != &other) {
                cudaFree(ptr);
                ptr = other.ptr;
                size = other.size;
                other.ptr = nullptr;
                other.size = 0;
            }
            return *this;
        }
        ~cuda_ptr() { cudaFree(ptr); }

        auto operator&() { return ptr; }

        auto copy_to(void* dest, cudaStream_t s) { return cudaMemcpyAsync(dest, ptr, size, cudaMemcpyDeviceToHost, s); }

        auto copy_from(const void* src, cudaStream_t s) { return cudaMemcpyAsync(ptr, src, size, cudaMemcpyHostToDevice, s); }
    };

    cuda_ptr d_val_idx{}, d_col_idx{}, d_row_ptr{};

    void acquire() {
        cusolverSpCreate(&handle);
        cudaStreamCreate(&stream);
        cusolverSpSetStream(handle, stream);
        cusparseCreateMatDescr(&descr);
        cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
        d_row_ptr = cuda_ptr(sizeof(int) * (this->n_rows + 1));
    }

    void release() const {
        if(handle) cusolverSpDestroy(handle);
        if(stream) cudaStreamDestroy(stream);
        if(descr) cusparseDestroyMatDescr(descr);
    }

    template<sp_d ET> void device_alloc(csr_form<ET, int>&& csr_mat) {
        d_val_idx = cuda_ptr(sizeof(ET) * csr_mat.n_elem);
        d_col_idx = cuda_ptr(sizeof(int) * csr_mat.n_elem);

        d_val_idx.copy_from(csr_mat.val_mem(), stream);
        d_col_idx.copy_from(csr_mat.col_mem(), stream);
        d_row_ptr.copy_from(csr_mat.row_mem(), stream);
    }

protected:
    using SparseMat<T>::direct_solve;

    int direct_solve(Mat<T>&, const Mat<T>&) override;

public:
    SparseMatCUDA(const uword in_row, const uword in_col, const uword in_elem = 0)
        : SparseMat<T>(in_row, in_col, in_elem) { acquire(); }

    SparseMatCUDA(const SparseMatCUDA& other)
        : SparseMat<T>(other) { acquire(); }

    SparseMatCUDA(SparseMatCUDA&&) = delete;
    SparseMatCUDA& operator=(const SparseMatCUDA&) = delete;
    SparseMatCUDA& operator=(SparseMatCUDA&&) = delete;

    ~SparseMatCUDA() override { release(); }

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SparseMatCUDA>(*this); }
};

template<sp_d T> int SparseMatCUDA<T>::direct_solve(Mat<T>& X, const Mat<T>& B) {
    if(!this->factored) {
        std::is_same_v<T, float> || Precision::MIXED == this->setting.precision ? device_alloc(csr_form<float, int>(this->triplet_mat)) : device_alloc(csr_form<double, int>(this->triplet_mat));

        this->factored = true;
    }

    const size_t n_rhs = (std::is_same_v<T, float> || Precision::MIXED == this->setting.precision ? sizeof(float) : sizeof(double)) * B.n_elem;

    cuda_ptr d_b{n_rhs}, d_x{n_rhs};

    int singularity;
    auto code = 0;

    if constexpr(std::is_same_v<T, float>) {
        d_b.copy_from(B.memptr(), stream);

        for(auto I = 0llu; I < B.n_elem; I += B.n_rows) code += cusolverSpScsrlsvqr(handle, int(this->n_rows), int(d_val_idx.size), descr, (float*)&d_val_idx, (int*)&d_row_ptr, (int*)&d_col_idx, (float*)&d_b + I, float(this->setting.tolerance), 3, (float*)&d_x + I, &singularity);

        X.set_size(arma::size(B));

        d_x.copy_to(X.memptr(), stream);

        cudaDeviceSynchronize();
    }
    else if(Precision::FULL == this->setting.precision) {
        d_b.copy_from(B.memptr(), stream);

        for(auto I = 0llu; I < B.n_elem; I += B.n_rows) code += cusolverSpDcsrlsvqr(handle, int(this->n_rows), int(d_val_idx.size), descr, (double*)&d_val_idx, (int*)&d_row_ptr, (int*)&d_col_idx, (double*)&d_b + I, this->setting.tolerance, 3, (double*)&d_x + I, &singularity);

        X.set_size(arma::size(B));

        d_x.copy_to(X.memptr(), stream);

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

            d_b.copy_from(residual.memptr(), stream);

            code = 0;
            for(auto I = 0llu; I < B.n_elem; I += B.n_rows) code += cusolverSpScsrlsvqr(handle, int(this->n_rows), int(d_val_idx.size), descr, (float*)&d_val_idx, (int*)&d_row_ptr, (int*)&d_col_idx, (float*)&d_b + I, float(this->setting.tolerance), 3, (float*)&d_x + I, &singularity);
            if(0 != code) break;

            d_x.copy_to(residual.memptr(), stream);

            cudaDeviceSynchronize();

            const mat incre = multiplier * conv_to<mat>::from(residual);

            X += incre;

            suanpan_debug("Mixed precision algorithm multiplier: {:.5E}.\n", multiplier = arma::norm(full_residual -= this->operator*(incre)));
        }
    }

    return 0 == code ? SUANPAN_SUCCESS : SUANPAN_FAIL;
}

#endif

//! @}
