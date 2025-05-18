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
 * todo: use cuDSS library instead
 *
 * @author tlc
 * @date 17/05/2025
 * @version 0.2.0
 * @file SparseMatCUDA.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef SPARSEMATCUDA_HPP
#define SPARSEMATCUDA_HPP

#include "../SparseMat.hpp"
#include "../csr_form.hpp"
#include "../cuda_ptr.hpp"

#include <cusolverSp.h>

template<sp_d T> class SparseMatCUDA final : public SparseMat<T> {
    cusolverSpHandle_t handle = nullptr;
    cudaStream_t stream = nullptr;
    cusparseMatDescr_t descr = nullptr;

    cuda_ptr d_val_idx{}, d_col_idx{}, d_row_ptr{};

    void init_config() {
        cusolverSpCreate(&handle);
        cudaStreamCreate(&stream);
        cusolverSpSetStream(handle, stream);
        cusparseCreateMatDescr(&descr);
        cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
        d_row_ptr = cuda_ptr(sizeof(int), static_cast<int>(this->n_rows) + 1);
        this->factored = false;
    }

    void release() const {
        if(handle) cusolverSpDestroy(handle);
        if(stream) cudaStreamDestroy(stream);
        if(descr) cusparseDestroyMatDescr(descr);
    }

    template<sp_d ET> void device_alloc(csr_form<ET, int>&& csr_mat) {
        d_val_idx = cuda_ptr(sizeof(ET), csr_mat.n_elem);
        d_col_idx = cuda_ptr(sizeof(int), csr_mat.n_elem);

        d_val_idx.copy_from(csr_mat.val_mem(), stream);
        d_col_idx.copy_from(csr_mat.col_mem(), stream);
        d_row_ptr.copy_from(csr_mat.row_mem(), stream);
    }

protected:
    using SparseMat<T>::direct_solve;

    int direct_solve(Mat<T>&, const Mat<T>&) override;

public:
    SparseMatCUDA(const uword in_row, const uword in_col, const uword in_elem = 0)
        : SparseMat<T>(in_row, in_col, in_elem) { init_config(); }

    SparseMatCUDA(const SparseMatCUDA& other)
        : SparseMat<T>(other) { init_config(); }

    SparseMatCUDA(SparseMatCUDA&&) = delete;
    SparseMatCUDA& operator=(const SparseMatCUDA&) = delete;
    SparseMatCUDA& operator=(SparseMatCUDA&&) = delete;

    ~SparseMatCUDA() override { release(); }

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SparseMatCUDA>(*this); }
};

template<sp_d T> int SparseMatCUDA<T>::direct_solve(Mat<T>& X, const Mat<T>& B) {
    const auto single_precision = std::is_same_v<T, float> || Precision::MIXED == this->setting.precision;

    if(!this->factored) {
        this->factored = true;
        single_precision ? device_alloc(csr_form<float, int>(this->triplet_mat)) : device_alloc(csr_form<double, int>(this->triplet_mat));
    }

    const size_t unit_size = single_precision ? sizeof(float) : sizeof(double);

    const cuda_ptr d_b{unit_size, static_cast<int>(B.n_elem)}, d_x{unit_size, static_cast<int>(B.n_elem)};

    int singularity;
    Col<int> code(B.n_cols);

    if constexpr(std::is_same_v<T, float>) {
        d_b.copy_from(B.memptr(), stream);

        for(auto I = 0llu, J = 0llu; I < B.n_elem; I += B.n_rows, ++J) code[J] = cusolverSpScsrlsvqr(handle, static_cast<int>(this->n_rows), d_val_idx.size, descr, d_val_idx.get<float>(), d_row_ptr.get(), d_col_idx.get(), d_b.get<float>(I), std::numeric_limits<float>::epsilon(), 3, d_x.get<float>(I), &singularity);

        if(0 == code.max()) {
            X.set_size(arma::size(B));
            d_x.copy_to(X.memptr(), stream);
        }
    }
    else if(Precision::FULL == this->setting.precision) {
        d_b.copy_from(B.memptr(), stream);

        for(auto I = 0llu, J = 0llu; I < B.n_elem; I += B.n_rows, ++J) code[J] = cusolverSpDcsrlsvqr(handle, static_cast<int>(this->n_rows), d_val_idx.size, descr, d_val_idx.get<double>(), d_row_ptr.get(), d_col_idx.get(), d_b.get<double>(I), std::numeric_limits<double>::epsilon(), 3, d_x.get<double>(I), &singularity);

        if(0 == code.max()) {
            X.set_size(arma::size(B));
            d_x.copy_to(X.memptr(), stream);
        }
    }
    else {
        X = arma::zeros(arma::size(B));

        mat full_residual = B;

        std::uint8_t counter{0};
        while(counter++ < this->setting.iterative_refinement) {
            const auto multiplier = norm(full_residual);
            if(multiplier < this->setting.tolerance) break;
            suanpan_debug("Mixed precision algorithm multiplier: {:.5E}.\n", multiplier);

            auto residual = conv_to<fmat>::from(full_residual / multiplier);
            d_b.copy_from(residual.memptr(), stream);

            for(auto I = 0llu, J = 0llu; I < B.n_elem; I += B.n_rows, ++J) code[J] = cusolverSpScsrlsvqr(handle, static_cast<int>(this->n_rows), d_val_idx.size, descr, d_val_idx.get<float>(), d_row_ptr.get(), d_col_idx.get(), d_b.get<float>(I), std::numeric_limits<float>::epsilon(), 3, d_x.get<float>(I), &singularity);
            if(0 != code.max()) break;

            d_x.copy_to(residual.memptr(), stream);
            full_residual = B - this->operator*(X += multiplier * conv_to<mat>::from(residual));
        }
    }

    return 0 == code.max() ? SUANPAN_SUCCESS : SUANPAN_FAIL;
}

#endif

//! @}
