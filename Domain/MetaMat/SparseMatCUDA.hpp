/*******************************************************************************
 * Copyright (C) 2017-2019 Theodore Chang
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
 * @author tlc
 * @date 14/01/2021
 * @version 0.1.0
 * @file SparseMatCUDA.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef SPARSEMATCUDA_HPP
#define SPARSEMATCUDA_HPP

#ifdef SUANPAN_CUDA

#include "SparseMat.hpp"
#include <cusparse.h>
#include <cusolverSp.h>

template<typename T> class SparseMatCUDA final : public SparseMat<T> {
	cusolverSpHandle_t handle = nullptr;
	cudaStream_t stream = nullptr;
	cusparseMatDescr_t descrA = nullptr;

	void* d_val_idx;
	void* d_col_idx;
	void* d_row_ptr;

	void *d_b, *d_x;

	void acquire();
	void release() const;

	void device_alloc(size_t, size_t, size_t, size_t);
	void device_dealloc() const;
public:
	SparseMatCUDA(uword, uword, uword = 0);
	SparseMatCUDA(const SparseMatCUDA&);
	SparseMatCUDA(SparseMatCUDA&&) noexcept = delete;
	SparseMatCUDA& operator=(const SparseMatCUDA&) = delete;
	SparseMatCUDA& operator=(SparseMatCUDA&&) noexcept = delete;
	~SparseMatCUDA() override;

	unique_ptr<MetaMat<T>> make_copy() override;

	int solve(Mat<T>&, const Mat<T>&) override;
};

template<typename T> void SparseMatCUDA<T>::acquire() {
	cusolverSpCreate(&handle);
	cudaStreamCreate(&stream);
	cusolverSpSetStream(handle, stream);
	cusparseCreateMatDescr(&descrA);
	cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO);
}

template<typename T> void SparseMatCUDA<T>::release() const {
	if(handle) cusolverSpDestroy(handle);
	if(stream) cudaStreamDestroy(stream);
	if(descrA) cusparseDestroyMatDescr(descrA);
}

template<typename T> void SparseMatCUDA<T>::device_alloc(const size_t n_val, const size_t n_col, const size_t n_row, const size_t n_rhs) {
	cudaMalloc(&d_val_idx, n_val);
	cudaMalloc(&d_col_idx, n_col);
	cudaMalloc(&d_row_ptr, n_row);
	cudaMalloc(&d_b, n_rhs);
	cudaMalloc(&d_x, n_rhs);
}

template<typename T> void SparseMatCUDA<T>::device_dealloc() const {
	if(d_val_idx) cudaFree(d_val_idx);
	if(d_col_idx) cudaFree(d_col_idx);
	if(d_row_ptr) cudaFree(d_row_ptr);
	if(d_b) cudaFree(d_b);
	if(d_x) cudaFree(d_x);
}

template<typename T> SparseMatCUDA<T>::SparseMatCUDA(const uword in_row, const uword in_col, const uword in_elem)
	: SparseMat<T>(in_row, in_col, in_elem) { acquire(); }

template<typename T> SparseMatCUDA<T>::SparseMatCUDA(const SparseMatCUDA& other)
	: SparseMat<T>(other) { acquire(); }

template<typename T> SparseMatCUDA<T>::~SparseMatCUDA() { release(); }

template<typename T> unique_ptr<MetaMat<T>> SparseMatCUDA<T>::make_copy() { return std::make_unique<SparseMatCUDA<T>>(*this); }

template<typename T> int SparseMatCUDA<T>::solve(Mat<T>& X, const Mat<T>& B) {
	csr_form<T, int> csr_mat(this->triplet_mat);

	auto n_val = sizeof(double) * csr_mat.n_elem;
	auto n_col = sizeof(int) * csr_mat.n_elem;
	auto n_row = sizeof(int) * (csr_mat.n_rows + 1);
	auto n_rhs = sizeof(double) * B.n_elem;

	device_alloc(n_val, n_col, n_row, n_rhs);

	cudaMemcpyAsync(d_val_idx, csr_mat.val_idx, n_val, cudaMemcpyHostToDevice, stream);
	cudaMemcpyAsync(d_col_idx, csr_mat.col_idx, n_col, cudaMemcpyHostToDevice, stream);
	cudaMemcpyAsync(d_row_ptr, csr_mat.row_ptr, n_row, cudaMemcpyHostToDevice, stream);

	cudaMemcpyAsync(d_b, B.memptr(), n_rhs, cudaMemcpyHostToDevice, stream);

	int singularity;

	cudaDeviceSynchronize();

	for(auto I = 0; I < B.n_elem; I += B.n_rows) cusolverSpDcsrlsvqr(handle, csr_mat.n_rows, csr_mat.n_elem, descrA, (double*)d_val_idx, (int*)d_row_ptr, (int*)d_col_idx, (double*)d_b + I, this->tolerance, 0, (double*)d_x + I, &singularity);

	cudaDeviceSynchronize();

	X.set_size(arma::size(B));

	cudaMemcpy(X.memptr(), d_x, n_rhs, cudaMemcpyDeviceToHost);

	device_dealloc();

	return SUANPAN_SUCCESS;
}

#endif

#endif

//! @}
