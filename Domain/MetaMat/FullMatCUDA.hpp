/*******************************************************************************
 * Copyright (C) 2017-2021 Theodore Chang
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

#include "FullMat.hpp"
#include <cublas_v2.h>
#include <cusolverDn.h>
#include <cuda_runtime.h>

template<typename T> class FullMatCUDA final : public FullMat<T> {
	cusolverDnHandle_t handle = nullptr;
	cublasHandle_t cublasHandle = nullptr;
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

template<typename T> void FullMatCUDA<T>::acquire() {
	cusolverDnCreate(&handle);
	cublasCreate_v2(&cublasHandle);
	cudaStreamCreate(&stream);
	cusolverDnSetStream(handle, stream);
	cublasSetStream_v2(cublasHandle, stream);

	cudaMalloc(&info, sizeof(int));
	cudaMemset(info, 0, sizeof(int));
	cudaMalloc(&ipiv, sizeof(int) * this->n_rows);

	if(int bufferSize = 0; std::is_same<T, float>::value || Precision::MIXED == this->precision) {
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

template<typename T> void FullMatCUDA<T>::release() const {
	if(handle) cusolverDnDestroy(handle);
	if(cublasHandle) cublasDestroy_v2(cublasHandle);
	if(stream) cudaStreamDestroy(stream);

	if(info) cudaFree(info);
	if(d_A) cudaFree(d_A);
	if(buffer) cudaFree(buffer);
	if(ipiv) cudaFree(ipiv);
}

template<typename T> FullMatCUDA<T>::FullMatCUDA(const uword in_rows, const uword in_cols)
	: FullMat<T>(in_rows, in_cols) { acquire(); }

template<typename T> FullMatCUDA<T>::FullMatCUDA(const FullMatCUDA& other)
	: FullMat<T>(other) { acquire(); }

template<typename T> FullMatCUDA<T>::~FullMatCUDA() { release(); }

template<typename T> unique_ptr<MetaMat<T>> FullMatCUDA<T>::make_copy() { return make_unique<FullMatCUDA<T>>(*this); }

template<typename T> int FullMatCUDA<T>::solve(Mat<T>& X, Mat<T>&& B) { return solve(X, B); }

template<typename T> int FullMatCUDA<T>::solve(Mat<T>& X, const Mat<T>& B) {
	if(std::is_same<T, float>::value) {
		// pure float
		if(!this->factored) {
			cudaMemcpy(d_A, this->memptr(), sizeof(float) * this->n_elem, cudaMemcpyHostToDevice);

			cusolverDnSgetrf(handle, this->n_rows, this->n_cols, (float*)d_A, this->n_rows, (float*)buffer, ipiv, info);

			int h_info = 0;
			cudaMemcpy(&h_info, info, sizeof(int), cudaMemcpyDeviceToHost);
			if(0 != h_info) return h_info;

			this->factored = true;
		}

		const size_t btye_size = sizeof(float) * B.n_elem;

		void* d_x = nullptr;
		cudaMalloc(&d_x, btye_size);
		cudaMemcpy(d_x, B.memptr(), btye_size, cudaMemcpyHostToDevice);
		cusolverDnSgetrs(handle, CUBLAS_OP_N, this->n_rows, B.n_cols, (float*)d_A, this->n_rows, ipiv, (float*)d_x, this->n_rows, info);
		cudaDeviceSynchronize();

		X.set_size(arma::size(B));

		cudaMemcpy(X.memptr(), d_x, btye_size, cudaMemcpyDeviceToHost);

		if(d_x) cudaFree(d_x);
	}
	else if(Precision::MIXED == this->precision) {
		// mixed precision
		if(!this->factored) {
			this->s_memory = this->to_float();

			cudaMemcpy(d_A, this->s_memory.memptr(), sizeof(float) * this->s_memory.n_elem, cudaMemcpyHostToDevice);

			cusolverDnSgetrf(handle, this->n_rows, this->n_cols, (float*)d_A, this->n_rows, (float*)buffer, ipiv, info);

			int h_info = 0;
			cudaMemcpy(&h_info, info, sizeof(int), cudaMemcpyDeviceToHost);
			if(0 != h_info) return h_info;

			this->factored = true;
		}

		const size_t btye_size = sizeof(float) * B.n_elem;

		void* d_x = nullptr;
		cudaMalloc(&d_x, btye_size);

		X = arma::zeros(B.n_rows, B.n_cols);

		mat full_residual = B;

		auto multiplier = 1.;

		auto counter = 0;
		while(++counter < 20) {
			auto residual = conv_to<fmat>::from(full_residual / multiplier);

			cudaMemcpy(d_x, residual.memptr(), btye_size, cudaMemcpyHostToDevice);
			cusolverDnSgetrs(handle, CUBLAS_OP_N, this->n_rows, B.n_cols, (float*)d_A, this->n_rows, ipiv, (float*)d_x, this->n_rows, info);
			cudaDeviceSynchronize();

			cudaMemcpy(residual.memptr(), d_x, btye_size, cudaMemcpyDeviceToHost);

			const vec incre = multiplier * conv_to<mat>::from(residual);

			X += incre;

			suanpan_debug("mixed precision algorithm multiplier: %.5E\n", multiplier = arma::norm(full_residual -= this->operator*(incre)));

			if(multiplier < this->tolerance) break;
		}

		if(d_x) cudaFree(d_x);
	}
	else {
		// pure double
		if(!this->factored) {
			cudaMemcpy(d_A, this->memptr(), sizeof(double) * this->n_elem, cudaMemcpyHostToDevice);

			cusolverDnDgetrf(handle, this->n_rows, this->n_cols, (double*)d_A, this->n_rows, (double*)buffer, ipiv, info);

			int h_info = 0;
			cudaMemcpy(&h_info, info, sizeof(int), cudaMemcpyDeviceToHost);
			if(0 != h_info) return h_info;

			this->factored = true;
		}

		const size_t btye_size = sizeof(double) * B.n_elem;

		void* d_x = nullptr;
		cudaMalloc(&d_x, btye_size);
		cudaMemcpy(d_x, B.memptr(), btye_size, cudaMemcpyHostToDevice);
		cusolverDnDgetrs(handle, CUBLAS_OP_N, this->n_rows, B.n_cols, (double*)d_A, this->n_rows, ipiv, (double*)d_x, this->n_rows, info);
		cudaDeviceSynchronize();

		X.set_size(arma::size(B));

		cudaMemcpy(X.memptr(), d_x, btye_size, cudaMemcpyDeviceToHost);

		if(d_x) cudaFree(d_x);
	}

	return 0;
}

#endif

#endif

//! @}
