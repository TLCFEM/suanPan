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

#include <cusolverSp.h>
#include <cusparse.h>

template<typename T> class SparseMatCUDA final : public SparseMat<T> {
	cusolverSpHandle_t handle = nullptr;
	cusparseMatDescr_t descrA = nullptr;

	void cu_create();
	void cu_destory() const;
public:
	SparseMatCUDA();
	SparseMatCUDA(uword, uword, uword = 0);
	~SparseMatCUDA() override;

	unique_ptr<MetaMat<T>> make_copy() override;

	int solve(Mat<T>&, const Mat<T>&) override;
};

template<typename T> void SparseMatCUDA<T>::cu_create() {
	cusolverSpCreate(&handle);
	cusparseCreateMatDescr(&descrA);
	cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO);
	cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
}

template<typename T> void SparseMatCUDA<T>::cu_destory() const {
	cusparseDestroyMatDescr(descrA);
	cusolverSpDestroy(handle);
}

template<typename T> SparseMatCUDA<T>::SparseMatCUDA()
	: SparseMat<T>() { cu_create(); }

template<typename T> SparseMatCUDA<T>::SparseMatCUDA(const uword in_row, const uword in_col, const uword in_elem)
	: SparseMat<T>(in_row, in_col, in_elem) { cu_create(); }

template<typename T> SparseMatCUDA<T>::~SparseMatCUDA() { cu_destory(); }

template<typename T> unique_ptr<MetaMat<T>> SparseMatCUDA<T>::make_copy() { return make_unique<SparseMatCUDA<T>>(*this); }

template<typename T> int SparseMatCUDA<T>::solve(Mat<T>& X, const Mat<T>& B) {
	csr_form<T, int> csr_mat(this->triplet_mat);

	X.set_size(size(B));

	int singularity;

	if(cusolverStatus_t solver_info; std::is_same<T, float>::value) {
		using E = float;
		for(auto I = 0; I < B.n_cols; ++I) {
			solver_info = cusolverSpScsrlsvluHost(handle, csr_mat.n_rows, csr_mat.c_size, descrA, (E*)csr_mat.val_idx, csr_mat.row_ptr, csr_mat.col_idx, (E*)B.colptr(I), this->tolerance, 2, (E*)X.colptr(I), &singularity);
			if(CUSOLVER_STATUS_SUCCESS != solver_info) {
				suanpan_error("error code %u returned during CUDA solving.\n", static_cast<unsigned>(solver_info));
				return SUANPAN_FAIL;
			}
		}
	}
	else if(std::is_same<T, double>::value) {
		using E = double;
		for(auto I = 0; I < B.n_cols; ++I) {
			solver_info = cusolverSpDcsrlsvluHost(handle, csr_mat.n_rows, csr_mat.c_size, descrA, (E*)csr_mat.val_idx, csr_mat.row_ptr, csr_mat.col_idx, (E*)B.colptr(I), this->tolerance, 2, (E*)X.colptr(I), &singularity);
			if(CUSOLVER_STATUS_SUCCESS != solver_info) {
				suanpan_error("error code %u returned during CUDA solving.\n", static_cast<unsigned>(solver_info));
				return SUANPAN_FAIL;
			}
		}
	}

	return SUANPAN_SUCCESS;
}

#endif

#endif

//! @}
