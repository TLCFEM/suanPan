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

#ifndef SPARSEMATCUDA_HPP
#define SPARSEMATCUDA_HPP

#ifdef SUANPAN_CUDA

#include <cusolverSp.h>
#include <cusparse.h>

template<typename T> class SparseMatCUDA final : public SparseMat<T> {
	static void cu_destory(cusolverSpHandle_t);
public:
	using SparseMat<T>::tolerance;
	using SparseMat<T>::triplet_mat;

	SparseMatCUDA();
	SparseMatCUDA(uword, uword, uword = 0);

	unique_ptr<MetaMat<T>> make_copy() override;

	int solve(Mat<T>&, const Mat<T>&) override;
};

template<typename T> void SparseMatCUDA<T>::cu_destory(const cusolverSpHandle_t handle) { if(const auto solver_info = cusolverSpDestroy(handle); CUSOLVER_STATUS_SUCCESS != solver_info) suanpan_error("error code %u returned by CUDA.\n", static_cast<unsigned>(solver_info)); }

template<typename T> SparseMatCUDA<T>::SparseMatCUDA()
	: SparseMat<T>() {}

template<typename T> SparseMatCUDA<T>::SparseMatCUDA(const uword in_row, const uword in_col, const uword in_elem)
	: SparseMat<T>(in_row, in_col, in_elem) {}

template<typename T> unique_ptr<MetaMat<T>> SparseMatCUDA<T>::make_copy() { return make_unique<SparseMatCUDA<T>>(*this); }

template<typename T> int SparseMatCUDA<T>::solve(Mat<T>& out_mat, const Mat<T>& in_mat) {
	cusolverSpHandle_t handle = nullptr;
	auto solver_info = cusolverSpCreate(&handle);

	if(CUSOLVER_STATUS_SUCCESS != solver_info) {
		suanpan_error("error code %u returned during initialisation.\n", static_cast<unsigned>(solver_info));
		cu_destory(handle);
		return SUANPAN_FAIL;
	}

	cusparseMatDescr_t descrA = nullptr;
	if(const auto sparse_info = cusparseCreateMatDescr(&descrA); CUSPARSE_STATUS_SUCCESS != sparse_info) {
		suanpan_error("error code %u returned during initialisation.\n", static_cast<unsigned>(sparse_info));
		cu_destory(handle);
		return SUANPAN_FAIL;
	}

	csr_form<T, int> csr_mat(triplet_mat);

	out_mat.set_size(size(in_mat));

	int singularity;

	solver_info = cusolverSpDcsrlsvluHost(handle, csr_mat.n_rows, csr_mat.c_size, descrA, csr_mat.val_idx, csr_mat.row_ptr, csr_mat.col_idx, in_mat.memptr(), tolerance, 2, out_mat.memptr(), &singularity);

	if(-1 == singularity) singularity = SUANPAN_SUCCESS;
	else {
		suanpan_error("the matrix is singular.\n");
		singularity = SUANPAN_FAIL;
	}

	if(CUSOLVER_STATUS_SUCCESS != solver_info) {
		suanpan_error("error code %u returned during CUDA solving.\n", static_cast<unsigned>(solver_info));
		singularity = SUANPAN_FAIL;
	}

	cu_destory(handle);

	return singularity;
}

#endif

#endif

//! @}
