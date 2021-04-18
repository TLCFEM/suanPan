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
 * @class SparseMatMUMPS
 * @brief A SparseMatMUMPS class that holds matrices.
 *
 * * MUMPS uses int.
 *
 * @author tlc
 * @date 14/08/2020
 * @version 0.1.0
 * @file SparseMatMUMPS.hpp
 * @addtogroup MetaMat
 * @{
 */

#ifndef SPARSEMATMUMPS_HPP
#define SPARSEMATMUMPS_HPP

#include <mumps/dmumps_c.h>

template<typename T> class SparseMatMUMPS final : public SparseMat<T> {
	DMUMPS_STRUC_C mumps_job{};

	s32_vec l_irn, l_jrn;
public:
	using SparseMat<T>::SparseMat;

	unique_ptr<MetaMat<T>> make_copy() override;

	int solve(Mat<T>&, const Mat<T>&) override;

	[[nodiscard]] int sign_det() const override;
};

template<typename T> unique_ptr<MetaMat<T>> SparseMatMUMPS<T>::make_copy() { return make_unique<SparseMatMUMPS<T>>(*this); }

template<typename T> int SparseMatMUMPS<T>::solve(Mat<T>& X, const Mat<T>& B) {
	this->triplet_mat.csc_condense();

	X = B;

	mumps_job.comm_fortran = -987654;
	mumps_job.par = 1;
	mumps_job.sym = 0;
	mumps_job.job = -1;
	dmumps_c(&mumps_job);

	mumps_job.n = static_cast<int>(this->triplet_mat.n_rows);
	mumps_job.nnz = static_cast<int64_t>(this->triplet_mat.c_size);
	mumps_job.nrhs = static_cast<int>(B.n_cols);
	mumps_job.lrhs = static_cast<int>(B.n_rows);

	l_irn.set_size(mumps_job.nnz);
	l_jrn.set_size(mumps_job.nnz);

#ifdef SUANPAN_MT
	tbb::parallel_for(0, static_cast<int>(mumps_job.nnz), [&](const int I) {
		l_irn[I] = static_cast<int>(this->triplet_mat.row_idx[I] + 1);
		l_jrn[I] = static_cast<int>(this->triplet_mat.col_idx[I] + 1);
	});
#else
	for(auto I = 0; I < mumps_job.nnz; ++I) {
		l_irn[I] = static_cast<int>(this->triplet_mat.row_idx[I] + 1);
		l_jrn[I] = static_cast<int>(this->triplet_mat.col_idx[I] + 1);
	}
#endif

	mumps_job.irn = l_irn.memptr();
	mumps_job.jcn = l_jrn.memptr();
	mumps_job.a = this->triplet_mat.val_idx;
	mumps_job.rhs = X.memptr();

	mumps_job.icntl[0] = -1;
	mumps_job.icntl[1] = -1;
	mumps_job.icntl[2] = -1;
	mumps_job.icntl[3] = 0;
	mumps_job.icntl[32] = 1; // determinant

	mumps_job.job = 6;
	dmumps_c(&mumps_job);

	mumps_job.job = -2;
	dmumps_c(&mumps_job);

	return mumps_job.info[0];
}

template<typename T> int SparseMatMUMPS<T>::sign_det() const { return mumps_job.rinfog[11] < 0. ? -1 : 1; }

#endif

//! @}
