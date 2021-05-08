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

#include "SparseMat.hpp"
#include <mumps/dmumps_c.h>

template<typename T> class SparseMatBaseMUMPS : public SparseMat<T> {
	const int sym = 0;

	DMUMPS_STRUC_C mumps_job{sym, 1, -1, -987654};

	s32_vec l_irn, l_jrn;

	void alloc();
	void dealloc();
	void run();
public:
	SparseMatBaseMUMPS(uword, uword, uword, int);
	SparseMatBaseMUMPS(const SparseMatBaseMUMPS&);
	SparseMatBaseMUMPS(SparseMatBaseMUMPS&&) noexcept = delete;
	SparseMatBaseMUMPS& operator=(const SparseMatBaseMUMPS&) = delete;
	SparseMatBaseMUMPS& operator=(SparseMatBaseMUMPS&&) noexcept = delete;
	~SparseMatBaseMUMPS() override;

	void zeros() override;

	int solve(Mat<T>&, Mat<T>&&) override;
	int solve(Mat<T>&, const Mat<T>&) override;

	[[nodiscard]] int sign_det() const override;
};

template<typename T> void SparseMatBaseMUMPS<T>::alloc() {
	if(this->factored) return;

	dealloc();

	this->factored = true;

	this->triplet_mat.csc_condense();

	mumps_job.job = -1;
	dmumps_c(&mumps_job);

	mumps_job.n = static_cast<int>(this->triplet_mat.n_rows);
	mumps_job.nnz = static_cast<int64_t>(this->triplet_mat.n_elem);

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

	mumps_job.icntl[0] = -1;
	mumps_job.icntl[1] = -1;
	mumps_job.icntl[2] = -1;
	mumps_job.icntl[3] = 0;
	mumps_job.icntl[19] = 0; // dense rhs
	mumps_job.icntl[32] = 1; // determinant

	mumps_job.job = 4;
	dmumps_c(&mumps_job);
}

template<typename T> void SparseMatBaseMUMPS<T>::dealloc() {
	if(3 != mumps_job.job) return;
	mumps_job.job = -2;
	dmumps_c(&mumps_job);
}

template<typename T> void SparseMatBaseMUMPS<T>::run() {
	mumps_job.job = 3;
	dmumps_c(&mumps_job);
}

template<typename T> SparseMatBaseMUMPS<T>::SparseMatBaseMUMPS(const uword in_row, const uword in_col, const uword in_elem, const int in_sym)
	: SparseMat<T>(in_row, in_col, in_elem)
	, sym(in_sym) {}

template<typename T> SparseMatBaseMUMPS<T>::SparseMatBaseMUMPS(const SparseMatBaseMUMPS& other)
	: SparseMat<T>(other)
	, sym(other.sym)
	, mumps_job{other.sym, 1, -1, -987654}
	, l_irn(other.l_irn)
	, l_jrn(other.l_jrn) {}

template<typename T> SparseMatBaseMUMPS<T>::~SparseMatBaseMUMPS() { dealloc(); }

template<typename T> void SparseMatBaseMUMPS<T>::zeros() {
	SparseMat<T>::zeros();
	dealloc();
}

template<typename T> int SparseMatBaseMUMPS<T>::solve(Mat<T>& X, Mat<T>&& B) {
	alloc();

	mumps_job.rhs = B.memptr();
	mumps_job.lrhs = static_cast<int>(B.n_rows);
	mumps_job.nrhs = static_cast<int>(B.n_cols);

	run();

	X = std::move(B);

	return mumps_job.info[0];
}

template<typename T> int SparseMatBaseMUMPS<T>::solve(Mat<T>& X, const Mat<T>& B) {
	alloc();

	X = B;

	mumps_job.rhs = X.memptr();
	mumps_job.lrhs = static_cast<int>(X.n_rows);
	mumps_job.nrhs = static_cast<int>(X.n_cols);

	run();

	return mumps_job.info[0];
}

template<typename T> int SparseMatBaseMUMPS<T>::sign_det() const { return mumps_job.rinfog[11] < 0. ? -1 : 1; }

template<typename T> class SparseMatMUMPS final : public SparseMatBaseMUMPS<T> {
public:
	SparseMatMUMPS(uword, uword, uword = 0);

	unique_ptr<MetaMat<T>> make_copy() override;
};

template<typename T> SparseMatMUMPS<T>::SparseMatMUMPS(const uword in_row, const uword in_col, const uword in_elem)
	: SparseMatBaseMUMPS<T>(in_row, in_col, in_elem, 0) {}

template<typename T> unique_ptr<MetaMat<T>> SparseMatMUMPS<T>::make_copy() { return std::make_unique<SparseMatMUMPS<T>>(*this); }

template<typename T> class SparseSymmMatMUMPS final : public SparseMatBaseMUMPS<T> {
public:
	SparseSymmMatMUMPS(uword, uword, uword = 0);

	unique_ptr<MetaMat<T>> make_copy() override;
};

template<typename T> SparseSymmMatMUMPS<T>::SparseSymmMatMUMPS(const uword in_row, const uword in_col, const uword in_elem)
	: SparseMatBaseMUMPS<T>(in_row, in_col, in_elem, 0) {}

template<typename T> unique_ptr<MetaMat<T>> SparseSymmMatMUMPS<T>::make_copy() { return std::make_unique<SparseSymmMatMUMPS<T>>(*this); }

#endif

//! @}
