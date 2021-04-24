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
 * @class SparseMatPARDISO
 * @brief A SparseMatPARDISO class that holds matrices.
 *
 * TODO: improve performance by storing factorization and resusing it
 *
 * @author tlc
 * @date 20/01/2021
 * @version 0.1.0
 * @file SparseMatPARDISO.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef SPARSEMATPARDISO_HPP
#define SPARSEMATPARDISO_HPP

#ifdef SUANPAN_MKL

#include "SparseMat.hpp"
#include <mkl_pardiso.h>

template<typename T> class SparseMatPARDISO final : public SparseMat<T> {
	int maxfct = 1;
	int mnum = 1;
	int mtype = 1;
	int msglvl = 0;

	void* pt[64];
public:
	using SparseMat<T>::SparseMat;

	unique_ptr<MetaMat<T>> make_copy() override;

	int solve(Mat<T>&, const Mat<T>&) override;
};

template<typename T> unique_ptr<MetaMat<T>> SparseMatPARDISO<T>::make_copy() { return std::make_unique<SparseMatPARDISO<T>>(*this); }

template<typename T> int SparseMatPARDISO<T>::solve(Mat<T>& X, const Mat<T>& B) {
	X.set_size(B.n_rows, B.n_cols);

	csr_form<T, int> csr_mat(this->triplet_mat);

	auto n = static_cast<int>(B.n_rows);
	auto nrhs = static_cast<int>(B.n_cols);
	int error;

	std::vector iparm(64, 0);

	pardisoinit(pt, &mtype, iparm.data());

	iparm[34] = 1; // zero-based indexing
	if(std::is_same<T, float>::value) iparm[27] = 1;

	auto phase = 12;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, (void*)csr_mat.val_idx, csr_mat.row_ptr, csr_mat.col_idx, nullptr, &nrhs, iparm.data(), &msglvl, nullptr, nullptr, &error);

	phase = 33;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, (void*)csr_mat.val_idx, csr_mat.row_ptr, csr_mat.col_idx, nullptr, &nrhs, iparm.data(), &msglvl, (void*)B.memptr(), (void*)X.memptr(), &error);

	phase = -1;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, nullptr, csr_mat.row_ptr, csr_mat.col_idx, nullptr, &nrhs, iparm.data(), &msglvl, nullptr, nullptr, &error);

	return 0 == error ? SUANPAN_SUCCESS : SUANPAN_FAIL;
}

#endif

#endif

//! @}
