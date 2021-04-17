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

#include <mkl_pardiso.h>

template<typename T> class SparseMatPARDISO final : public SparseMat<T> {
public:
	using SparseMat<T>::SparseMat;
	using SparseMat<T>::triplet_mat;

	unique_ptr<MetaMat<T>> make_copy() override;

	int solve(Mat<T>&, const Mat<T>&) override;

	void print() override;
};

template<typename T> unique_ptr<MetaMat<T>> SparseMatPARDISO<T>::make_copy() { return make_unique<SparseMatPARDISO<T>>(*this); }

template<typename T> int SparseMatPARDISO<T>::solve(Mat<T>& out_mat, const Mat<T>& in_mat) {
	out_mat.set_size(size(in_mat));

	csr_form<T, int> csr_mat(triplet_mat);

	auto maxfct = 1;
	auto mnum = 1;
	auto mtype = 11;
	auto n = static_cast<int>(in_mat.n_rows);
	auto nrhs = static_cast<int>(in_mat.n_cols);
	auto msglvl = 0;
	int error;

	podarray<int> pt(64), iparm(64);
	const podarray<int> perm(n);

	pardisoinit((void*)pt.memptr(), &mtype, iparm.memptr());

	iparm(34) = 1; // zero-based indexing
	if(std::is_same<T, float>::value) iparm(27) = 1;

	auto phase = 13;
	pardiso((void*)pt.memptr(), &maxfct, &mnum, &mtype, &phase, &n, (void*)csr_mat.val_idx, csr_mat.row_ptr, csr_mat.col_idx, perm.mem, &nrhs, iparm.memptr(), &msglvl, (void*)in_mat.memptr(), (void*)out_mat.memptr(), &error);

	phase = -1;
	pardiso((void*)pt.memptr(), &maxfct, &mnum, &mtype, &phase, &n, (void*)csr_mat.val_idx, csr_mat.row_ptr, csr_mat.col_idx, perm.mem, &nrhs, iparm.memptr(), &msglvl, (void*)in_mat.memptr(), (void*)out_mat.memptr(), &error);

	return 0 == error ? SUANPAN_SUCCESS : SUANPAN_FAIL;
}

template<typename T> void SparseMatPARDISO<T>::print() { triplet_mat.print(); }

#endif

#endif

//! @}
