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
 * @class SparseMatSuperLU
 * @brief A SparseMatSuperLU class that holds matrices.
 *
 * @author tlc
 * @date 14/08/2020
 * @version 0.1.0
 * @file SparseMatSuperLU.hpp
 * @addtogroup MetaMat
 * @{
 */

#ifndef SPARSEMATSUPERLU_HPP
#define SPARSEMATSUPERLU_HPP

#include "SparseMat.hpp"

#ifdef SUANPAN_SUPERLUMT

extern int SUANPAN_NUM_THREADS;

typedef enum { SLU_NC, SLU_NCP, SLU_NR, SLU_SC, SLU_SCP, SLU_SR, SLU_DN, SLU_NR_loc } Stype_t;

typedef enum { SLU_S, SLU_D, SLU_C, SLU_Z } Dtype_t;

typedef enum { SLU_GE, SLU_TRLU, SLU_TRUU, SLU_TRL, SLU_TRU, SLU_SYL, SLU_SYU, SLU_HEL, SLU_HEU } Mtype_t;

typedef struct {
	Stype_t Stype;
	Dtype_t Dtype;
	Mtype_t Mtype;
	int nrow;
	int ncol;
	void* Store;
} SuperMatrix;

extern "C" {
void dCreate_CompCol_Matrix(SuperMatrix*, int, int, int, double*, int*, int*, Stype_t, Dtype_t, Mtype_t);
void dCreate_Dense_Matrix(SuperMatrix*, int, int, double*, int, Stype_t, Dtype_t, Mtype_t);
void pdgssv(int, SuperMatrix*, int*, int*, SuperMatrix*, SuperMatrix*, SuperMatrix*, int*);
void Destroy_SuperMatrix_Store(SuperMatrix*);
void Destroy_CompCol_NCP(SuperMatrix*);
void Destroy_SuperNode_SCP(SuperMatrix*);
void get_perm_c(int, SuperMatrix*, int*);
}

#endif

template<typename T> class SparseMatSuperLU final : public SparseMat<T> {
#ifdef SUANPAN_SUPERLUMT
	podarray<int> perm_r, perm_c;
#endif
public:
	using SparseMat<T>::SparseMat;

	unique_ptr<MetaMat<T>> make_copy() override;

	int solve(Mat<T>&, const Mat<T>&) override;
};

template<typename T> unique_ptr<MetaMat<T>> SparseMatSuperLU<T>::make_copy() { return std::make_unique<SparseMatSuperLU<T>>(*this); }

template<typename T> int SparseMatSuperLU<T>::solve(Mat<T>& out_mat, const Mat<T>& in_mat) {
#ifdef SUANPAN_SUPERLUMT
	out_mat = in_mat;

	const csc_form<double, int> csc_mat(this->triplet_mat);

	SuperMatrix A, L, U, B;

	dCreate_CompCol_Matrix(&A, csc_mat.n_rows, csc_mat.n_cols, csc_mat.n_elem, csc_mat.val_idx, csc_mat.row_idx, csc_mat.col_ptr, SLU_NC, SLU_D, SLU_GE);
	dCreate_Dense_Matrix(&B, static_cast<int>(out_mat.n_rows), static_cast<int>(out_mat.n_cols), out_mat.memptr(), static_cast<int>(out_mat.n_rows), SLU_DN, SLU_D, SLU_GE);

	perm_r.set_size(csc_mat.n_rows);
	perm_c.set_size(csc_mat.n_cols);

	get_perm_c(1, &A, perm_c.memptr());

	auto flag = 0;

	pdgssv(SUANPAN_NUM_THREADS, &A, perm_c.memptr(), perm_r.memptr(), &L, &U, &B, &flag);

	Destroy_SuperMatrix_Store(&A);
	Destroy_SuperMatrix_Store(&B);
	Destroy_SuperNode_SCP(&L);
	Destroy_CompCol_NCP(&U);

	return flag;
#else
	csc_form<T, uword> csc_mat(this->triplet_mat);

	const uvec row_idx(csc_mat.row_idx, csc_mat.n_elem, false, false);
	const uvec col_ptr(csc_mat.col_ptr, csc_mat.n_cols + 1, false, false);
	const Col<T> val_idx(csc_mat.val_idx, csc_mat.n_elem, false, false);

	return spsolve(out_mat, SpMat<T>(row_idx, col_ptr, val_idx, csc_mat.n_rows, csc_mat.n_cols), in_mat) ? SUANPAN_SUCCESS : SUANPAN_FAIL;
#endif
}

#endif

//! @}
