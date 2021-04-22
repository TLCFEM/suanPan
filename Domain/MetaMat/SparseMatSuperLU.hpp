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

// ReSharper disable CppCStyleCast
#ifndef SPARSEMATSUPERLU_HPP
#define SPARSEMATSUPERLU_HPP

#include "SparseMat.hpp"
#include "csc_form.hpp"
#include <superlu-mt/superlu-mt.h>

template<typename T> class SparseMatSuperLU final : public SparseMat<T> {
	SuperMatrix A{}, L{}, U{};

	superlu_array_wrangler<int> perm_r;
	superlu_array_wrangler<int> perm_c;

	superlu::superlu_options_t options;

	superlu_stat_wrangler stat;

	void* t_val = nullptr;
	int* t_row = nullptr;
	int* t_col = nullptr;

	bool allocated = false;

#ifdef SUANPAN_SUPERLUMT
	const int ordering_num = 1;
#endif

	template<typename ET> void alloc_supermatrix(csc_form<ET, int>&&);
	void dealloc_supermatrix();

	int solve_trs(Mat<T>&, Mat<T>&&);
	int solve_trs(Mat<T>&, const Mat<T>&);
public:
	SparseMatSuperLU(uword, uword, uword = 0);
	SparseMatSuperLU(const SparseMatSuperLU&);
	SparseMatSuperLU(SparseMatSuperLU&&) noexcept = delete;
	SparseMatSuperLU& operator=(const SparseMatSuperLU&) = delete;
	SparseMatSuperLU& operator=(SparseMatSuperLU&&) noexcept = delete;
	~SparseMatSuperLU() override;

	void zeros() override;

	unique_ptr<MetaMat<T>> make_copy() override;

	int solve(Mat<T>&, Mat<T>&&) override;
	int solve(Mat<T>&, const Mat<T>&) override;
};

template<typename T> template<typename ET> void SparseMatSuperLU<T>::alloc_supermatrix(csc_form<ET, int>&& in) {
	dealloc_supermatrix();

	t_val = in.val_idx;
	t_row = in.row_idx;
	t_col = in.col_ptr;

	in.val_idx = nullptr;
	in.row_idx = nullptr;
	in.col_ptr = nullptr;

	if(std::is_same<ET, double>::value) {
		using E = double;
		dCreate_CompCol_Matrix(&A, in.n_rows, in.n_cols, in.n_elem, (E*)t_val, t_row, t_col, Stype_t::SLU_NC, Dtype_t::SLU_D, Mtype_t::SLU_GE);
	}
	else {
		using E = float;
		sCreate_CompCol_Matrix(&A, in.n_rows, in.n_cols, in.n_elem, (E*)t_val, t_row, t_col, Stype_t::SLU_NC, Dtype_t::SLU_S, Mtype_t::SLU_GE);
	}

	allocated = true;
}

template<typename T> void SparseMatSuperLU<T>::dealloc_supermatrix() {
	if(!allocated) return;

	if(std::is_same<T, float>::value || Precision::MIXED == this->precision) delete[] static_cast<float*>(t_val);
	else delete[] static_cast<double*>(t_val);

	delete[]t_row;
	delete[]t_col;

	Destroy_SuperMatrix_Store(&A);
#ifdef SUANPAN_SUPERLUMT
	Destroy_SuperNode_SCP(&L);
	Destroy_CompCol_NCP(&U);
#else
	Destroy_SuperNode_Matrix(&L);
	Destroy_CompCol_Matrix(&U);
#endif
	allocated = false;
}

template<typename T> SparseMatSuperLU<T>::SparseMatSuperLU(const uword in_row, const uword in_col, const uword in_elem)
	: SparseMat<T>(in_row, in_col, in_elem)
	, perm_r(this->n_rows + 1)
	, perm_c(this->n_cols + 1) { set_default_options(&options); }

template<typename T> SparseMatSuperLU<T>::SparseMatSuperLU(const SparseMatSuperLU& other)
	: SparseMat<T>(other)
	, perm_r(this->n_rows + 1)
	, perm_c(this->n_cols + 1) { set_default_options(&options); }

template<typename T> SparseMatSuperLU<T>::~SparseMatSuperLU() { dealloc_supermatrix(); }

template<typename T> void SparseMatSuperLU<T>::zeros() {
	SparseMat<T>::zeros();
	dealloc_supermatrix();
}

template<typename T> unique_ptr<MetaMat<T>> SparseMatSuperLU<T>::make_copy() { return std::make_unique<SparseMatSuperLU<T>>(*this); }

template<typename T> int SparseMatSuperLU<T>::solve(Mat<T>& out_mat, const Mat<T>& in_mat) {
	if(this->factored) return solve_trs(out_mat, in_mat);

	this->factored = true;

	auto flag = 0;

	if(std::is_same<T, float>::value || Precision::FULL == this->precision) {
		alloc_supermatrix(csc_form<T, int>(this->triplet_mat));

		out_mat = in_mat;

		superlu_supermatrix_wrangler B;
		sp_auxlib::wrap_to_supermatrix(B.get_ref(), out_mat);

#ifdef SUANPAN_SUPERLUMT
		get_perm_c(ordering_num, &A, perm_c.get_ptr());
		if(std::is_same<T, float>::value) psgssv(SUANPAN_NUM_THREADS, &A, perm_c.get_ptr(), perm_r.get_ptr(), &L, &U, B.get_ptr(), &flag);
		else pdgssv(SUANPAN_NUM_THREADS, &A, perm_c.get_ptr(), perm_r.get_ptr(), &L, &U, B.get_ptr(), &flag);
#else
		superlu::gssv<T>(&options, &A, perm_c.get_ptr(), perm_r.get_ptr(), &L, &U, B.get_ptr(), stat.get_ptr(), &flag);
#endif

		return flag;
	}

	alloc_supermatrix(csc_form<float, int>(this->triplet_mat));

	const fmat f_mat(arma::size(in_mat));

	superlu_supermatrix_wrangler B;
	sp_auxlib::wrap_to_supermatrix(B.get_ref(), f_mat);

#ifdef SUANPAN_SUPERLUMT
	get_perm_c(ordering_num, &A, perm_c.get_ptr());
	psgssv(SUANPAN_NUM_THREADS, &A, perm_c.get_ptr(), perm_r.get_ptr(), &L, &U, B.get_ptr(), &flag);
#else
	superlu::gssv<float>(&options, &A, perm_c.get_ptr(), perm_r.get_ptr(), &L, &U, B.get_ptr(), stat.get_ptr(), &flag);
#endif

	return 0 == flag ? solve_trs(out_mat, in_mat) : flag;
}

template<typename T> int SparseMatSuperLU<T>::solve_trs(Mat<T>& out_mat, const Mat<T>& in_mat) {
	auto flag = 0;

	if(std::is_same<T, float>::value || Precision::FULL == this->precision) {
		out_mat = in_mat;

		superlu_supermatrix_wrangler B;
		sp_auxlib::wrap_to_supermatrix(B.get_ref(), out_mat);

		superlu::gstrs<T>(options.Trans, &L, &U, perm_c.get_ptr(), perm_r.get_ptr(), B.get_ptr(), stat.get_ptr(), &flag);

		return flag;
	}

	out_mat.zeros(arma::size(in_mat));

	mat full_residual = in_mat;

	auto multiplier = 1.;

	auto counter = 0;
	while(++counter < 20) {
		auto residual = conv_to<fmat>::from(full_residual / multiplier);

		superlu_supermatrix_wrangler B;
		sp_auxlib::wrap_to_supermatrix(B.get_ref(), residual);

		superlu::gstrs<float>(options.Trans, &L, &U, perm_c.get_ptr(), perm_r.get_ptr(), B.get_ptr(), stat.get_ptr(), &flag);

		if(0 != flag) break;

		const mat incre = multiplier * conv_to<mat>::from(residual);

		out_mat += incre;

		suanpan_debug("mixed precision algorithm multiplier: %.5E\n", multiplier = norm(full_residual -= this->operator*(incre)));

		if(multiplier < this->tolerance) break;
	}

	return flag;
}

template<typename T> int SparseMatSuperLU<T>::solve(Mat<T>& out_mat, Mat<T>&& in_mat) {
	if(this->factored) return solve_trs(out_mat, std::forward<Mat<T>>(in_mat));

	this->factored = true;

	auto flag = 0;

	if(std::is_same<T, float>::value || Precision::FULL == this->precision) {
		alloc_supermatrix(csc_form<T, int>(this->triplet_mat));

		superlu_supermatrix_wrangler B;
		sp_auxlib::wrap_to_supermatrix(B.get_ref(), in_mat);

#ifdef SUANPAN_SUPERLUMT
		get_perm_c(ordering_num, &A, perm_c.get_ptr());
		if(std::is_same<T, float>::value) psgssv(SUANPAN_NUM_THREADS, &A, perm_c.get_ptr(), perm_r.get_ptr(), &L, &U, B.get_ptr(), &flag);
		else pdgssv(SUANPAN_NUM_THREADS, &A, perm_c.get_ptr(), perm_r.get_ptr(), &L, &U, B.get_ptr(), &flag);
#else
		superlu::gssv<T>(&options, &A, perm_c.get_ptr(), perm_r.get_ptr(), &L, &U, B.get_ptr(), stat.get_ptr(), &flag);
#endif

		out_mat = std::move(in_mat);

		return flag;
	}

	alloc_supermatrix(csc_form<float, int>(this->triplet_mat));

	const fmat f_mat(arma::size(in_mat));

	superlu_supermatrix_wrangler B;
	sp_auxlib::wrap_to_supermatrix(B.get_ref(), f_mat);

#ifdef SUANPAN_SUPERLUMT
	get_perm_c(ordering_num, &A, perm_c.get_ptr());
	psgssv(SUANPAN_NUM_THREADS, &A, perm_c.get_ptr(), perm_r.get_ptr(), &L, &U, B.get_ptr(), &flag);
#else
	superlu::gssv<float>(&options, &A, perm_c.get_ptr(), perm_r.get_ptr(), &L, &U, B.get_ptr(), stat.get_ptr(), &flag);
#endif

	return 0 == flag ? solve_trs(out_mat, std::forward<Mat<T>>(in_mat)) : flag;
}

template<typename T> int SparseMatSuperLU<T>::solve_trs(Mat<T>& out_mat, Mat<T>&& in_mat) {
	auto flag = 0;

	if(std::is_same<T, float>::value || Precision::FULL == this->precision) {
		superlu_supermatrix_wrangler B;
		sp_auxlib::wrap_to_supermatrix(B.get_ref(), in_mat);

		superlu::gstrs<T>(options.Trans, &L, &U, perm_c.get_ptr(), perm_r.get_ptr(), B.get_ptr(), stat.get_ptr(), &flag);

		out_mat = std::move(in_mat);

		return flag;
	}

	out_mat.zeros(arma::size(in_mat));

	auto multiplier = 1.;

	auto counter = 0;
	while(++counter < 20) {
		auto residual = conv_to<fmat>::from(in_mat / multiplier);

		superlu_supermatrix_wrangler B;
		sp_auxlib::wrap_to_supermatrix(B.get_ref(), residual);

		superlu::gstrs<float>(options.Trans, &L, &U, perm_c.get_ptr(), perm_r.get_ptr(), B.get_ptr(), stat.get_ptr(), &flag);

		if(0 != flag) break;

		const mat incre = multiplier * conv_to<mat>::from(residual);

		out_mat += incre;

		suanpan_debug("mixed precision algorithm multiplier: %.5E\n", multiplier = norm(in_mat -= this->operator*(incre)));

		if(multiplier < this->tolerance) break;
	}

	return flag;
}
#endif

//! @}
