/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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
#include "SparseMat.hpp"

template<sp_d T> class SparseMatBaseMUMPS : public SparseMat<T> {
    const int sym = 0;

    DMUMPS_STRUC_C mumps_job{sym, 1, -1, -987654};

    s32_vec l_irn, l_jrn;

    int alloc();
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

    int direct_solve(Mat<T>&, Mat<T>&&) override;
    int direct_solve(Mat<T>&, const Mat<T>&) override;

    [[nodiscard]] int sign_det() const override;
};

template<sp_d T> int SparseMatBaseMUMPS<T>::alloc() {
    if(this->factored) return 0;

    dealloc();

    this->factored = true;

    this->triplet_mat.csc_condense();

    mumps_job.job = -1;
    dmumps_c(&mumps_job);

    mumps_job.n = static_cast<int>(this->triplet_mat.n_rows);
    mumps_job.nnz = static_cast<int64_t>(this->triplet_mat.n_elem);

    l_irn.set_size(mumps_job.nnz);
    l_jrn.set_size(mumps_job.nnz);

    suanpan_for(0, static_cast<int>(mumps_job.nnz), [&](const int I) {
        l_irn[I] = static_cast<int>(this->triplet_mat.row(I) + 1);
        l_jrn[I] = static_cast<int>(this->triplet_mat.col(I) + 1);
    });

    mumps_job.irn = l_irn.memptr();
    mumps_job.jcn = l_jrn.memptr();
    mumps_job.a = this->triplet_mat.val_mem();

    mumps_job.icntl[0] = -1;
    mumps_job.icntl[1] = -1;
    mumps_job.icntl[2] = -1;
    mumps_job.icntl[3] = 0;
    mumps_job.icntl[13] = 100;
    mumps_job.icntl[19] = 0; // dense rhs
    mumps_job.icntl[32] = 1; // determinant
    mumps_job.icntl[34] = 1; // BLR

    mumps_job.job = 4;
    dmumps_c(&mumps_job);

    if(0 != mumps_job.info[0]) suanpan_error("factorization fails with code %d.\n", mumps_job.info[0]);

    return mumps_job.info[0];
}

template<sp_d T> void SparseMatBaseMUMPS<T>::dealloc() {
    if(3 != mumps_job.job) return;
    mumps_job.job = -2;
    dmumps_c(&mumps_job);
}

template<sp_d T> void SparseMatBaseMUMPS<T>::run() {
    mumps_job.job = 3;
    dmumps_c(&mumps_job);
}

template<sp_d T> SparseMatBaseMUMPS<T>::SparseMatBaseMUMPS(const uword in_row, const uword in_col, const uword in_elem, const int in_sym)
    : SparseMat<T>(in_row, in_col, in_elem)
    , sym(in_sym) {}

template<sp_d T> SparseMatBaseMUMPS<T>::SparseMatBaseMUMPS(const SparseMatBaseMUMPS& other)
    : SparseMat<T>(other)
    , sym(other.sym)
    , mumps_job{other.sym, 1, -1, -987654}
    , l_irn(other.l_irn)
    , l_jrn(other.l_jrn) {}

template<sp_d T> SparseMatBaseMUMPS<T>::~SparseMatBaseMUMPS() { dealloc(); }

template<sp_d T> void SparseMatBaseMUMPS<T>::zeros() {
    SparseMat<T>::zeros();
    dealloc();
}

template<sp_d T> int SparseMatBaseMUMPS<T>::direct_solve(Mat<T>& X, Mat<T>&& B) {
    if(const auto code = alloc(); 0 != code) return code;

    mumps_job.rhs = B.memptr();
    mumps_job.lrhs = static_cast<int>(B.n_rows);
    mumps_job.nrhs = static_cast<int>(B.n_cols);

    run();

    X = std::move(B);

    return mumps_job.info[0];
}

template<sp_d T> int SparseMatBaseMUMPS<T>::direct_solve(Mat<T>& X, const Mat<T>& B) {
    if(const auto code = alloc(); 0 != code) return code;

    X = B;

    mumps_job.rhs = X.memptr();
    mumps_job.lrhs = static_cast<int>(X.n_rows);
    mumps_job.nrhs = static_cast<int>(X.n_cols);

    run();

    return mumps_job.info[0];
}

template<sp_d T> int SparseMatBaseMUMPS<T>::sign_det() const {
    if(IterativeSolver::NONE != this->setting.iterative_solver) throw invalid_argument("analysis requires the sign of determinant but iterative solver does not support it");
    return mumps_job.rinfog[11] < 0. ? -1 : 1;
}

template<sp_d T> class SparseMatMUMPS final : public SparseMatBaseMUMPS<T> {
public:
    SparseMatMUMPS(uword, uword, uword = 0);

    unique_ptr<MetaMat<T>> make_copy() override;
};

template<sp_d T> SparseMatMUMPS<T>::SparseMatMUMPS(const uword in_row, const uword in_col, const uword in_elem)
    : SparseMatBaseMUMPS<T>(in_row, in_col, in_elem, 0) {}

template<sp_d T> unique_ptr<MetaMat<T>> SparseMatMUMPS<T>::make_copy() { return std::make_unique<SparseMatMUMPS<T>>(*this); }

template<sp_d T> class SparseSymmMatMUMPS final : public SparseMatBaseMUMPS<T> {
public:
    SparseSymmMatMUMPS(uword, uword, uword = 0);

    unique_ptr<MetaMat<T>> make_copy() override;
};

template<sp_d T> SparseSymmMatMUMPS<T>::SparseSymmMatMUMPS(const uword in_row, const uword in_col, const uword in_elem)
    : SparseMatBaseMUMPS<T>(in_row, in_col, in_elem, 0) {}

template<sp_d T> unique_ptr<MetaMat<T>> SparseSymmMatMUMPS<T>::make_copy() { return std::make_unique<SparseSymmMatMUMPS<T>>(*this); }

#endif

//! @}
