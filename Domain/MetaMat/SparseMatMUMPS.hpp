/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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
#include <mumps/smumps_c.h>
#include "SparseMat.hpp"

template<sp_d T> class SparseMatBaseMUMPS : public SparseMat<T> {
    const int sym = 0;

    DMUMPS_STRUC_C dmumps_job{sym, 1, -1, -987654};
    SMUMPS_STRUC_C smumps_job{sym, 1, -1, -987654};

    triplet_form<float, int> s_mat;

    s32_vec l_irn, l_jrn;

    template<bool convert, typename ST,std::invocable<ST*> F, typename COO> int alloc(ST& mumps_job, F& mumps_c, COO& triplet_mat) {
        if(this->factored) return 0;

        dealloc(mumps_job, mumps_c);

        this->factored = true;

        mumps_job.job = -1;
        mumps_c(&mumps_job);

        triplet_mat.csc_condense();

        if constexpr(convert) {
            l_irn.set_size(triplet_mat.n_elem);
            l_jrn.set_size(triplet_mat.n_elem);

            suanpan_for(0, static_cast<int>(triplet_mat.n_elem), [&](const int I) {
                l_irn[I] = static_cast<int>(triplet_mat.row(I) + 1);
                l_jrn[I] = static_cast<int>(triplet_mat.col(I) + 1);
            });

            mumps_job.irn = l_irn.memptr();
            mumps_job.jcn = l_jrn.memptr();
        }
        else {
            mumps_job.irn = triplet_mat.row_mem();
            mumps_job.jcn = triplet_mat.col_mem();
        }

        mumps_job.a = triplet_mat.val_mem();

        mumps_job.n = static_cast<int>(triplet_mat.n_rows);
        mumps_job.nnz = static_cast<int64_t>(triplet_mat.n_elem);

        mumps_job.icntl[0] = -1;
        mumps_job.icntl[1] = -1;
        mumps_job.icntl[2] = -1;
        mumps_job.icntl[3] = 0;
        mumps_job.icntl[13] = 100;
        mumps_job.icntl[19] = 0; // dense rhs
        mumps_job.icntl[32] = 1; // determinant
        mumps_job.icntl[34] = 1; // BLR

        mumps_job.job = 4;
        mumps_c(&mumps_job);

        if(0 != mumps_job.info[0])
            suanpan_error("Error code {} received.\n", mumps_job.info[0]);

        return mumps_job.info[0];
    }

    template<typename ST,std::invocable<ST*> F> static void dealloc(ST& mumps_job, F& mumps_c) {
        if(3 != mumps_job.job) return;
        mumps_job.job = -2;
        mumps_c(&mumps_job);
    }

    template<typename ST,std::invocable<ST*> F> static void run(ST& mumps_job, F& mumps_c) {
        mumps_job.job = 3;
        mumps_c(&mumps_job);
    }

protected:
    int direct_solve(Mat<T>& X, Mat<T>&& B) override {
        int INFO;

        if constexpr(std::is_same_v<T, float>) {
            if(0 != (INFO = alloc<true>(smumps_job, smumps_c, this->triplet_mat))) return INFO;

            smumps_job.rhs = B.memptr();
            smumps_job.lrhs = static_cast<int>(B.n_rows);
            smumps_job.nrhs = static_cast<int>(B.n_cols);

            run(smumps_job, smumps_c);

            X = std::move(B);

            INFO = smumps_job.info[0];
        }
        else if(Precision::FULL == this->setting.precision) {
            if(0 != (INFO = alloc<true>(dmumps_job, dmumps_c, this->triplet_mat))) return INFO;

            dmumps_job.rhs = B.memptr();
            dmumps_job.lrhs = static_cast<int>(B.n_rows);
            dmumps_job.nrhs = static_cast<int>(B.n_cols);

            run(dmumps_job, dmumps_c);

            X = std::move(B);

            INFO = dmumps_job.info[0];
        }
        else {
            if(!this->factored) s_mat = triplet_form<float, int>(this->triplet_mat, SparseBase::ONE, false);

            if(0 != (INFO = alloc<false>(smumps_job, smumps_c, s_mat))) return INFO;

            INFO = this->mixed_trs(X, std::forward<Mat<T>>(B), [&](fmat& residual) {
                smumps_job.rhs = residual.memptr();
                smumps_job.lrhs = static_cast<int>(residual.n_rows);
                smumps_job.nrhs = static_cast<int>(residual.n_cols);

                run(smumps_job, smumps_c);

                return smumps_job.info[0];
            });
        }

        return INFO;
    }

    int direct_solve(Mat<T>& X, const Mat<T>& B) override { return this->direct_solve(X, Mat<T>(B)); }

public:
    SparseMatBaseMUMPS(const uword in_row, const uword in_col, const uword in_elem, const int in_sym)
        : SparseMat<T>(in_row, in_col, in_elem)
        , sym(in_sym) {}

    SparseMatBaseMUMPS(const SparseMatBaseMUMPS& other)
        : SparseMat<T>(other)
        , sym(other.sym)
        , dmumps_job{other.sym, 1, -1, -987654}
        , smumps_job{other.sym, 1, -1, -987654}
        , l_irn(other.l_irn)
        , l_jrn(other.l_jrn) {}

    SparseMatBaseMUMPS(SparseMatBaseMUMPS&&) noexcept = delete;
    SparseMatBaseMUMPS& operator=(const SparseMatBaseMUMPS&) = delete;
    SparseMatBaseMUMPS& operator=(SparseMatBaseMUMPS&&) noexcept = delete;

    ~SparseMatBaseMUMPS() override {
        dealloc(dmumps_job, dmumps_c);
        dealloc(smumps_job, smumps_c);
    }

    void zeros() override {
        SparseMat<T>::zeros();
        dealloc(dmumps_job, dmumps_c);
        dealloc(smumps_job, smumps_c);
    }

    [[nodiscard]] int sign_det() const override {
        if(IterativeSolver::NONE != this->setting.iterative_solver) throw invalid_argument("analysis requires the sign of determinant but iterative solver does not support it");

        int det_sign;

        if constexpr(std::is_same_v<T, float>) det_sign = smumps_job.rinfog[11] < 0. ? -1 : 1;
        else if(Precision::FULL == this->setting.precision) det_sign = dmumps_job.rinfog[11] < 0. ? -1 : 1;
        else det_sign = smumps_job.rinfog[11] < 0. ? -1 : 1;

        return det_sign;
    }
};

template<sp_d T> class SparseMatMUMPS final : public SparseMatBaseMUMPS<T> {
public:
    SparseMatMUMPS(const uword in_row, const uword in_col, const uword in_elem = 0)
        : SparseMatBaseMUMPS<T>(in_row, in_col, in_elem, 0) {}

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SparseMatMUMPS>(*this); }
};

template<sp_d T> class SparseSymmMatMUMPS final : public SparseMatBaseMUMPS<T> {
public:
    SparseSymmMatMUMPS(const uword in_row, const uword in_col, const uword in_elem = 0)
        : SparseMatBaseMUMPS<T>(in_row, in_col, in_elem, 0) {}

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SparseSymmMatMUMPS>(*this); }
};

#endif

//! @}
