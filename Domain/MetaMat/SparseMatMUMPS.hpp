/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
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

// ReSharper disable CppClangTidyClangDiagnosticMissingFieldInitializers
// ReSharper disable CppCStyleCast
#ifndef SPARSEMATMUMPS_HPP
#define SPARSEMATMUMPS_HPP

#include "SparseMat.hpp"

#include <mumps/dmumps_c.h>
#include <mumps/smumps_c.h>

namespace mumps {
    template<typename> struct struc {};
    template<> struct struc<double> {
        using struct_type = DMUMPS_STRUC_C;
        using entry_type = double;
        static auto mumps_c(DMUMPS_STRUC_C* ptr) { return dmumps_c(ptr); }
    };
    template<> struct struc<float> {
        using struct_type = SMUMPS_STRUC_C;
        using entry_type = float;
        static auto mumps_c(SMUMPS_STRUC_C* ptr) { return smumps_c(ptr); }
    };
} // namespace mumps

template<sp_d T> class SparseMatBaseMUMPS : public SparseMat<T> {
    using struct_t = typename mumps::struc<T>::struct_type;
    using entry_t = typename mumps::struc<T>::entry_type;

    const la_it sym;

    struct_t id{sym, 1, -1, -987654};

    triplet_form<T, la_it> coo_mat;

    auto perform_job(const int job) {
        id.job = job;
        mumps::struc<T>::mumps_c(&id);
    }

    auto init_config() {
        id.icntl[3] = 0;   // level of printing for error, warning, and diagnostic messages
        id.icntl[9] = 2;   // iterative refinement to the computed solution
        id.icntl[13] = 50; // percentage increase in the estimated working space
        id.icntl[19] = 0;  // dense rhs
        id.icntl[32] = 1;  // determinant
        id.icntl[34] = 1;  // BLR
    }

protected:
    int direct_solve(Mat<T>& X, Mat<T>&& B) override {
        if(!this->factored) {
            if(0 == sym) coo_mat = triplet_form<T, la_it>(this->triplet_mat, SparseBase::ONE, false);
            else {
                auto lower_mat = this->triplet_mat.lower();
                coo_mat = triplet_form<T, la_it>(lower_mat, SparseBase::ONE, false);
            }

            id.n = coo_mat.n_rows;
            id.nnz = coo_mat.n_elem;
            id.irn = coo_mat.row_mem();
            id.jcn = coo_mat.col_mem();
            id.a = (entry_t*)coo_mat.val_mem();

            perform_job(4);

            this->factored = true;
        }

        id.lrhs = B.n_rows;
        id.nrhs = B.n_cols;
        id.rhs = (entry_t*)B.memptr();

        perform_job(3);

        if(id.infog[0] < 0) return SUANPAN_FAIL;

        X = std::move(B);

        return SUANPAN_SUCCESS;
    }

    int direct_solve(Mat<T>& X, const Mat<T>& B) override { return this->direct_solve(X, Mat<T>(B)); }

public:
    SparseMatBaseMUMPS(const uword in_row, const uword in_col, const uword in_elem, const int in_sym)
        : SparseMat<T>(in_row, in_col, in_elem)
        , sym(in_sym) {
        perform_job(-1);
        init_config();
    }

    SparseMatBaseMUMPS(const SparseMatBaseMUMPS& other)
        : SparseMat<T>(other)
        , sym(other.sym) {
        perform_job(-1);
        init_config();
        this->factored = false;
    }

    SparseMatBaseMUMPS(SparseMatBaseMUMPS&&) noexcept = delete;
    SparseMatBaseMUMPS& operator=(const SparseMatBaseMUMPS&) = delete;
    SparseMatBaseMUMPS& operator=(SparseMatBaseMUMPS&&) noexcept = delete;

    ~SparseMatBaseMUMPS() override { perform_job(-2); }

    [[nodiscard]] int sign_det() const override { return id.rinfog[11] < entry_t{0} ? -1 : 1; }
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
        : SparseMatBaseMUMPS<T>(in_row, in_col, in_elem, 2) {}

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SparseSymmMatMUMPS>(*this); }
};

#endif

//! @}
