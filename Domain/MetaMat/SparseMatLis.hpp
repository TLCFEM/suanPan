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
 * @class SparseMatLis
 * @brief A SparseMatLis class that holds matrices.
 *
 * @author tlc
 * @date 15/06/2023
 * @version 0.1.0
 * @file SparseMatLis.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef SPARSEMATLIS_HPP
#define SPARSEMATLIS_HPP

#include "SparseMat.hpp"
#include "csr_form.hpp"

#include <lis/lislib.h>

template<sp_d T> class SparseMatLis final : public SparseMat<T> {
    class lis_vector final {
        LIS_VECTOR v{};

        bool is_set = false;

        auto unset() {
            if(is_set) lis_vector_unset(v);
            is_set = false;
        }

    public:
        explicit lis_vector(const uword n) {
            lis_vector_create(0, &v);
            lis_vector_set_size(v, static_cast<LIS_INT>(n), 0);
        }

        lis_vector(const lis_vector&) = delete;
        lis_vector(lis_vector&&) noexcept = delete;
        lis_vector& operator=(const lis_vector&) = delete;
        lis_vector& operator=(lis_vector&&) noexcept = delete;

        ~lis_vector() {
            unset();
            lis_vector_destroy(v);
        }

        auto set(LIS_SCALAR* value) {
            unset();
            is_set = true;
            lis_vector_set(v, value);
            return v;
        }
    };

    class lis_matrix final {
        LIS_MATRIX a_mat{};

        bool is_set = false;

        auto unset() {
            if(is_set) lis_matrix_unset(a_mat);
            lis_matrix_destroy(a_mat);
            is_set = false;
        }

    public:
        explicit lis_matrix(csr_form<LIS_SCALAR, LIS_INT>& A) { set(A); }

        lis_matrix(const lis_matrix&) = delete;
        lis_matrix(lis_matrix&&) noexcept = delete;
        lis_matrix& operator=(const lis_matrix&) = delete;
        lis_matrix& operator=(lis_matrix&&) noexcept = delete;

        ~lis_matrix() { unset(); }

        auto get() const { return a_mat; }

        auto set(csr_form<LIS_SCALAR, LIS_INT>& A) {
            unset();
            lis_matrix_create(0, &a_mat);
            lis_matrix_set_size(a_mat, A.n_rows, 0);
            lis_matrix_set_csr(A.n_elem, A.row_mem(), A.col_mem(), A.val_mem(), a_mat);
            lis_matrix_assemble(a_mat);
            is_set = true;
        }
    };

    class lis_solver final {
        LIS_SOLVER solver{};

    public:
        lis_solver() { lis_solver_create(&solver); }

        lis_solver(const lis_solver&) { lis_solver_create(&solver); }
        lis_solver(lis_solver&&) noexcept = delete;
        lis_solver& operator=(const lis_solver&) = delete;
        lis_solver& operator=(lis_solver&&) noexcept = delete;

        ~lis_solver() { lis_solver_destroy(solver); }

        auto set_option(const char* option) const { return lis_solver_set_option(option, solver); }

        auto solve(const LIS_MATRIX A, const LIS_VECTOR B, const LIS_VECTOR X) const { return lis_solve(A, B, X, solver); }
    };

    lis_solver solver;

protected:
    using SparseMat<T>::direct_solve;
    using SparseMat<T>::setting;

    int direct_solve(Mat<T>&, const Mat<T>&) override;

public:
    using SparseMat<T>::SparseMat;

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SparseMatLis>(*this); }
};

template<sp_d T> int SparseMatLis<T>::direct_solve(Mat<T>& X, const Mat<T>& B) {
    X.set_size(B.n_rows, B.n_cols);

    csr_form<double, LIS_INT> csr_mat(this->triplet_mat, SparseBase::ZERO, true);

    lis_matrix a(csr_mat);
    lis_vector b(B.n_rows), x(B.n_rows);

    solver.set_option(setting.lis_options.c_str());

    LIS_INT info = 0;
    for(uword I = 0; I < B.n_cols; ++I) {
        info = solver.solve(a.get(), b.set((double*)B.colptr(I)), x.set((double*)X.colptr(I)));
        if(0 != info) break;
    }

    return info;
}

#endif

//! @}
