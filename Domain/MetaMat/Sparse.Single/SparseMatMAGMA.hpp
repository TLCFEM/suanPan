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
 * @class SparseMatMAGMA
 * @brief A SparseMatMAGMA class that holds matrices.
 *
 * TODO: improve performance by storing factorization and reusing it
 *
 * @author tlc
 * @date 24/02/2023
 * @version 0.1.0
 * @file SparseMatMAGMA.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
// ReSharper disable CppClangTidyBugproneBranchClone
// ReSharper disable CppClangTidyClangDiagnosticMissingFieldInitializers
#ifndef SPARSEMATMAGMA_HPP
#define SPARSEMATMAGMA_HPP

#include "SparseMat.hpp"

#include <variant>

using magma_opt_t = std::variant<magma_dopts, magma_sopts>;

template<sp_d T> class SparseMatMAGMA final : public SparseMat<T> {
    magma_dopts dopts;
    magma_sopts sopts;

    magma_queue_t queue{};

    csr_form<T, magma_index_t> csr_mat;

    magma_s_matrix A_f{Magma_CSR, Magma_CPU};
    magma_s_matrix dA_f{Magma_CSR, Magma_DEV};
    magma_s_matrix b_f{Magma_DENSE, Magma_CPU};
    magma_s_matrix db_f{Magma_DENSE, Magma_DEV};
    magma_s_matrix x_f{Magma_DENSE, Magma_CPU};
    magma_s_matrix dx_f{Magma_DENSE, Magma_DEV};

    magma_d_matrix A_d{Magma_CSR, Magma_CPU};
    magma_d_matrix dA_d{Magma_CSR, Magma_DEV};
    magma_d_matrix b_d{Magma_DENSE, Magma_CPU};
    magma_d_matrix db_d{Magma_DENSE, Magma_DEV};
    magma_d_matrix x_d{Magma_DENSE, Magma_CPU};
    magma_d_matrix dx_d{Magma_DENSE, Magma_DEV};

protected:
    int direct_solve(Mat<T>& out_mat, const Mat<T>& in_mat) override { return this->direct_solve(out_mat, Mat<T>(in_mat)); }

    int direct_solve(Mat<T>&, Mat<T>&&) override;

public:
    SparseMatMAGMA(const uword in_row, const uword in_col, const magma_opt_t& in_opts, const uword in_elem = 0)
        : SparseMat<T>(in_row, in_col, in_elem) {
        magma_init();
        if(SUANPAN_VERBOSE) magma_print_environment();
        magma_queue_create(0, &queue);
        if constexpr(std::is_same_v<T, float>) {
            sopts = std::get<magma_sopts>(in_opts);
            magma_ssolverinfo_init(&sopts.solver_par, &sopts.precond_par, queue);
        }
        else {
            dopts = std::get<magma_dopts>(in_opts);
            magma_dsolverinfo_init(&dopts.solver_par, &dopts.precond_par, queue);
        }
    }

    SparseMatMAGMA(const SparseMatMAGMA&);
    SparseMatMAGMA(SparseMatMAGMA&&) noexcept = delete;
    SparseMatMAGMA& operator=(const SparseMatMAGMA&) = delete;
    SparseMatMAGMA& operator=(SparseMatMAGMA&&) noexcept = delete;

    ~SparseMatMAGMA() override {
        magma_smfree(&dx_f, queue);
        magma_smfree(&db_f, queue);
        magma_smfree(&dA_f, queue);
        magma_smfree(&x_f, queue);
        magma_smfree(&b_f, queue);
        magma_smfree(&A_f, queue);
        magma_dmfree(&dx_d, queue);
        magma_dmfree(&db_d, queue);
        magma_dmfree(&dA_d, queue);
        magma_dmfree(&x_d, queue);
        magma_dmfree(&b_d, queue);
        magma_dmfree(&A_d, queue);
        magma_queue_destroy(queue);
        magma_finalize();
    }

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SparseMatMAGMA>(*this); }
};

template<sp_d T> SparseMatMAGMA<T>::SparseMatMAGMA(const SparseMatMAGMA& other)
    : SparseMat<T>(other)
    , dopts(other.dopts)
    , sopts(other.sopts)
    , csr_mat(other.csr_mat) {
    magma_queue_create(0, &queue);
    if constexpr(std::is_same_v<T, float>) {
        magma_ssolverinfo_init(&sopts.solver_par, &sopts.precond_par, queue);
        if(this->factored) {
            magma_scsrset(csr_mat.n_rows, csr_mat.n_cols, csr_mat.row_mem(), csr_mat.col_mem(), csr_mat.val_mem(), &A_f, queue);
            magma_smtransfer(A_f, &dA_f, Magma_CPU, Magma_DEV, queue);
        }
    }
    else {
        magma_dsolverinfo_init(&dopts.solver_par, &dopts.precond_par, queue);
        if(this->factored) {
            magma_dcsrset(csr_mat.n_rows, csr_mat.n_cols, csr_mat.row_mem(), csr_mat.col_mem(), csr_mat.val_mem(), &A_d, queue);
            magma_dmtransfer(A_d, &dA_d, Magma_CPU, Magma_DEV, queue);
        }
    }
}

template<sp_d T> int SparseMatMAGMA<T>::direct_solve(Mat<T>& X, Mat<T>&& B) {
    auto b_rows = static_cast<magma_index_t>(B.n_rows), b_cols = static_cast<magma_index_t>(B.n_cols);

    if constexpr(std::is_same_v<T, float>) {
        if(!this->factored) {
            csr_mat = csr_form<float, magma_index_t>(this->triplet_mat, SparseBase::ZERO);
            this->factored = true;
            magma_scsrset(csr_mat.n_rows, csr_mat.n_cols, csr_mat.row_mem(), csr_mat.col_mem(), csr_mat.val_mem(), &A_f, queue);
            magma_smtransfer(A_f, &dA_f, Magma_CPU, Magma_DEV, queue);
        }

        magma_svset(b_rows, b_cols, (float*)B.memptr(), &b_f, queue); // NOLINT(clang-diagnostic-cast-qual)
        magma_smtransfer(b_f, &db_f, Magma_CPU, Magma_DEV, queue);

        magma_s_precondsetup(dA_f, db_f, &sopts.solver_par, &sopts.precond_par, queue);

        magma_svinit(&dx_f, Magma_DEV, b_rows, b_cols, T(0), queue);

        magma_s_solver(dA_f, db_f, &dx_f, &sopts, queue);

        magma_smtransfer(dx_f, &x_f, Magma_DEV, Magma_CPU, queue);

        X = std::forward<Mat<float>>(B);
        magma_svcopy(x_f, &b_rows, &b_cols, X.memptr(), queue);
    }
    else {
        if(!this->factored) {
            csr_mat = csr_form<double, magma_index_t>(this->triplet_mat, SparseBase::ZERO);
            this->factored = true;
            magma_dcsrset(csr_mat.n_rows, csr_mat.n_cols, csr_mat.row_mem(), csr_mat.col_mem(), csr_mat.val_mem(), &A_d, queue);
            magma_dmtransfer(A_d, &dA_d, Magma_CPU, Magma_DEV, queue);
        }

        magma_dvset(b_rows, b_cols, (double*)B.memptr(), &b_d, queue); // NOLINT(clang-diagnostic-cast-qual)
        magma_dmtransfer(b_d, &db_d, Magma_CPU, Magma_DEV, queue);

        magma_d_precondsetup(dA_d, db_d, &dopts.solver_par, &dopts.precond_par, queue);

        magma_dvinit(&dx_d, Magma_DEV, b_rows, b_cols, T(0), queue);

        magma_d_solver(dA_d, db_d, &dx_d, &dopts, queue);

        magma_dmtransfer(dx_d, &x_d, Magma_DEV, Magma_CPU, queue);

        X = std::forward<Mat<double>>(B);
        magma_dvcopy(x_d, &b_rows, &b_cols, X.memptr(), queue);
    }

    return SUANPAN_SUCCESS;
}

#endif

//! @}
