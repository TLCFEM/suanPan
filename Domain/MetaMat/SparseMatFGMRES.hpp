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
 * @class SparseMatFGMRES
 * @brief A SparseMatFGMRES class that holds matrices.
 *
 * @author tlc
 * @date 27/03/2025
 * @version 0.2.0
 * @file SparseMatFGMRES.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef SPARSEMATBASEFGMRES_HPP
#define SPARSEMATBASEFGMRES_HPP

#ifdef SUANPAN_MKL

#include "SparseMat.hpp"
#include "csr_form.hpp"

#include <mkl_rci.h>

template<sp_d T> class SparseMatBaseFGMRES : public SparseMat<T> {
    MKL_INT ipar[128]{};
    double dpar[128]{};

    podarray<double> work;

protected:
    using SparseMat<T>::direct_solve;

    int direct_solve(Mat<T>&, const Mat<T>&) override;

public:
    SparseMatBaseFGMRES(const uword in_row, const uword in_col, const uword in_elem)
        : SparseMat<T>(in_row, in_col, in_elem) {}

    SparseMatBaseFGMRES(const SparseMatBaseFGMRES&) = default;
    SparseMatBaseFGMRES(SparseMatBaseFGMRES&&) noexcept = delete;
    SparseMatBaseFGMRES& operator=(const SparseMatBaseFGMRES&) = delete;
    SparseMatBaseFGMRES& operator=(SparseMatBaseFGMRES&&) noexcept = delete;

    ~SparseMatBaseFGMRES() override = default;
};

template<sp_d T> int SparseMatBaseFGMRES<T>::direct_solve(Mat<T>& X, const Mat<T>& B) {
    const auto N = static_cast<MKL_INT>(B.n_rows);

    const auto restart = std::min(150, N);

    work.zeros((2 * restart + 1) * N + restart * (restart + 9) / 2 + 1);

    X = B;

    const auto precond = this->triplet_mat.diag();

    csr_form<T, int> csr_mat(this->triplet_mat);

    MKL_INT request;

    dfgmres_init(&N, nullptr, nullptr, &request, ipar, dpar, work.memptr());
    if(0 != request) return request;

    ipar[8] = 1;
    ipar[9] = 0;
    ipar[10] = 1; // use preconditioner
    ipar[11] = 1;
    dpar[0] = this->setting.tolerance;

    for(auto I = 0llu; I < B.n_cols; ++I) {
        while(true) {
            dfgmres(&N, (double*)X.colptr(I), (double*)B.colptr(I), &request, ipar, dpar, work.memptr());
            if(-1 == request || -10 == request || -11 == request || -12 == request) {
                suanpan_error("Error code {} received.\n", request);
                return -1;
            }
            if(0 == request || 4 == request && dpar[6] <= dpar[0]) {
                MKL_INT counter;
                dfgmres_get(&N, (double*)X.colptr(I), (double*)B.colptr(I), &request, ipar, dpar, work.memptr(), &counter);
                suanpan_debug("Converged in {} iterations.\n", counter);
                break;
            }
            if(1 == request) {
                const vec xn(&work[ipar[21] - 1llu], N);
                // ReSharper disable once CppInitializedValueIsAlwaysRewritten
                // ReSharper disable once CppEntityAssignedButNoRead
                vec yn(&work[ipar[22] - 1llu], N, false, true);
                // ReSharper disable once CppDFAUnusedValue
                yn = csr_mat * xn;
            }
            else if(3 == request) {
                const vec xn(&work[ipar[21] - 1llu], N);
                // ReSharper disable once CppInitializedValueIsAlwaysRewritten
                // ReSharper disable once CppEntityAssignedButNoRead
                vec yn(&work[ipar[22] - 1llu], N, false, true);
                // ReSharper disable once CppDFAUnusedValue
                yn = xn / precond;
            }
        }
    }

    return request;
}

template<sp_d T> class SparseMatFGMRES final : public SparseMatBaseFGMRES<T> {
public:
    SparseMatFGMRES(const uword in_row, const uword in_col, const uword in_elem = 0)
        : SparseMatBaseFGMRES<T>(in_row, in_col, in_elem) {}

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SparseMatFGMRES>(*this); }
};

template<sp_d T> using SparseSymmMatFGMRES = SparseMatFGMRES<T>;

#endif

#endif

//! @}
