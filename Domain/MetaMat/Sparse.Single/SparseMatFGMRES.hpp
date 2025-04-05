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
#ifndef SPARSEMATFGMRES_HPP
#define SPARSEMATFGMRES_HPP

#include "../SparseMat.hpp"
#include "../csr_form.hpp"

#include <mkl_rci.h>

template<sp_d T> class SparseMatFGMRES final : public SparseMat<T> {
    MKL_INT ipar[128]{};
    double dpar[128]{};

    podarray<double> work;

protected:
    using SparseMat<T>::direct_solve;

    int direct_solve(Mat<T>&, const Mat<T>&) override;

public:
    SparseMatFGMRES(const uword in_row, const uword in_col, const uword in_elem = 0)
        : SparseMat<T>(in_row, in_col, in_elem) {}

    SparseMatFGMRES(const SparseMatFGMRES& other)
        : SparseMat<T>(other) {}

    SparseMatFGMRES(SparseMatFGMRES&&) noexcept = delete;
    SparseMatFGMRES& operator=(const SparseMatFGMRES&) = delete;
    SparseMatFGMRES& operator=(SparseMatFGMRES&&) noexcept = delete;

    ~SparseMatFGMRES() override = default;

    unique_ptr<MetaMat<T>> make_copy() override { return std::make_unique<SparseMatFGMRES>(*this); }
};

template<sp_d T> int SparseMatFGMRES<T>::direct_solve(Mat<T>& X, const Mat<T>& B) {
    const auto N = static_cast<MKL_INT>(B.n_rows);
    // ReSharper disable once CppRedundantCastExpression
    const auto R = std::min(static_cast<MKL_INT>(150), N);

    work.zeros((2 * R + 1) * N + R * (R + 9) / 2 + 1);

    X = B;

    const auto precond = this->triplet_mat.diag();

    csr_form<T, int> csr_mat(this->triplet_mat);

    MKL_INT info;

    dfgmres_init(&N, nullptr, nullptr, &info, ipar, dpar, work.memptr());
    if(0 != info) return info;

    ipar[8] = 1;
    ipar[9] = 0;
    ipar[10] = 1; // use preconditioner
    ipar[11] = 1;
    dpar[0] = this->setting.tolerance;

    for(auto I = 0llu; I < B.n_cols; ++I) {
        while(true) {
            dfgmres(&N, (double*)X.colptr(I), (double*)B.colptr(I), &info, ipar, dpar, work.memptr()); // NOLINT(clang-diagnostic-cast-qual)
            if(-1 == info || -10 == info || -11 == info || -12 == info) {
                suanpan_error("Error code {} received.\n", info);
                return -1;
            }
            if(0 == info || 4 == info && dpar[6] <= dpar[0]) {
                MKL_INT counter;
                dfgmres_get(&N, (double*)X.colptr(I), (double*)B.colptr(I), &info, ipar, dpar, work.memptr(), &counter); // NOLINT(clang-diagnostic-cast-qual)
                suanpan_debug("Converged in {} iterations.\n", counter);
                break;
            }
            if(1 == info) {
                const vec xn(&work[ipar[21] - 1], N);
                // ReSharper disable once CppInitializedValueIsAlwaysRewritten
                // ReSharper disable once CppEntityAssignedButNoRead
                vec yn(&work[ipar[22] - 1], N, false, true);
                // ReSharper disable once CppDFAUnusedValue
                yn = csr_mat * xn;
            }
            else if(3 == info) {
                const vec xn(&work[ipar[21] - 1], N);
                // ReSharper disable once CppInitializedValueIsAlwaysRewritten
                // ReSharper disable once CppEntityAssignedButNoRead
                vec yn(&work[ipar[22] - 1], N, false, true);
                // ReSharper disable once CppDFAUnusedValue
                yn = xn / precond;
            }
        }
    }

    return info;
}

#endif

//! @}
