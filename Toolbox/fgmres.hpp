
/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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
 * @file fgmres.hpp
 * @brief Defensive C++ wrapper around the oneMKL reverse-communication FGMRES
 *        routines (dfgmres_init, dfgmres_check, dfgmres, dfgmres_get).
 *
 * oneMKL keeps the solver state in user-owned ipar/dpar/tmp arrays and asks the
 * caller, via RCI_request, to perform the matrix-vector product (request 1) and
 * to apply the inverse preconditioner (request 3).  This wrapper owns those
 * arrays, validates the setup, services both callbacks, retrieves the solution
 * and maps every oneMKL status to a documented return code.  The reverse
 * communication loop is bounded so a bad request stream cannot spin forever.
 *
 * Reference: Developer Reference for Intel oneAPI Math Kernel Library for C,
 * 2024.2, "dfgmres" and "FGMRES Interface Description".
 *
 * @author tlc
 * @date 16/05/2025
 * @version 0.2.0
 * @addtogroup Utility
 * @{
 */

#ifndef FGMRES_HPP
#define FGMRES_HPP

#include <mkl_rci.h>
#include <suanPan.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <vector>

/**
 * Solve A*x = b with the oneMKL FGMRES solver, using diagonal preconditioning.
 *
 * @param kernel    matrix-like object with n_rows and kernel * mat -> mat
 * @param diagonal  diagonal of A, used as the preconditioner (length n_rows)
 * @param left      solution vector x; on return contains the approximated solution
 * @param right     right-hand side b; kept unchanged (ipar[12] = 0)
 * @param tolerance relative residual tolerance (dpar[0])
 * @param itercount optional output for the number of performed iterations
 *
 * @return 0 on convergence; -1 invalid input; -2 init failed; -3 check failed;
 *         -4 numerical error (oneMKL error or non-finite solution);
 *         -5 request limit exceeded; -6 get failed;
 *         -7 iteration limit reached without convergence.
 */
template<typename T> requires requires(const T& a, const mat& b) { { a* b } -> std::same_as<mat>; a.n_rows; }
int fgmres_solve(const T& kernel, const vec& diagonal, double* left, double* right, const double tolerance, MKL_INT* itercount = nullptr) {
    const auto n = static_cast<MKL_INT>(kernel.n_rows);
    const auto r = std::min(MKL_INT{150}, n);

    if(n <= 0 || nullptr == left || nullptr == right || !std::isfinite(tolerance) || tolerance <= 0. || diagonal.n_elem != static_cast<uword>(n)) {
        suanpan_error("FGMRES: invalid arguments.\n");
        return -1;
    }

    for(const auto value : diagonal)
        if(0. == value || !std::isfinite(value)) {
            suanpan_error("FGMRES: preconditioner diagonal must be finite and nonzero.\n");
            return -1;
        }

    // oneMKL requires ((2*r + 1)*n + r*(r + 9)/2 + 1) double slots in tmp.
    std::vector<double> work(static_cast<std::size_t>((2 * r + 1) * n + r * (r + 9) / 2 + 1), 0.);
    std::array<MKL_INT, 128> ipar{};
    std::array<double, 128> dpar{};
    MKL_INT info{};

    // Pass x and b as well: a non-null initial guess makes the call well-defined.
    dfgmres_init(&n, left, right, &info, ipar.data(), dpar.data(), work.data());
    if(0 != info) {
        suanpan_error("FGMRES: initialization failed with code {}.\n", info);
        return -2;
    }

    ipar[8] = 1;  // residual stopping test
    ipar[9] = 0;  // no user-defined stopping test
    ipar[10] = 1; // use preconditioner
    ipar[11] = 1; // automatic zero-norm test
    ipar[12] = 0; // write the solution into x and keep b untouched
    dpar[0] = tolerance;

    // Recommended after changing ipar/dpar; it also fixes ipar[15:20] if needed.
    dfgmres_check(&n, left, right, &info, ipar.data(), dpar.data(), work.data());
    if(0 != info) {
        // -1001/-1010/-1011 only report warnings or adjusted parameters.
        if(info < -1011) {
            suanpan_error("FGMRES: parameter check failed with code {}.\n", info);
            return -3;
        }
        suanpan_debug("FGMRES: parameter check reported code {}.\n", info);
    }

    const auto limit = std::max(MKL_INT{32}, 10 * ipar[4] + 32);
    bool converged = false;
    bool exhausted = false;

    for(MKL_INT count = 0; count < limit && !converged && !exhausted; ++count) {
        dfgmres(&n, left, right, &info, ipar.data(), dpar.data(), work.data());

        switch(info) {
        case 0:
            converged = true;
            break;
        case 1: // y = A*x
        case 3: { // y = M^-1*x
            if(ipar[21] < 1 || ipar[22] < 1 || static_cast<std::size_t>(ipar[21] + n - 1) > work.size() || static_cast<std::size_t>(ipar[22] + n - 1) > work.size()) {
                suanpan_error("FGMRES: invalid work offsets returned by oneMKL.\n");
                return -4;
            }
            const vec xn(&work[ipar[21] - 1], static_cast<uword>(n), false, true);
            vec yn(&work[ipar[22] - 1], static_cast<uword>(n), false, true);
            if(1 == info) {
                const mat product = kernel * xn;
                if(product.n_rows != static_cast<uword>(n) || product.n_cols != 1) {
                    suanpan_error("FGMRES: matrix-vector product has unexpected size.\n");
                    return -4;
                }
                yn = product;
            }
            else yn = xn / diagonal;
            break;
        }
        case 2: // user stopping test, normally disabled
            if(dpar[4] <= dpar[3]) converged = true;
            break;
        case 4: // zero-norm check, compare with dpar[7], not the residual tolerance
            if(dpar[6] <= dpar[7]) converged = true;
            break;
        case -1: // iteration limit reached before convergence
            exhausted = true;
            break;
        default:
            suanpan_error("FGMRES: oneMKL returned error code {}.\n", info);
            return -4;
        }
    }

    if(!converged && !exhausted) {
        suanpan_error("FGMRES: exceeded {} reverse-communication requests.\n", limit);
        return -5;
    }

    MKL_INT counter{};
    dfgmres_get(&n, left, right, &info, ipar.data(), dpar.data(), work.data(), &counter);
    if(nullptr != itercount) *itercount = counter;
    if(0 != info) {
        suanpan_error("FGMRES: could not retrieve the solution (code {}).\n", info);
        return -6;
    }

    for(uword i = 0; i < static_cast<uword>(n); ++i)
        if(!std::isfinite(left[i])) {
            suanpan_error("FGMRES: solution contains non-finite values.\n");
            return -4;
        }

    if(!converged) {
        suanpan_debug("FGMRES: stopped after {} iterations without convergence.\n", counter);
        return -7;
    }

    suanpan_debug("FGMRES: converged in {} iterations.\n", counter);
    return 0;
}

#endif

//! @}
