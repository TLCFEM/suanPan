
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

#ifndef FGMRES_HPP
#define FGMRES_HPP

#include <mkl_rci.h>
#include <suanPan.h>

template<typename T> requires requires(const T& a, const mat& b) { { a * b } -> std::same_as<mat>; a.n_rows; } int fgmres_solve(const T& kernel, const vec& diagonal, double* left, double* right, const double tolerance) {
    const auto N = static_cast<MKL_INT>(kernel.n_rows);
    // ReSharper disable once CppRedundantCastExpression
    const auto R = std::min(static_cast<MKL_INT>(150), N);

    std::vector work((2 * R + 1) * N + R * (R + 9) / 2 + 1, 0.);

    MKL_INT ipar[128]{};
    double dpar[128]{};

    MKL_INT info{};

    dfgmres_init(&N, nullptr, nullptr, &info, ipar, dpar, work.data());

    ipar[8] = 1;  // residual stopping test
    ipar[9] = 0;  // no user-defined stopping test
    ipar[10] = 1; // use preconditioner
    ipar[11] = 1; // automatic test
    dpar[0] = tolerance;

    while(true) {
        dfgmres(&N, left, right, &info, ipar, dpar, work.data());
        if(-1 == info || -10 == info || -11 == info || -12 == info) {
            suanpan_error("Error code {} received.\n", info);
            return -1;
        }
        if(0 == info || 4 == info && dpar[6] <= dpar[0]) {
            MKL_INT counter;
            dfgmres_get(&N, left, right, &info, ipar, dpar, work.data(), &counter);
            suanpan_debug("Converged in {} iterations.\n", counter);
            return 0;
        }
        const vec xn(&work[ipar[21] - 1], N);
        // ReSharper disable once CppInitializedValueIsAlwaysRewritten
        // ReSharper disable once CppEntityAssignedButNoRead
        vec yn(&work[ipar[22] - 1], N, false, true);
        // ReSharper disable once CppDFAUnusedValue
        if(1 == info) yn = kernel * xn;
        // ReSharper disable once CppDFAUnusedValue
        else if(3 == info) yn = xn / diagonal;
    }
}

#endif
