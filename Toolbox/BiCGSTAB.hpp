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

#ifndef BICGSTAB_HPP
#define BICGSTAB_HPP

#include <Toolbox/utility.h>
#include "Preconditioner.hpp"

template<class Matrix, IsPreconditioner Preconditioner, sp_d data_t> int BiCGSTAB(const Matrix& A, Col<data_t>& x, const Col<data_t>& b, const Preconditioner& conditioner, int& max_iteration, data_t& tolerance) {
    constexpr auto ZERO = data_t(0);
    constexpr auto ONE = data_t(1);

    data_t norm_b = arma::norm(b);
    if(suanpan::approx_equal(norm_b, ZERO)) norm_b = ONE;

    if(x.empty()) x = conditioner.apply(b);

    Col<data_t> r = b - A * x;
    const auto initial_r = r;

    data_t residual = arma::norm(r) / norm_b;
    suanpan_debug("BiCGSTAB solver local residual: %.4e.\n", residual);
    if(residual < tolerance) {
        tolerance = residual;
        max_iteration = 0;
        return 0;
    }

    auto pre_rho = ZERO, alpha = ZERO, omega = ZERO;
    Col<data_t> v, p;

    for(auto i = 1; i <= max_iteration; ++i) {
        const auto rho = arma::dot(initial_r, r);
        if(suanpan::approx_equal(rho, ZERO)) {
            tolerance = residual;
            max_iteration = i;
            return SUANPAN_FAIL;
        }

        if(1 == i) p = r;
        else p = r + rho / pre_rho * alpha / omega * (p - omega * v);

        const auto phat = conditioner.apply(p);
        v = A * phat;
        alpha = rho / arma::dot(initial_r, v);
        const Col<data_t> s = r - alpha * v;

        suanpan_debug("BiCGSTAB solver local residual: %.4e.\n", residual = arma::norm(s) / norm_b);
        if(residual < tolerance) {
            x += alpha * phat;
            tolerance = residual;
            max_iteration = i;
            return SUANPAN_SUCCESS;
        }

        const auto shat = conditioner.apply(s);
        const Col<data_t> t = A * shat;
        omega = arma::dot(t, s) / arma::dot(t, t);
        x += alpha * phat + omega * shat;
        r = s - omega * t;

        pre_rho = rho;

        suanpan_debug("BiCGSTAB solver local residual: %.4e.\n", residual = arma::norm(r) / norm_b);
        if(residual < tolerance) {
            tolerance = residual;
            max_iteration = i;
            return SUANPAN_SUCCESS;
        }

        if(suanpan::approx_equal(omega, ZERO)) {
            tolerance = residual;
            max_iteration = i;
            return SUANPAN_FAIL;
        }
    }

    tolerance = residual;
    return SUANPAN_FAIL;
}

#endif
