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

#ifndef GMRES_HPP
#define GMRES_HPP

#include <Toolbox/utility.h>
#include "Preconditioner.hpp"

template<sp_d data_t, CanEvaluate<data_t> System, IsPreconditioner<data_t> Preconditioner> int GMRES(const System& system, Col<data_t>& x, const Col<data_t>& b, const Preconditioner& conditioner, const int m, int& max_iteration, data_t& tolerance) {
    constexpr auto ZERO = data_t(0);
    constexpr auto ONE = data_t(1);

    auto generate_rotation = [](const data_t dx, const data_t dy, data_t& cs, data_t& sn) -> void {
        if(suanpan::approx_equal(dy, ZERO)) {
            cs = ONE;
            sn = ZERO;
        }
        else if(std::fabs(dy) > std::fabs(dx)) {
            const data_t fraction = dx / dy;
            sn = ONE / std::sqrt(ONE + fraction * fraction);
            cs = fraction * sn;
        }
        else {
            const data_t fraction = dy / dx;
            cs = ONE / std::sqrt(ONE + fraction * fraction);
            sn = fraction * cs;
        }
    };

    auto apply_rotation = [](data_t& dx, data_t& dy, const data_t cs, const data_t sn) -> void {
        const data_t factor = cs * dx + sn * dy;
        dy = cs * dy - sn * dx;
        dx = factor;
    };

    if(x.empty()) x = conditioner.apply(b);

    const auto mp = m + 1;

    Mat<data_t> hessenberg(mp, m, fill::zeros);

    auto counter = 1;
    data_t beta, residual;
    Col<data_t> s(mp, fill::none), cs(mp, fill::none), sn(mp, fill::none), r;

    auto norm_b = arma::norm(conditioner.apply(b));
    if(suanpan::approx_equal(norm_b, ZERO)) norm_b = ONE;

    auto stop_criterion = [&] {
        residual = (beta = arma::norm(r = conditioner.apply(b - system.evaluate(x)))) / norm_b;
        suanpan_debug("GMRES solver local residual: %.4e.\n", residual);
        if(residual > tolerance) return SUANPAN_FAIL;
        tolerance = residual;
        max_iteration = counter;
        return SUANPAN_SUCCESS;
    };

    if(SUANPAN_SUCCESS == stop_criterion()) return SUANPAN_SUCCESS;

    Mat<data_t> v(b.n_rows, mp, fill::none);

    auto update = [&](const int k) -> Col<data_t> {
        Col<data_t> y = s.head(k + 1llu);

        for(auto i = k; i >= 0; --i) {
            y(i) /= hessenberg(i, i);
            for(auto j = i - 1; j >= 0; --j) y(j) -= hessenberg(j, i) * y(i);
        }

        return v.head_cols(k + 1llu) * y;
    };

    while(counter <= max_iteration) {
        v.col(0) = r / beta;
        s.zeros();
        s(0) = beta;

        for(auto i = 0, j = 1; i < m && counter <= max_iteration; ++i, ++j, ++counter) {
            auto w = conditioner.apply(system.evaluate(v.col(i)));
            for(auto k = 0; k <= i; ++k) w -= (hessenberg(k, i) = arma::dot(w, v.col(k))) * v.col(k);
            v.col(j) = w / (hessenberg(j, i) = arma::norm(w));

            for(auto k = 0; k < i; ++k) apply_rotation(hessenberg(k, i), hessenberg(k + 1llu, i), cs(k), sn(k));

            generate_rotation(hessenberg(i, i), hessenberg(j, i), cs(i), sn(i));
            apply_rotation(hessenberg(i, i), hessenberg(j, i), cs(i), sn(i));
            apply_rotation(s(i), s(j), cs(i), sn(i));

            residual = std::fabs(s(j)) / norm_b;
            suanpan_debug("GMRES solver local residual: %.4e.\n", residual);
            if(residual < tolerance) {
                x += update(i);
                tolerance = residual;
                max_iteration = counter;
                return SUANPAN_SUCCESS;
            }
        }

        x += update(m - 1);
        if(SUANPAN_SUCCESS == stop_criterion()) return SUANPAN_SUCCESS;
    }

    tolerance = residual;
    return SUANPAN_FAIL;
}

#endif
