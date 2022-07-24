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

#ifndef ITERATIVESOLVER_HPP
#define ITERATIVESOLVER_HPP

#include <Toolbox/utility.h>
#include "SolverSetting.hpp"

template<typename T, typename data_t> concept HasEvaluate = requires(T* t, const Col<data_t>& x) { { t->evaluate(x) } -> std::convertible_to<Col<data_t>> ; };

template<sp_d data_t, HasEvaluate<data_t> System> int GMRES(const System* system, Col<data_t>& x, const Col<data_t>& b, SolverSetting<data_t>& setting) {
    constexpr sp_d auto ZERO = data_t(0);
    constexpr sp_d auto ONE = data_t(1);

    const auto& conditioner = setting.preconditioner;

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

    if(x.empty()) x = conditioner->apply(b);
    else x.zeros(arma::size(b));

    const auto mp = setting.restart + 1;

    Mat<data_t> hessenberg(mp, setting.restart, fill::zeros);

    auto counter = 1;
    data_t beta, residual;
    Col<data_t> s(mp, fill::none), cs(mp, fill::none), sn(mp, fill::none), r;

    auto norm_b = arma::norm(conditioner->apply(b));
    if(suanpan::approx_equal(norm_b, ZERO)) norm_b = ONE;

    auto stop_criterion = [&] {
        residual = (beta = arma::norm(r = conditioner->apply(b - system->evaluate(x)))) / norm_b;
        suanpan_debug("GMRES solver local residual: %.4e.\n", residual);
        if(residual > setting.tolerance) return SUANPAN_FAIL;
        setting.tolerance = residual;
        setting.max_iteration = counter;
        return SUANPAN_SUCCESS;
    };

    if(SUANPAN_SUCCESS == stop_criterion()) return SUANPAN_SUCCESS;

    Mat<data_t> v(b.n_rows, mp, fill::none);

    auto update = [&](const int k) -> Col<data_t> {
        Col<data_t> y = s.head(k + 1llu);

        for(auto i = k; i >= 0; --i) {
            y(i) /= hessenberg(i, i);
            y.head(i) -= hessenberg.col(i).head(i) * y(i);
        }

        return v.head_cols(k + 1llu) * y;
    };

    while(counter <= setting.max_iteration) {
        v.col(0) = r / beta;
        s.zeros();
        s(0) = beta;

        for(auto i = 0, j = 1; i < setting.restart && counter <= setting.max_iteration; ++i, ++j, ++counter) {
            auto w = conditioner->apply(system->evaluate(v.col(i)));
            for(auto k = 0; k <= i; ++k) w -= (hessenberg(k, i) = arma::dot(w, v.col(k))) * v.col(k);
            v.col(j) = w / (hessenberg(j, i) = arma::norm(w));

            for(auto k = 0; k < i; ++k) apply_rotation(hessenberg(k, i), hessenberg(k + 1llu, i), cs(k), sn(k));

            generate_rotation(hessenberg(i, i), hessenberg(j, i), cs(i), sn(i));
            apply_rotation(hessenberg(i, i), hessenberg(j, i), cs(i), sn(i));
            apply_rotation(s(i), s(j), cs(i), sn(i));

            residual = std::fabs(s(j)) / norm_b;
            suanpan_debug("GMRES solver local residual: %.4e.\n", residual);
            if(residual < setting.tolerance) {
                x += update(i);
                setting.tolerance = residual;
                setting.max_iteration = counter;
                return SUANPAN_SUCCESS;
            }
        }

        x += update(setting.restart - 1);
        if(SUANPAN_SUCCESS == stop_criterion()) return SUANPAN_SUCCESS;
    }

    setting.tolerance = residual;
    return SUANPAN_FAIL;
}

template<sp_d data_t, HasEvaluate<data_t> System> int BiCGSTAB(const System* system, Col<data_t>& x, const Col<data_t>& b, SolverSetting<data_t>& setting) {
    constexpr sp_d auto ZERO = data_t(0);
    constexpr sp_d auto ONE = data_t(1);

    const auto& conditioner = setting.preconditioner;

    data_t norm_b = arma::norm(b);
    if(suanpan::approx_equal(norm_b, ZERO)) norm_b = ONE;

    if(x.empty()) x = conditioner->apply(b);
    else x.zeros(arma::size(b));

    Col<data_t> r = b - system->evaluate(x);
    const auto initial_r = r;

    data_t residual = arma::norm(r) / norm_b;
    suanpan_debug("BiCGSTAB solver local residual: %.4e.\n", residual);
    if(residual < setting.tolerance) {
        setting.tolerance = residual;
        setting.max_iteration = 0;
        return 0;
    }

    sp_d auto pre_rho = ZERO, alpha = ZERO, omega = ZERO;
    Col<data_t> v, p;

    for(auto i = 1; i <= setting.max_iteration; ++i) {
        const auto rho = arma::dot(initial_r, r);
        if(suanpan::approx_equal(rho, ZERO)) {
            setting.tolerance = residual;
            setting.max_iteration = i;
            return SUANPAN_FAIL;
        }

        if(1 == i) p = r;
        else p = r + rho / pre_rho * alpha / omega * (p - omega * v);

        const auto phat = conditioner->apply(p);
        v = system->evaluate(phat);
        alpha = rho / arma::dot(initial_r, v);
        const Col<data_t> s = r - alpha * v;

        suanpan_debug("BiCGSTAB solver local residual: %.4e.\n", residual = arma::norm(s) / norm_b);
        if(residual < setting.tolerance) {
            x += alpha * phat;
            setting.tolerance = residual;
            setting.max_iteration = i;
            return SUANPAN_SUCCESS;
        }

        const auto shat = conditioner->apply(s);
        const Col<data_t> t = system->evaluate(shat);
        omega = arma::dot(t, s) / arma::dot(t, t);
        x += alpha * phat + omega * shat;
        r = s - omega * t;

        pre_rho = rho;

        suanpan_debug("BiCGSTAB solver local residual: %.4e.\n", residual = arma::norm(r) / norm_b);
        if(residual < setting.tolerance) {
            setting.tolerance = residual;
            setting.max_iteration = i;
            return SUANPAN_SUCCESS;
        }

        if(suanpan::approx_equal(omega, ZERO)) {
            setting.tolerance = residual;
            setting.max_iteration = i;
            return SUANPAN_FAIL;
        }
    }

    setting.tolerance = residual;
    return SUANPAN_FAIL;
}

#endif
