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

#include "LinearHardeningNM.h"
#include <Toolbox/utility.h>

int LinearHardeningNM::compute_local_integration(vec& q, mat& jacobian, const bool yield_flagi, const bool yield_flagj) {
    trial_history = current_history;
    const vec current_beta(&current_history(0), d_size, false, true);
    const auto &ani = current_history(d_size), &anj = current_history(d_size + 1llu);

    vec beta(&trial_history(0), d_size, false, true);
    auto &ai = trial_history(d_size), &aj = trial_history(d_size + 1llu);
    auto &flagi = trial_history(d_size + 2llu), &flagj = trial_history(d_size + 3llu); // yield flag

    flagi = yield_flagi;
    flagj = yield_flagj;

    const vec trial_q = q = trial_resistance.head(d_size) / yield_diag;

    auto gamma = 0.;

    auto counter = 0u;
    auto ref_error = 1.;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("LinearHardeningNM cannot converge within %u iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        vec z(d_size, fill::zeros);
        mat pzpq(d_size, d_size, fill::zeros);
        vec pzpai(d_size, fill::zeros);
        vec pzpaj(d_size, fill::zeros);
        vec residual(g_size, fill::none);

        jacobian.eye(g_size, g_size);
        residual(ge).fill(0.);

        if(yield_flagi) {
            const vec si = q(ni) - beta(ni), hi = compute_h(ai);
            residual(ge) += compute_f(si, hi);

            const vec g = compute_df(si, hi);
            const vec dh = -si % compute_dh(ai) / hi;
            const mat dg = ti * compute_ddf(si, hi);
            z(ni) += g;
            pzpq += dg * ti.t();
            pzpai = dg * dh;
            jacobian(ge, gc).fill(dot(g, dh));
        }

        if(yield_flagj) {
            const vec sj = q(nj) - beta(nj), hj = compute_h(aj);
            residual(ge) += compute_f(sj, hj);

            const vec g = compute_df(sj, hj);
            const vec dh = -sj % compute_dh(aj) / hj;
            const mat dg = tj * compute_ddf(sj, hj);
            z(nj) += g;
            pzpq += dg * tj.t();
            pzpaj = dg * dh;
            jacobian(ge, gd).fill(dot(g, dh));
        }

        if(1 == counter) {
            gamma = residual(ge(0)) / dot(z, z);
            q = trial_q - gamma * z;
            continue;
        }

        const auto norm_zi = norm(z(ni));
        const auto norm_zj = norm(z(nj));

        residual(ga) = q - trial_q + gamma * z;
        residual(gc).fill(ai - ani - gamma * norm_zi);
        residual(gd).fill(aj - anj - gamma * norm_zj);

        jacobian(ga, ge) = z;
        jacobian(gc, ge).fill(-norm_zi);
        jacobian(gd, ge).fill(-norm_zj);
        jacobian(ge, ga) = z.t();
        jacobian(ge, ge).fill(0.);

        mat prpz(g_size, d_size, fill::zeros), dzdx(d_size, g_size, fill::zeros);

        prpz.rows(ga) = gamma * eye(d_size, d_size);
        prpz.rows(gc) = -gamma * (ti * normalise(z(ni))).t();
        prpz.rows(gd) = -gamma * (tj * normalise(z(nj))).t();

        dzdx.cols(ga) = pzpq;
        dzdx.cols(gc) = pzpai;
        dzdx.cols(gd) = pzpaj;

        if(has_kinematic) {
            residual(gb) = beta - current_beta - kinematic_modulus * gamma * z;

            jacobian(gb, ge) = -kinematic_modulus * z;
            jacobian(ge, gb) = -z.t();

            prpz.rows(gb) = -kinematic_modulus * gamma * eye(d_size, d_size);

            dzdx.cols(gb) = -pzpq;
        }

        const vec incre = solve(jacobian += prpz * dzdx, residual);

        auto error = norm(residual);
        if(2 == counter) ref_error = std::max(1., error);
        suanpan_debug("LinearHardeningNM local iteration error: %.5E.\n", error /= ref_error);
        if(norm(incre) <= tolerance && error <= tolerance) return SUANPAN_SUCCESS;

        q -= incre(ga);
        if(has_kinematic) beta -= incre(gb);
        ai -= incre(gc(0));
        aj -= incre(gd(0));
        gamma -= incre(ge(0));
    }
}

vec LinearHardeningNM::compute_h(const double alpha) const { return {n_size, fill::value(std::max(datum::eps, 1. + isotropic_modulus * alpha))}; }

vec LinearHardeningNM::compute_dh(const double alpha) const { return compute_h(alpha).transform([&](const double h) { return suanpan::approx_equal(h, datum::eps) ? 0. : isotropic_modulus; }); }

LinearHardeningNM::LinearHardeningNM(const unsigned T, const double EEA, const double EEIS, const double HH, const double KK, const double LD, vec&& YF)
    : NonlinearNM(T, EEA, EEIS, !suanpan::approx_equal(KK, 0.), LD, std::forward<vec>(YF))
    , isotropic_modulus(HH)
    , kinematic_modulus(KK) {}

LinearHardeningNM::LinearHardeningNM(const unsigned T, const double EEA, const double EEIS, const double EEIW, const double HH, const double KK, const double LD, vec&& YF)
    : NonlinearNM(T, EEA, EEIS, EEIW, !suanpan::approx_equal(KK, 0.), LD, std::forward<vec>(YF))
    , isotropic_modulus(HH)
    , kinematic_modulus(KK) {}
