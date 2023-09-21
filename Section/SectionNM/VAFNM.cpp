/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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

#include "VAFNM.h"
#include <Toolbox/utility.h>

int VAFNM::compute_local_integration(vec& q, mat& jacobian) {
    trial_history = current_history;
    const vec current_beta(&current_history(0), d_size);
    const auto &ani = current_history(d_size), &anj = current_history(d_size + 1llu);

    vec beta(&trial_history(0), d_size, false, true);
    auto &ai = trial_history(d_size), &aj = trial_history(d_size + 1llu);
    auto &flagi = trial_history(d_size + 2llu), &flagj = trial_history(d_size + 3llu); // yield flag

    const vec trial_q = q = trial_resistance.head(d_size) / yield_diag;

    flagi = 0.;
    flagj = 0.;

    auto gamma = 0.;

    auto counter = 0u;
    auto ref_error = 1.;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        vec z(d_size, fill::zeros);
        mat pzpq(d_size, d_size, fill::zeros);
        vec pzpai(d_size, fill::zeros);
        vec pzpaj(d_size, fill::zeros);
        vec residual(g_size, fill::none);

        jacobian.eye(g_size, g_size);
        residual(ge).fill(0.);

        {
            const vec si = q(ni) - beta(ni), hi = compute_h(ai);
            if(const auto fi = compute_f(si, hi); fi > 0. || static_cast<bool>(flagi)) {
                flagi = 1.;
                residual(ge) += fi;

                const vec g = compute_df(si, hi);
                const vec dh = -si % compute_dh(ai) / hi;
                const mat dg = ti * compute_ddf(si, hi);
                z(ni) += g;
                pzpq += dg * ti.t();
                pzpai = dg * dh;
                jacobian(ge, gc).fill(dot(g, dh));
            }
        }

        {
            const vec sj = q(nj) - beta(nj), hj = compute_h(aj);
            if(const auto fj = compute_f(sj, hj); fj > 0. || static_cast<bool>(flagj)) {
                flagj = 1.;
                residual(ge) += fj;

                const vec g = compute_df(sj, hj);
                const vec dh = -sj % compute_dh(aj) / hj;
                const mat dg = tj * compute_ddf(sj, hj);
                z(nj) += g;
                pzpq += dg * tj.t();
                pzpaj = dg * dh;
                jacobian(ge, gd).fill(dot(g, dh));
            }
        }

        const vec m = normalise(z);

        const auto norm_mi = norm(m(ni));
        const auto norm_mj = norm(m(nj));

        if(1u == counter) {
            gamma = .5 * residual(ge(0)) / dot(m, z);
            q -= gamma * m;
            if(has_kinematic) {
                beta += gamma * kin_modulus % m;
                beta /= 1. + kin_base * gamma;
            }
            ai += gamma * norm_mi;
            aj += gamma * norm_mj;
            continue;
        }

        residual(ga) = q - trial_q + gamma * m;
        residual(gc).fill(ai - ani - gamma * norm_mi);
        residual(gd).fill(aj - anj - gamma * norm_mj);

        jacobian(ga, ge) = m;
        jacobian(gc, ge).fill(-norm_mi);
        jacobian(gd, ge).fill(-norm_mj);
        jacobian(ge, ga) = z.t();
        jacobian(ge, ge).fill(0.);

        mat prpz(g_size, d_size, fill::zeros), dzdx(d_size, g_size, fill::zeros);

        prpz.rows(ga) = gamma * eye(d_size, d_size);
        prpz.rows(gc) = -gamma * (ti * normalise(m(ni))).t();
        prpz.rows(gd) = -gamma * (tj * normalise(m(nj))).t();

        dzdx.cols(ga) = pzpq;
        dzdx.cols(gc) = pzpai;
        dzdx.cols(gd) = pzpaj;

        if(has_kinematic) {
            residual(gb) = (1. + kin_base * gamma) % beta - current_beta - gamma * kin_modulus % m;

            jacobian(gb, gb) += diagmat(kin_base * gamma);
            jacobian(gb, ge) = kin_base % beta - kin_modulus % m;
            jacobian(ge, gb) = -z.t();

            prpz.rows(gb) = diagmat(-gamma * kin_modulus);

            dzdx.cols(gb) = -pzpq;
        }

        const vec incre = solve(jacobian += prpz * (eye(d_size, d_size) - m * m.t()) / norm(z) * dzdx, residual);

        const auto error = inf_norm(incre);
        if(2u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || inf_norm(residual) < tolerance) return SUANPAN_SUCCESS;

        q -= incre(ga);
        if(has_kinematic) beta -= incre(gb);
        ai -= incre(gc(0));
        aj -= incre(gd(0));
        gamma -= incre(ge(0));
    }
}

vec VAFNM::compute_h(const double alpha) const { return {n_size, fill::value(std::max(datum::eps, 1. + iso_modulus * alpha + iso_saturation - iso_saturation * exp(-iso_decay * alpha)))}; }

vec VAFNM::compute_dh(const double alpha) const { return compute_h(alpha).transform([&](const double h) { return suanpan::approx_equal(h, datum::eps) ? 0. : iso_modulus + iso_saturation * iso_decay * exp(-iso_decay * alpha); }); }

VAFNM::VAFNM(const unsigned T, const double EEA, const double EEIS, const double HH, const double HS, const double HD, vec&& KK, vec&& KB, const double LD, vec&& YF)
    : NonlinearNM(T, EEA, EEIS, any(KK) || any(KB), LD, std::forward<vec>(YF))
    , iso_modulus(HH)
    , iso_saturation(HS)
    , iso_decay(HD)
    , kin_modulus{KK(0), KK(1), KK(1)}
    , kin_base{KB(0), KB(1), KB(1)} {}

VAFNM::VAFNM(const unsigned T, const double EEA, const double EEIS, const double EEIW, const double HH, const double HS, const double HD, vec&& KK, vec&& KB, const double LD, vec&& YF)
    : NonlinearNM(T, EEA, EEIS, EEIW, any(KK) || any(KB), LD, std::forward<vec>(YF))
    , iso_modulus(HH)
    , iso_saturation(HS)
    , iso_decay(HD)
    , kin_modulus{KK(0), KK(1), KK(1), KK(2), KK(2)}
    , kin_base{KB(0), KB(1), KB(1), KB(2), KB(2)} {}
