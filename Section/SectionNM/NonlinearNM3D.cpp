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

// ReSharper disable IdentifierTypo
#include "NonlinearNM3D.h"
#include <Recorder/OutputType.h>
#include <Toolbox/utility.h>

const uvec NonlinearNM3D::sa{0, 1, 2, 3, 4};
const uvec NonlinearNM3D::sb{5, 6, 7, 8, 9};
const uvec NonlinearNM3D::sai{0, 1, 3};
const uvec NonlinearNM3D::saj{0, 2, 4};
const uvec NonlinearNM3D::sbi = sai + 5;
const uvec NonlinearNM3D::sbj = saj + 5;

void NonlinearNM3D::initialize_weight(const vec&, const double k) {
    plastic_weight = eye(5, 5);
    kin_weight = k * initial_stiffness(sa, sa);
}

NonlinearNM3D::NonlinearNM3D(const unsigned T, const double EEA, const double EEIS, const double EEIW, const double KK, const double LD)
    : SectionNM3D(T, EEA, EEIS, EEIW, LD)
    , has_kinematic(!suanpan::approx_equal(KK, 0.)) { initialize_history(9); }

int NonlinearNM3D::update_trial_status(const vec& t_deformation) {
    const vec incre_deformation = (trial_deformation = t_deformation) - current_deformation;

    if(norm(incre_deformation) <= datum::eps) return SUANPAN_SUCCESS;

    trial_resistance = current_resistance + (trial_stiffness = initial_stiffness) * incre_deformation;

    trial_history = current_history;
    const vec current_back_resistance(&current_history(0), 5, false, true);
    const auto& current_alphai = current_history(5);
    const auto& current_alphaj = current_history(6);
    vec back_resistance(&trial_history(0), 5, false, true);
    auto& alphai = trial_history(5);
    auto& alphaj = trial_history(6);
    auto& flagi = trial_history(7);
    auto& flagj = trial_history(8);

    auto fi = 0., fj = 0., hi = 0., hj = 0.;
    vec gamma(2, fill::zeros), gi, gj, eni, enj, keni, kenj;
    mat gedni, gednj, gkedni, gkednj;

    const auto update_variable = [&] {
        const vec bi = trial_resistance(sai) - back_resistance(sai);
        const vec bj = trial_resistance(saj) - back_resistance(saj);

        fi = compute_f(bi);
        fj = compute_f(bj);
        gi = compute_df(bi);
        gj = compute_df(bj);
        const vec pi = plastic_weight.cols(sai) * gi;
        const vec pj = plastic_weight.cols(saj) * gj;
        const auto mi = norm(pi);
        const auto mj = norm(pj);
        const vec ni = pi / mi;
        const vec nj = pj / mj;
        const mat pni = gamma(0) / mi * (eye(5, 5) - ni * ni.t()) * plastic_weight.cols(sai) * compute_ddf(bi);
        const mat pnj = gamma(1) / mj * (eye(5, 5) - nj * nj.t()) * plastic_weight.cols(saj) * compute_ddf(bj);
        eni = initial_stiffness(sa, sa) * ni;
        enj = initial_stiffness(sa, sa) * nj;
        gedni = initial_stiffness(sa, sa) * pni;
        gednj = initial_stiffness(sa, sa) * pnj;
        keni = kin_weight * ni;
        kenj = kin_weight * nj;
        gkedni = kin_weight * pni;
        gkednj = kin_weight * pnj;
        alphai = current_alphai + std::max(0., gamma(0));
        alphaj = current_alphaj + std::max(0., gamma(1));
        hi = compute_h(alphai);
        hj = compute_h(alphaj);
    };

    update_variable();

    if(fi < hi && fj < hj) return SUANPAN_SUCCESS;

    const vec trial_q = trial_resistance;

    vec residual(local_size), incre;
    mat jacobian;

    auto counter = 0;
    auto ref_error = 1.;
    auto reduced = false;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("NonlinearNM3D cannot converge within %u iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        jacobian = eye(local_size, local_size);

        // resistance
        residual(sa) = trial_resistance(sa) - trial_q(sa);

        // back resistance
        if(has_kinematic) residual(sb) = back_resistance - current_back_resistance;

        // gamma
        if(gamma(0) <= datum::eps && fi < hi) residual(si).fill(0.);
        else {
            flagi = 1.;

            // resistance
            residual(sa) += gamma(0) * eni;

            jacobian(sa, sai) += gedni;
            jacobian(sa, si) += eni;

            // back resistance
            if(has_kinematic) {
                jacobian(sa, sbi) -= gedni;

                residual(sb) -= gamma(0) * keni;

                jacobian(sb, sai) -= gkedni;
                jacobian(sb, sbi) += gkedni;
                jacobian(sb, si) -= keni;
            }

            residual(si).fill(fi - hi);

            jacobian(si, sai) = gi.t();
            if(has_kinematic) jacobian(si, sbi) = -gi.t();
            jacobian(si, si).fill(-compute_dh(alphai));
        }

        if(gamma(1) <= datum::eps && fj < hj) residual(sj).fill(0.);
        else {
            flagj = 1.;

            // resistance
            residual(sa) += gamma(1) * enj;

            jacobian(sa, saj) += gednj;
            jacobian(sa, sj) += enj;

            // back resistance
            if(has_kinematic) {
                jacobian(sa, sbj) -= gednj;

                residual(sb) -= gamma(1) * kenj;

                jacobian(sb, saj) -= gkednj;
                jacobian(sb, sbj) += gkednj;
                jacobian(sb, sj) -= kenj;
            }

            residual(sj).fill(fj - hj);

            jacobian(sj, saj) = gj.t();
            if(has_kinematic) jacobian(sj, sbj) = -gj.t();
            jacobian(sj, sj).fill(-compute_dh(alphaj));
        }

        if(local_size == rank(jacobian) + 1llu) {
            reduced = true;
            if(!solve(incre, jacobian(spa, spa), residual.head(local_size - 1llu))) return SUANPAN_FAIL;
            incre.resize(local_size);
            incre(si) *= .5;
            incre(sj) = incre(si);
        }
        else {
            reduced = false;
            if(!solve(incre, jacobian, residual)) return SUANPAN_FAIL;
        }

        auto error = norm(residual);

        if(1 == counter) ref_error = std::max(1., error);
        suanpan_debug("NonlinearNM3D local iteration error: %.5E.\n", error /= ref_error);
        if(error <= tolerance && norm(incre) <= tolerance) break;

        trial_resistance.head(5) -= incre(sa);
        if(has_kinematic) back_resistance -= incre(sb);
        gamma -= incre.tail(2);

        update_variable();
    }

    if(reduced) trial_stiffness(sa, sa) = solve(jacobian(spa, spa), resize(initial_stiffness(sa, sa), local_size - 1llu, 5)).eval().head_rows(5);
    else trial_stiffness(sa, sa) = solve(jacobian, resize(initial_stiffness(sa, sa), local_size, 5)).eval().rows(sa);

    return SUANPAN_SUCCESS;
}

vector<vec> NonlinearNM3D::record(const OutputType P) {
    if(P == OutputType::YF) return {current_history.tail(2)};

    return Section::record(P);
}
