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
#include "NonlinearNM.h"
#include <Recorder/OutputType.h>

NonlinearNM::NonlinearNM(const unsigned T, const double EEA, const double EEIS, const bool KK, const double LD, vec&& YF)
    : DataNonlinearNM{EEA, EEIS, 0., std::forward<vec>(YF)}
    , SectionNM(T, SectionType::NM2D)
    , yield_diag{yield_force(0), yield_force(1), yield_force(1)}
    , ti([] {
        mat transi(3, 2, fill::zeros);
        transi(0, 0) = .5;
        transi(1, 1) = 1.;
        return transi;
    }())
    , tj([] {
        mat transj(3, 2, fill::zeros);
        transj(0, 0) = .5;
        transj(2, 1) = 1.;
        return transj;
    }())
    , si{0, 1}
    , sj{0, 2}
    , gi(KK ? uvec{0, 1, 3, 4, 6} : uvec{0, 1, 3})
    , gj(KK ? uvec{0, 2, 3, 5, 7} : uvec{0, 2, 4})
    , ga{0, 1, 2}
    , gb{3, 4, 5}
    , gc(KK ? uvec{6, 7} : uvec{3, 4})
    , g_size(KK ? 8 : 5)
    , has_kinematic(KK)
    , sa{0, 1}
    , sb{2, 3}
    , sc{has_kinematic ? 4u : 2u}
    , n_size(2)
    , j_size(has_kinematic ? 5 : 3) { access::rw(linear_density) = LD; }

NonlinearNM::NonlinearNM(const unsigned T, const double EEA, const double EEIS, const double EEIW, const bool KK, const double LD, vec&& YF)
    : DataNonlinearNM{EEA, EEIS, EEIW, std::forward<vec>(YF)}
    , SectionNM(T, SectionType::NM3D)
    , yield_diag{yield_force(0), yield_force(1), yield_force(1), yield_force(2), yield_force(2)}
    , ti([] {
        mat transi(5, 3, fill::zeros);
        transi(0, 0) = .5;
        transi(1, 1) = transi(3, 2) = 1.;
        return transi;
    }())
    , tj([] {
        mat transj(5, 3, fill::zeros);
        transj(0, 0) = .5;
        transj(2, 1) = transj(4, 2) = 1.;
        return transj;
    }())
    , si{0, 1, 3}
    , sj{0, 2, 4}
    , gi(KK ? uvec{0, 1, 3, 5, 6, 8, 10} : uvec{0, 1, 3, 5})
    , gj(KK ? uvec{0, 2, 4, 5, 7, 9, 11} : uvec{0, 2, 4, 6})
    , ga{0, 1, 2, 3, 4}
    , gb{5, 6, 7, 8, 9}
    , gc(KK ? uvec{10, 11} : uvec{5, 6})
    , g_size(KK ? 12 : 7)
    , has_kinematic(KK)
    , sa{0, 1, 2}
    , sb{3, 4, 5}
    , sc{has_kinematic ? 6u : 3u}
    , n_size(3)
    , j_size(has_kinematic ? 7 : 4) { access::rw(linear_density) = LD; }

int NonlinearNM::initialize(const shared_ptr<DomainBase>&) {
    if(SectionType::NM2D == section_type) initial_stiffness.zeros(3, 3);
    else {
        initial_stiffness.zeros(6, 6);

        initial_stiffness(5, 5) = 1E3 * EA;

        initial_stiffness(3, 3) = initial_stiffness(4, 4) = 2. * (initial_stiffness(3, 4) = initial_stiffness(4, 3) = 2. * EIW);
    }

    initial_stiffness(0, 0) = EA;

    initial_stiffness(1, 1) = initial_stiffness(2, 2) = 2. * (initial_stiffness(1, 2) = initial_stiffness(2, 1) = 2. * EIS);

    trial_stiffness = current_stiffness = initial_stiffness;

    initialize_history(2 * n_size + 3);

    return SUANPAN_SUCCESS;
}

int NonlinearNM::update_trial_status(const vec& t_deformation) {
    const vec incre_deformation = (trial_deformation = t_deformation) - current_deformation;

    if(norm(incre_deformation) <= datum::eps) return SUANPAN_SUCCESS;

    trial_history = current_history;
    const vec current_beta(&current_history(0), 2llu * n_size - 1, false, true);

    vec beta(&trial_history(0), 2llu * n_size - 1, false, true);
    vec alpha(&trial_history(2llu * n_size - 1), 2, false, true);
    auto& flagi = trial_history(2llu * n_size + 1);
    auto& flagj = trial_history(2llu * n_size + 2);

    trial_resistance = current_resistance + initial_stiffness * incre_deformation;

    const vec trial_q = ti * (trial_resistance(si) / yield_force) + tj * (trial_resistance(sj) / yield_force);
    vec q = trial_q;

    mat jacobian;
    vec residual(g_size, fill::none), gamma(2, fill::zeros);

    auto counter = 0u;
    auto ref_error = 1.;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("NonlinearNM cannot converge within %u iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        jacobian.eye(g_size, g_size);

        residual(ga) = q - trial_q;
        if(has_kinematic) residual(gb) = beta - current_beta;
        residual(gc).fill(0.);

        mat t_jacobian;
        vec t_residual;

        if(update_nodal_quantity(t_jacobian, t_residual, gamma(0), q(si), beta(si), alpha(0))) {
            flagi = 1.;

            jacobian(gc(0), gc(0)) = 0.;
            jacobian(gi, gi) += t_jacobian;
            residual(gi) += t_residual;
        }
        if(update_nodal_quantity(t_jacobian, t_residual, gamma(1), q(sj), beta(sj), alpha(1))) {
            flagj = 1.;

            jacobian(gc(1), gc(1)) = 0.;
            jacobian(gj, gj) += t_jacobian;
            residual(gj) += t_residual;
        }

        const vec incre = solve(jacobian, residual);

        auto error = norm(residual);
        if(1 == counter) ref_error = std::max(1., error);
        suanpan_debug("NonlinearNM local iteration error: %.5E.\n", error /= ref_error);
        if(norm(incre) <= tolerance && error <= tolerance) break;

        q -= incre(ga);
        if(has_kinematic) beta -= incre(gb);
        gamma -= incre(gc);
        alpha -= incre(gc);
    }

    if(SectionType::NM2D == section_type) {
        trial_resistance = yield_diag % q;

        mat left(g_size, 3, fill::zeros);
        left.rows(ga) = diagmat(1. / yield_diag) * initial_stiffness;

        trial_stiffness = diagmat(yield_diag) * solve(jacobian, left).eval().head_rows(3);
    }
    else {
        trial_resistance.head(5) = yield_diag % q;

        mat left(g_size, 5, fill::zeros);
        left.rows(ga) = diagmat(1. / yield_diag) * initial_stiffness(0, 0, size(5, 5));

        trial_stiffness(0, 0, size(5, 5)) = diagmat(yield_diag) * solve(jacobian, left).eval().head_rows(5);
    }

    return SUANPAN_SUCCESS;
}

vector<vec> NonlinearNM::record(const OutputType P) {
    if(P == OutputType::YF) return {current_history.tail(2)};

    return Section::record(P);
}
