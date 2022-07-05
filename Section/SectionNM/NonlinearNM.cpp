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
        transi(0, 0) = transi(1, 1) = 1.;
        return transi;
    }())
    , tj([] {
        mat transj(3, 2, fill::zeros);
        transj(0, 0) = transj(2, 1) = 1.;
        return transj;
    }())
    , has_kinematic(KK)
    , n_size(2)
    , g_size(has_kinematic ? 9 : 6)
    , ni{0, 1}
    , nj{0, 2}
    , ga{0, 1, 2}
    , gb(has_kinematic ? uvec{6, 7, 8} : uvec{})
    , gc{3}
    , gd{4}
    , ge{5} { access::rw(linear_density) = LD; }

NonlinearNM::NonlinearNM(const unsigned T, const double EEA, const double EEIS, const double EEIW, const bool KK, const double LD, vec&& YF)
    : DataNonlinearNM{EEA, EEIS, EEIW, std::forward<vec>(YF)}
    , SectionNM(T, SectionType::NM3D)
    , yield_diag{yield_force(0), yield_force(1), yield_force(1), yield_force(2), yield_force(2)}
    , ti([] {
        mat transi(5, 3, fill::zeros);
        transi(0, 0) = transi(1, 1) = transi(3, 2) = 1.;
        return transi;
    }())
    , tj([] {
        mat transj(5, 3, fill::zeros);
        transj(0, 0) = transj(2, 1) = transj(4, 2) = 1.;
        return transj;
    }())
    , has_kinematic(KK)
    , n_size(3)
    , g_size(has_kinematic ? 13 : 8)
    , ni{0, 1, 3}
    , nj{0, 2, 4}
    , ga{0, 1, 2, 3, 4}
    , gb(has_kinematic ? uvec{8, 9, 10, 11, 12} : uvec{})
    , gc{5}
    , gd{6}
    , ge{7} { access::rw(linear_density) = LD; }

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

    initialize_history(d_size + 4);

    return SUANPAN_SUCCESS;
}

int NonlinearNM::update_trial_status(const vec& t_deformation) {
    const vec incre_deformation = (trial_deformation = t_deformation) - current_deformation;

    if(norm(incre_deformation) <= datum::eps) return SUANPAN_SUCCESS;

    trial_resistance = current_resistance + (trial_stiffness = initial_stiffness) * incre_deformation;

    const vec trial_q = trial_resistance.head(d_size) / yield_diag;

    const vec current_beta(&current_history(0), d_size, false, true);
    const vec bni = current_beta(ni), bnj = current_beta(nj);
    const auto &ani = current_history(d_size), &anj = current_history(d_size + 1llu);

    const auto fi = compute_f(trial_q(ni) - bni, compute_h(ani));
    const auto fj = compute_f(trial_q(nj) - bnj, compute_h(anj));

    if(fi <= 0. && fj <= 0.) return SUANPAN_SUCCESS;

    vec q;
    mat jacobian;

    if(SUANPAN_SUCCESS != compute_local_integration(q, jacobian, fi > 0., fj > 0.)) return SUANPAN_FAIL;

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
    if(P == OutputType::HIST) return {current_history};

    return Section::record(P);
}
