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
    , has_kinematic(KK)
    , n_size(2)
    , g_size(has_kinematic ? 12 : 9)
    , si{0, 1}
    , sj{0, 2}
    , ga{0, 1, 2}
    , gc{3, 4, 5}
    , gd{6}
    , ge{7}
    , gf{8}
    , gb(has_kinematic ? uvec{9, 10, 11} : uvec{}) { access::rw(linear_density) = LD; }

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
    , has_kinematic(KK)
    , n_size(3)
    , g_size(has_kinematic ? 18 : 13)
    , si{0, 1, 3}
    , sj{0, 2, 4}
    , ga{0, 1, 2, 3, 4}
    , gc{5, 6, 7, 8, 9}
    , gd{10}
    , ge{11}
    , gf{12}
    , gb(has_kinematic ? uvec{13, 14, 15, 16, 17} : uvec{}) { access::rw(linear_density) = LD; }

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

    initialize_history(d_size + 3);

    return SUANPAN_SUCCESS;
}

int NonlinearNM::update_trial_status(const vec& t_deformation) {
    const vec incre_deformation = (trial_deformation = t_deformation) - current_deformation;

    if(norm(incre_deformation) <= datum::eps) return SUANPAN_SUCCESS;

    return SUANPAN_SUCCESS;

    trial_history = current_history;
    const vec current_beta(&current_history(0), d_size - 1, false, true);
    const vec bni = current_beta(si);
    const vec bnj = current_beta(sj);
    const auto& ani = current_history(d_size - 1);
    const auto& anj = current_history(d_size);

    vec beta(&trial_history(0), d_size - 1, false, true);
    auto& alphai = trial_history(d_size - 1);
    auto& alphaj = trial_history(d_size);
    auto& flagi = trial_history(d_size + 1);
    auto& flagj = trial_history(d_size + 2);

    trial_resistance = current_resistance + (trial_stiffness = initial_stiffness) * incre_deformation;

    const vec trial_q = trial_resistance / yield_diag;

    auto compute_nm_surface = [&](const vec& q, const vec& beta, const double alpha) {
        return compute_f(q - beta, compute_h(alpha));
    };

    const auto fi = compute_nm_surface(trial_q(si), beta(si), alphai);
    const auto fj = compute_nm_surface(trial_q(sj), beta(sj), alphaj);

    if(fi < 0. && fj < 0.) return SUANPAN_SUCCESS;

    vec q = trial_q;

    mat jacobian;
    vec residual(g_size, fill::none), e(d_size - 1, fill::zeros);
    auto gamma = 0.;

    auto compute_yield_function = [](const double fi, const double fj) {
        return std::max(0., fi) + std::max(0., fj);
    };

    auto counter = 0u;
    auto ref_error = 1.;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("NonlinearNM cannot converge within %u iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        jacobian.eye(g_size, g_size);

        vec z(d_size - 1, fill::zeros);

        residual(ga) = q - trial_q + e;
        if(has_kinematic) residual(gb) = beta - current_beta - e;
        residual(gc) = e - gamma * z;
        residual(gd).fill(alphai - gamma * norm(ti * z));
        residual(ge).fill(alphaj - gamma * norm(tj * z));
        residual(gf).fill(compute_yield_function(fi, fj));

        const vec incre = solve(jacobian, residual);

        auto error = norm(residual);
        if(1 == counter) ref_error = std::max(1., error);
        suanpan_debug("NonlinearNM local iteration error: %.5E.\n", error /= ref_error);
        if(norm(incre) <= tolerance && error <= tolerance) break;

        q -= incre(ga);
        if(has_kinematic) beta -= incre(gb);
        e -= incre(gc);
        alphai -= incre(gd(0));
        alphaj -= incre(ge(0));
        gamma -= incre(gf(0));
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
    if(P == OutputType::HIST) return {current_history};

    return Section::record(P);
}
