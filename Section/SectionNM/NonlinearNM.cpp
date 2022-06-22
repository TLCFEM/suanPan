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
#include <Toolbox/utility.h>

bool NonlinearNM::update_nodal_quantity(mat& jacobian, vec& residual, const double gm, const vec& q, const vec& b, const double alpha, const vec& trial_q, const vec& bn) const {
    jacobian = eye(j_size, j_size);
    residual = zeros(j_size);

    const vec s = q - b;
    const auto f = compute_f(s, alpha);

    if(gm <= datum::eps && f < 0.) return false;

    const auto g = compute_df(s, alpha);
    const double n = norm(g);
    const vec z = g / n;

    const vec ez = diagmat(elastic_diag) * z;
    const mat gedz = gm / n * (diagmat(elastic_diag) - diagmat(elastic_diag) * z * z.t()) * compute_ddf(s, alpha);

    residual(sa) = q - trial_q + gm * ez;
    residual(sc).fill(f);

    jacobian(sa, sa) += gedz;
    jacobian(sa, sc) += ez;

    jacobian(sc, sa) = g.t();
    jacobian(sc, sc).fill(compute_dh(s, alpha));

    if(has_kinematic) {
        residual(sb) = b - bn - kinematic_modulus * gm * ez;

        jacobian(sa, sb) -= gedz;

        jacobian(sb, sa) -= kinematic_modulus * gedz;
        jacobian(sb, sb) += kinematic_modulus * gedz;
        jacobian(sb, sc) -= kinematic_modulus * ez;

        jacobian(sc, sb) = -g.t();
    }

    return true;
}

NonlinearNM::NonlinearNM(const unsigned T, const double EEA, const double EEIS, const double KK, const double LD)
    : DataNonlinearNM{EEA, EEIS, 0., KK}
    , SectionNM(T, SectionType::NM2D)
    , has_kinematic(!suanpan::approx_equal(kinematic_modulus, 0.))
    , si{0, 1}
    , sj{0, 2}
    , sa{0, 1}
    , sb{2, 3}
    , sc{has_kinematic ? 4u : 2u}
    , n_size(2)
    , j_size(has_kinematic ? 5 : 3)
    , elastic_diag{EA, EIS}
    , border([&] {
        mat t_border = zeros(g_size);
        t_border(j_size) = -(t_border(0) = 1.);
        return t_border;
    }())
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
    , rabbit([&] {
        mat scale(3, 3, fill::zeros);
        scale(1, 1) = scale(2, 2) = 4.;
        scale(0, 0) = scale(1, 2) = scale(2, 1) = 2.;

        mat left(g_size, 3, fill::zeros);

        left.rows(sa) = diagmat(elastic_diag) * ti.t() * scale;
        left.rows(sa + j_size) = diagmat(elastic_diag) * tj.t() * scale;

        return left;
    }()) { access::rw(linear_density) = LD; }

NonlinearNM::NonlinearNM(const unsigned T, const double EEA, const double EEIS, const double EEIW, const double KK, const double LD)
    : DataNonlinearNM{EEA, EEIS, EEIW, KK}
    , SectionNM(T, SectionType::NM3D)
    , has_kinematic(!suanpan::approx_equal(kinematic_modulus, 0.))
    , si{0, 1, 3}
    , sj{0, 2, 4}
    , sa{0, 1, 2}
    , sb{3, 4, 5}
    , sc{has_kinematic ? 6u : 3u}
    , n_size(3)
    , j_size(has_kinematic ? 7 : 4)
    , elastic_diag{EA, EIS, EIW}
    , border([&] {
        mat t_border = zeros(g_size);
        t_border(j_size) = -(t_border(0) = 1.);
        return t_border;
    }())
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
    , rabbit([&] {
        mat scale(5, 5, fill::zeros);
        scale(1, 1) = scale(2, 2) = scale(3, 3) = scale(4, 4) = 4.;
        scale(0, 0) = scale(1, 2) = scale(2, 1) = scale(3, 4) = scale(4, 3) = 2.;

        mat left(g_size, 5, fill::zeros);

        left.rows(sa) = diagmat(elastic_diag) * ti.t() * scale;
        left.rows(sa + j_size) = diagmat(elastic_diag) * tj.t() * scale;

        return left;
    }()) { access::rw(linear_density) = LD; }

int NonlinearNM::initialize(const shared_ptr<DomainBase>&) {
    if(SectionType::NM2D == section_type) {
        initial_stiffness.zeros(3, 3);

        initial_stiffness(1, 1) = initial_stiffness(2, 2) = 2. * (initial_stiffness(1, 2) = initial_stiffness(2, 1) = 2. * EIS);
    }
    else {
        initial_stiffness.zeros(6, 6);

        initial_stiffness(5, 5) = 1E3 * EA;

        initial_stiffness(1, 1) = initial_stiffness(2, 2) = 2. * (initial_stiffness(1, 2) = initial_stiffness(2, 1) = 2. * EIS);
        initial_stiffness(3, 3) = initial_stiffness(4, 4) = 2. * (initial_stiffness(3, 4) = initial_stiffness(4, 3) = 2. * EIW);
    }

    initial_stiffness(0, 0) = EA;

    trial_stiffness = current_stiffness = initial_stiffness;

    initialize_history(2 * n_size + 4);

    return SUANPAN_SUCCESS;
}

int NonlinearNM::update_trial_status(const vec& t_deformation) {
    const vec incre_deformation = (trial_deformation = t_deformation) - current_deformation;

    if(norm(incre_deformation) <= datum::eps) return SUANPAN_SUCCESS;

    trial_resistance = current_resistance + (trial_stiffness = initial_stiffness) * incre_deformation;

    trial_history = current_history;
    const vec current_betai(&current_history(0), n_size, false, true);
    const vec current_betaj(&current_history(n_size), n_size, false, true);

    vec betai(&trial_history(0), n_size, false, true);
    vec betaj(&trial_history(n_size), n_size, false, true);
    auto& alphai = trial_history(2llu * n_size);
    auto& alphaj = trial_history(2llu * n_size + 1);
    auto& flagi = trial_history(2llu * n_size + 2);
    auto& flagj = trial_history(2llu * n_size + 3);

    vec qi = trial_resistance(si);
    vec qj = trial_resistance(sj);

    mat jacobian(g_size, g_size, fill::zeros);
    vec residual(g_size), gamma(2, fill::zeros);

    auto counter = 0u;
    auto ref_error = 1.;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("NonlinearNM2D cannot converge within %u iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        mat t_jacobian;
        vec t_residual;

        if(update_nodal_quantity(t_jacobian, t_residual, gamma(0), qi, betai, alphai, trial_resistance(si), current_betai)) flagi = 1.;
        jacobian(0, 0, size(j_size, j_size)) = t_jacobian;
        residual.head(j_size) = t_residual;

        if(update_nodal_quantity(t_jacobian, t_residual, gamma(1), qj, betaj, alphaj, trial_resistance(sj), current_betaj)) flagj = 1.;
        jacobian(j_size, j_size, size(j_size, j_size)) = t_jacobian;
        residual.tail(j_size) = t_residual;

        const vec ra = solve(jacobian, residual);
        const vec rb = solve(jacobian, border);

        const vec incre = ra - dot(border, ra) / dot(border, rb) * rb;

        auto error = norm(residual);
        if(1 == counter) ref_error = std::max(1., error);
        suanpan_debug("NonlinearNM2D local iteration error: %.5E.\n", error /= ref_error);
        if(error <= tolerance && norm(incre) <= tolerance) break;

        qi -= incre(sa);
        qj -= incre(sa + j_size);
        if(has_kinematic) {
            betai -= incre(sb);
            betaj -= incre(sb + j_size);
        }
        gamma(0) -= incre(sc(0));
        gamma(1) -= incre(sc(0) + j_size);
        alphai -= incre(sc(0));
        alphaj -= incre(sc(0) + j_size);
    }

    if(const mat right = [&] {
        const mat ra = solve(jacobian, rabbit);
        const vec rb = solve(jacobian, border);
        return (ra - rb * border.t() * ra / dot(rb, border)).eval();
    }(); SectionType::NM2D == section_type) {
        trial_resistance = ti * qi + tj * qj;

        trial_stiffness = ti * right.rows(sa) + tj * right.rows(sa + j_size);
    }
    else {
        trial_resistance.head(5) = ti * qi + tj * qj;

        trial_stiffness(0, 0, size(5, 5)) = ti * right.rows(sa) + tj * right.rows(sa + j_size);
    }

    return SUANPAN_SUCCESS;
}

vector<vec> NonlinearNM::record(const OutputType P) {
    if(P == OutputType::YF) return {current_history.tail(2)};

    return Section::record(P);
}
