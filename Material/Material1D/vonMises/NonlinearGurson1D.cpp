/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
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

#include "NonlinearGurson1D.h"

#include <Recorder/OutputType.h>
#include <Toolbox/tensor.h>

NonlinearGurson1D::NonlinearGurson1D(const unsigned T, const double E, const double V, const double Q1, const double Q2, const double FN, const double SN, const double EN, const double R)
    : DataNonlinearGurson1D{E, V, Q1, Q2, FN, SN, EN}
    , Material1D(T, R) {}

int NonlinearGurson1D::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    initialize_history(2);

    return SUANPAN_SUCCESS;
}

double NonlinearGurson1D::get_parameter(const ParameterType P) const { return material_property(elastic_modulus, poissons_ratio)(P); }

int NonlinearGurson1D::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    trial_history = current_history;
    auto& pe = trial_history(0); // equivalent plastic strain
    auto& f = trial_history(1);  // volume fraction
    const auto& current_pe = current_history(0);
    const auto& current_f = current_history(1);

    const auto& trial_s = trial_stress(0);
    auto s = trial_s;

    mat jacobian(4, 4);
    vec incre, residual(4);
    auto gamma = 0.;

    auto counter = 0u;
    auto ref_error = 1.;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        const auto hardening = compute_hardening(pe);
        const auto &k = hardening(0), &dk = hardening(1);
        const auto hyper_term = .5 * q2 * s / k;
        const auto cosh_term = std::cosh(hyper_term);
        const auto sinh_term = std::sinh(hyper_term);
        const auto an = para_b * std::exp(-.5 * std::pow((pe - en) / sn, 2.));

        const auto diff_pe = pe - current_pe, diff_s = s - trial_s;

        residual(0) = s * s + k * k * (f * q1 * (2. * cosh_term - q1 * f) - 1.);

        if(1 == counter && residual(0) < 0.) return SUANPAN_SUCCESS;

        residual(1) = (1. - f) * k * diff_pe - 2. * gamma * s * s + s * diff_s / nine_bulk;
        residual(2) = f - current_f + (1. - f) * diff_s / three_bulk - an * diff_pe;
        residual(3) = diff_s + nine_bulk * gamma * q1 * q2 * f * k * sinh_term;

        jacobian(0, 0) = 0.;
        jacobian(0, 1) = (4. * q1 * f * k * cosh_term - q1 * q2 * f * s * sinh_term - 2. * k * (q1 * q1 * f * f + 1.)) * dk;
        jacobian(0, 2) = 2. * k * k * q1 * (cosh_term - q1 * f);
        jacobian(0, 3) = 2. * s + q1 * q2 * f * k * sinh_term;
        jacobian(1, 0) = -2. * s * s;
        jacobian(1, 1) = (1. - f) * (dk * diff_pe + k);
        jacobian(1, 2) = -k * diff_pe;
        jacobian(1, 3) = (s + diff_s) / nine_bulk - 4. * s * gamma;
        jacobian(2, 0) = 0.;
        jacobian(2, 1) = an / sn / sn * (pe - en) * diff_pe - an;
        jacobian(2, 2) = 1. - diff_s / three_bulk;
        jacobian(2, 3) = (1. - f) / three_bulk;
        jacobian(3, 0) = nine_bulk * q1 * q2 * f * k * sinh_term;
        jacobian(3, 1) = nine_bulk * q1 * q2 * gamma * f * (sinh_term - hyper_term * cosh_term) * dk;
        jacobian(3, 2) = nine_bulk * q1 * q2 * gamma * k * sinh_term;
        jacobian(3, 3) = 1. + .5 * nine_bulk * q1 * q2 * q2 * gamma * f * cosh_term;

        if(!solve(incre, jacobian, residual)) return SUANPAN_FAIL;

        const auto error = inf_norm(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || ((error < tolerance || inf_norm(residual) < tolerance) && counter > 5u)) break;

        gamma -= incre(0);
        pe = std::max(pe - incre(1), 0.);
        f = suanpan::clamp_unit(f - incre(2)); // avoid overshoot
        s -= incre(3);
    }

    vec left, right(4);

    right(0) = 0.;
    right(1) = s / nine_bulk;
    right(2) = (1. - f) / three_bulk;
    right(3) = 1.;

    if(!solve(left, jacobian, right)) return SUANPAN_FAIL;

    trial_stress = s;

    trial_stiffness = left(3) * elastic_modulus;

    return SUANPAN_SUCCESS;
}

int NonlinearGurson1D::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int NonlinearGurson1D::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int NonlinearGurson1D::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

std::vector<vec> NonlinearGurson1D::record(const OutputType P) {
    if(P == OutputType::VF) return {vec{current_history(1)}};

    return Material1D::record(P);
}
