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

#include "NonlinearHoffman.h"
#include <Toolbox/tensor.h>

constexpr double NonlinearHoffman::four_third = 4. / 3.;
const double NonlinearHoffman::root_two_third = sqrt(.5 * four_third);
constexpr unsigned NonlinearHoffman::max_iteration = 20;

double NonlinearHoffman::compute_yield_function(const vec& t_stress) const {
    const auto &S1 = t_stress(0), &S2 = t_stress(1), &S3 = t_stress(2);
    const auto &S4 = t_stress(3), &S5 = t_stress(4), &S6 = t_stress(5);

    const auto T1 = pow(S1 - S2, 2.), T2 = pow(S2 - S3, 2.), T3 = pow(S3 - S1, 2.);

    return C1 * T1 + C2 * T2 + C3 * T3 + C4 * S4 * S4 + C5 * S5 * S5 + C6 * S6 * S6 + C7 * S1 + C8 * S2 + C9 * S3;
}

NonlinearHoffman::NonlinearHoffman(const unsigned T, vec&& EE, vec&& VV, vec&& SS, const double R)
    : DataNonlinearHoffman{std::forward<vec>(EE), std::forward<vec>(VV), std::forward<vec>(SS)}
    , Material3D(T, R) {
    // S(0) = \sigma_{11}^t    S(1) = \sigma_{11}^c
    // S(2) = \sigma_{22}^t    S(3) = \sigma_{22}^c
    // S(4) = \sigma_{33}^t    S(5) = \sigma_{33}^c
    // S(6) = \sigma_{12}^0    S(7) = \sigma_{23}^0    S(8) = \sigma_{13}^0

    proj_a.zeros(6, 6);
    proj_b.zeros(6, 1);

    const auto T1 = 1. / yield_stress(0) / yield_stress(1);
    const auto T2 = 1. / yield_stress(2) / yield_stress(3);
    const auto T3 = 1. / yield_stress(4) / yield_stress(5);

    proj_b(0) = C7 = (yield_stress(1) - yield_stress(0)) * (proj_a(0, 0) = T1);
    proj_b(1) = C8 = (yield_stress(3) - yield_stress(2)) * (proj_a(1, 1) = T2);
    proj_b(2) = C9 = (yield_stress(5) - yield_stress(4)) * (proj_a(2, 2) = T3);

    proj_a(0, 1) = proj_a(1, 0) = -(C1 = .5 * (T1 + T2 - T3));
    proj_a(1, 2) = proj_a(2, 1) = -(C2 = .5 * (T2 + T3 - T1));
    proj_a(2, 0) = proj_a(0, 2) = -(C3 = .5 * (T3 + T1 - T2));
    proj_a(3, 3) = C4 = 1. / yield_stress(6) / yield_stress(6);
    proj_a(4, 4) = C5 = 1. / yield_stress(7) / yield_stress(7);
    proj_a(5, 5) = C6 = 1. / yield_stress(8) / yield_stress(8);
    proj_a *= 2.;
}

int NonlinearHoffman::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::orthotropic_stiffness(modulus, ratio);

    inv_stiffness = inv(initial_stiffness);

    initialize_history(1);

    return SUANPAN_SUCCESS;
}

double NonlinearHoffman::get_parameter(const ParameterType P) const {
    if(ParameterType::DENSITY == P) return density;
    return 0.;
}

int NonlinearHoffman::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= tolerance) return SUANPAN_SUCCESS;

    trial_history = current_history;
    auto& plastic_strain = trial_history(0);
    const auto& current_plastic_strain = current_history(0);

    const vec e_strain = incre_strain + solve(initial_stiffness, current_stress);

    auto gamma = 0., ref_error = 1.;

    auto counter = 0u;
    while(true) {
        if(++counter == max_iteration) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        const mat hessian = inv(inv_stiffness + gamma * proj_a);
        const vec n = proj_a * (trial_stress = hessian * (e_strain - gamma * proj_b)) + proj_b;
        const auto eta = root_two_third * tensor::strain::norm(n);
        const auto k = compute_k(plastic_strain = current_plastic_strain + gamma * eta);
        const auto residual = compute_yield_function(trial_stress) - k * k;

        if(1u == counter && residual < 0.) {
            trial_stiffness = initial_stiffness;
            return SUANPAN_SUCCESS;
        }

        const auto dk = compute_dk(plastic_strain);
        const auto factor = 2. * k * dk * eta;
        const auto gradient = dot(n, hessian * n) + factor;

        const auto incre = residual / gradient;

        auto error = fabs(residual);
        if(1u == counter) ref_error = std::max(1., error);
        suanpan_debug("Local plasticity iteration error: {:.5E}.\n", error /= ref_error);
        if(error <= tolerance || fabs(incre) <= tolerance) {
            const rowvec nts = (n.t() - four_third * k * dk * gamma / eta * (n % tensor::strain::norm_weight).t() * proj_a) * hessian;
            trial_stiffness = hessian * (eye(6, 6) - n * nts / (dot(nts, n) + factor));
            return SUANPAN_SUCCESS;
        }

        gamma += incre;
    }
}

int NonlinearHoffman::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int NonlinearHoffman::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int NonlinearHoffman::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void NonlinearHoffman::print() {
    suanpan_info("A 3D nonlinear hardening model using Hoffman yielding criterion.\n");
}
