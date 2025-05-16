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

#include "NonlinearCamClay.h"

#include <Toolbox/tensor.h>

const double NonlinearCamClay::sqrt_three_two = sqrt(1.5);
const mat NonlinearCamClay::unit_dev_tensor = tensor::unit_deviatoric_tensor4();

NonlinearCamClay::NonlinearCamClay(const unsigned T, const double E, const double V, const double B, const double M, const double P, const double R)
    : DataNonlinearCamClay{fabs(E), fabs(V), B * B, fabs(M), fabs(P)}
    , Material3D(T, R) { access::rw(tolerance) = 1E-13; }

int NonlinearCamClay::initialize(const shared_ptr<DomainBase>&) {
    current_stiffness = trial_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic_modulus, poissons_ratio);

    initialize_history(1);

    return SUANPAN_SUCCESS;
}

double NonlinearCamClay::get_parameter(const ParameterType P) const { return material_property(elastic_modulus, poissons_ratio)(P); }

int NonlinearCamClay::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    trial_history = current_history;
    auto& alpha = trial_history(0);
    const auto& current_alpha = current_history(0);

    auto trial_s = tensor::dev(trial_stress);
    const auto trial_q = sqrt_three_two * tensor::stress::norm(trial_s);
    const auto p = tensor::mean3(trial_stress);

    auto ini_f = 0., gamma = 0.;

    vec2 residual(fill::none), incre(fill::none);
    mat22 jacobian(fill::none);

    auto counter = 0u;
    auto rel_error = 1.;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        const auto a = compute_a(alpha);
        const auto da = compute_da(alpha);
        const auto incre_alpha = alpha - current_alpha;
        const auto trial_p = p - bulk * incre_alpha;
        const auto rel_p = trial_p - pt + a;
        const auto square_b = rel_p >= 0. ? 1. : square_beta;
        const auto denom = square_m + six_shear * gamma;
        const auto square_qm = pow(m * trial_q / denom, 2.);

        residual(0) = rel_p * rel_p / square_b + square_qm - a * a;

        if(1u == counter) {
            if(residual(0) < 0.) return SUANPAN_SUCCESS;
            ini_f = std::max(1., residual(0)); // yield function can be very large, use relative error instead
        }

        residual(1) = incre_alpha - 2. * gamma / square_b * rel_p;

        jacobian(0, 0) = -2. * six_shear / denom * square_qm;
        jacobian(1, 0) = -2. * rel_p / square_b;
        jacobian(0, 1) = jacobian(1, 0) * (bulk - da) - 2. * a * da;
        jacobian(1, 1) = 1. - 2. * gamma / square_b * (da - bulk);

        if(!solve(incre, jacobian, residual, solve_opts::equilibrate)) return SUANPAN_FAIL;

        const auto error = inf_norm(incre);
        if(1u == counter) rel_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * rel_error || (inf_norm(residual) < tolerance * ini_f && counter > 5u)) {
            mat::fixed<6, 2> left(fill::none);

            rel_error = 2. * bulk / square_b; // reuse variable
            left.col(0) = rel_error * rel_p * tensor::unit_tensor2 + six_shear / denom * (trial_s *= square_m / denom);
            left.col(1) = -rel_error * gamma * tensor::unit_tensor2;

            trial_stress = trial_s + trial_p * tensor::unit_tensor2;

            rel_error = six_shear / denom; // reuse variable
            trial_stiffness -= 2. * shear * gamma * rel_error * unit_dev_tensor - join_rows(rel_error * trial_s, bulk * tensor::unit_tensor2) * solve(jacobian, left.t());

            return SUANPAN_SUCCESS;
        }

        gamma -= incre(0);
        alpha -= incre(1);
    }
}

int NonlinearCamClay::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int NonlinearCamClay::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int NonlinearCamClay::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}
