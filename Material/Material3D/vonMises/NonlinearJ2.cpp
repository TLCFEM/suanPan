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

#include "NonlinearJ2.h"

#include <Toolbox/tensor.h>

const double NonlinearJ2::root_two_third = std::sqrt(two_third);
const mat NonlinearJ2::unit_dev_tensor = tensor::unit_deviatoric_tensor4();

NonlinearJ2::NonlinearJ2(const unsigned T, const double E, const double V, const double R)
    : DataNonlinearJ2{E, V}
    , Material3D(T, R) {}

int NonlinearJ2::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic_modulus, poissons_ratio);

    initialize_history(7);

    return SUANPAN_SUCCESS;
}

double NonlinearJ2::get_parameter(const ParameterType P) const { return material_property(elastic_modulus, poissons_ratio)(P); }

int NonlinearJ2::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    trial_history = current_history;
    auto& plastic_strain = trial_history(0);
    vec back_stress(&trial_history(1), 6, false, true);

    const vec rel_stress = tensor::dev(trial_stress) - back_stress;
    const auto norm_rel_stress = tensor::stress::norm(rel_stress);

    double k, dk;

    auto update_isotropic_hardening = [&] {
        k = compute_k(plastic_strain);
        dk = 0.;
        k < 0. ? k = 0. : dk = compute_dk(plastic_strain);
    };

    update_isotropic_hardening();

    auto yield_func = norm_rel_stress - root_two_third * k;

    if(yield_func <= 0.) return SUANPAN_SUCCESS;

    const auto current_h = compute_h(plastic_strain);
    auto gamma = 0., incre_h = 0.;
    double denom;
    auto counter = 0u;
    auto ref_error = 1.;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        denom = double_shear + two_third * (dk + compute_dh(plastic_strain));
        const auto incre_gamma = yield_func / denom;
        const auto error = std::fabs(incre_gamma);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || ((error < tolerance || std::fabs(yield_func) < tolerance) && counter > 5u)) break;
        incre_h = compute_h(plastic_strain = current_history(0) + root_two_third * (gamma += incre_gamma)) - current_h;
        update_isotropic_hardening();
        yield_func = norm_rel_stress - double_shear * gamma - root_two_third * (k + incre_h);
    }

    back_stress += root_two_third * incre_h / norm_rel_stress * rel_stress;

    auto t_factor = double_shear * gamma / norm_rel_stress;
    trial_stress -= t_factor * rel_stress;

    t_factor *= double_shear;
    trial_stiffness += (t_factor - square_double_shear / denom) / norm_rel_stress / norm_rel_stress * rel_stress * rel_stress.t() - t_factor * unit_dev_tensor;

    return SUANPAN_SUCCESS;
}

int NonlinearJ2::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int NonlinearJ2::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int NonlinearJ2::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void NonlinearJ2::print() {
    suanpan_info("A 3D nonlinear hardening model using von-Mises yielding criterion.\n");
}
