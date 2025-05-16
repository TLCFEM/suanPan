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

#include "ArmstrongFrederick.h"

#include <Toolbox/tensor.h>

const double ArmstrongFrederick::root_three_two = sqrt(1.5);
const mat ArmstrongFrederick::unit_dev_tensor = tensor::unit_deviatoric_tensor4();

ArmstrongFrederick::ArmstrongFrederick(const unsigned T, DataArmstrongFrederick&& D, const double R)
    : DataArmstrongFrederick(std::move(D))
    , Material3D(T, R) {}

int ArmstrongFrederick::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = tensor::isotropic_stiffness(elastic_modulus, poissons_ratio);

    initialize_history(1 + 6 * size);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> ArmstrongFrederick::get_copy() { return std::make_unique<ArmstrongFrederick>(*this); }

double ArmstrongFrederick::get_parameter(const ParameterType P) const { return material_property(elastic_modulus, poissons_ratio)(P); }

int ArmstrongFrederick::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress = current_stress + (trial_stiffness = initial_stiffness) * incre_strain;

    trial_history = current_history;
    auto& p = trial_history(0);

    const auto trial_s = tensor::dev(trial_stress);

    auto eta = trial_s;
    for(unsigned I = 0; I < size; ++I) eta -= vec{&trial_history(1 + 6ull * I), 6, false, true};

    auto yield_func = root_three_two * tensor::stress::norm(eta) - std::max(0., yield + hardening * p + saturation * (1. - exp(-ms * p)));

    if(yield_func < 0.) return SUANPAN_SUCCESS;

    vec xi;
    auto gamma = 0.;
    double norm_xi, jacobian;

    auto counter = 0u;
    auto ref_error = 1.;
    while(true) {
        if(max_iteration == ++counter) {
            suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
            return SUANPAN_FAIL;
        }

        const auto exp_term = saturation * exp(-ms * p);

        auto k = yield + saturation + hardening * p - exp_term;
        auto dk = 0.;
        k < 0. ? k = 0. : dk = hardening + ms * exp_term;

        vec6 sum_a(fill::zeros);
        auto sum_b = 0.;
        for(unsigned I = 0; I < size; ++I) {
            const auto denom = 1. + b(I) * gamma;
            sum_a += vec{&trial_history(1 + 6ull * I), 6, false, true} / denom;
            sum_b += a(I) * gamma / denom;
        }

        norm_xi = tensor::stress::norm(xi = trial_s - sum_a);

        yield_func = root_three_two * (norm_xi - root_six_shear * gamma - sum_b) - k;

        sum_b = 0.;
        for(unsigned I = 0; I < size; ++I) sum_b += (b(I) / norm_xi * tensor::stress::double_contraction(xi, vec{&trial_history(1 + 6ull * I), 6, false, true}) - a(I)) * pow(1. + b(I) * gamma, -2.);

        jacobian = root_three_two * sum_b - three_shear - dk;

        const auto incre = yield_func / jacobian;
        const auto error = fabs(incre);
        if(1u == counter) ref_error = error;
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error < tolerance * ref_error || ((fabs(yield_func) < tolerance || error < tolerance) && counter > 5u)) break;

        gamma -= incre;
        p -= incre;
    }

    const vec u = xi / norm_xi;

    vec6 sum_c(fill::zeros);
    for(unsigned I = 0; I < size; ++I) {
        vec beta(&trial_history(1 + 6ull * I), 6, false, true);
        sum_c += b(I) * pow(1. + b(I) * gamma, -2.) * (beta - tensor::stress::double_contraction(u, beta) * u);
        beta = (beta + a(I) * gamma * u) / (1. + b(I) * gamma);
    }

    trial_stress -= root_six_shear * gamma * u;

    trial_stiffness += (root_six_shear * root_six_shear * gamma / jacobian / norm_xi * sum_c + root_six_shear * (root_six_shear / jacobian + double_shear * gamma / norm_xi) * u) * u.t() - double_shear * root_six_shear * gamma / norm_xi * unit_dev_tensor;

    return SUANPAN_SUCCESS;
}

int ArmstrongFrederick::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int ArmstrongFrederick::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int ArmstrongFrederick::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void ArmstrongFrederick::print() {
    suanpan_info("A 3D nonlinear hardening model using Armstrong-Frederick kinematic hardening rule.\n");
}
