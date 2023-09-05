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

#include "ConcreteK4.h"
#include <Toolbox/utility.h>

int ConcreteK4::compute_tension_branch() {
    auto& plastic_strain = trial_history(0);
    auto& kt = trial_history(1);

    const auto sign_sigma = suanpan::sign(trial_stress(0));

    double jacobian;

    auto counter = 0u;
    while(true) {
        if(max_iteration == ++counter) return SUANPAN_FAIL;

        const auto backbone = compute_tension_backbone(kt);
        const auto residual = fabs(trial_stress(0)) - backbone(0);

        if(1u == counter && residual <= 0.) {
            const auto damage = compute_tension_damage(kt);

            trial_stress *= 1. - damage(0);
            trial_stiffness *= 1. - damage(0);

            return SUANPAN_SUCCESS;
        }

        jacobian = elastic_modulus + backbone(1);
        const auto incre = residual / jacobian;
        const auto error = fabs(incre);
        suanpan_debug("Local tension iteration error: {:.5E}.\n", error);
        if(error < tolerance || fabs(residual) < tolerance) break;

        const auto incre_e = incre * sign_sigma;

        kt += incre;
        plastic_strain += incre_e;
        trial_stress -= elastic_modulus * incre_e;
    }

    trial_stiffness -= elastic_modulus / jacobian * elastic_modulus;

    const auto damage = compute_tension_damage(kt);

    trial_stiffness *= 1. - damage(0);
    trial_stiffness -= trial_stress * damage(1) * elastic_modulus / jacobian;

    trial_stress *= 1. - damage(0);

    return SUANPAN_SUCCESS;
}

int ConcreteK4::compute_compression_branch() {
    auto& plastic_strain = trial_history(0);
    auto& kc = trial_history(2);

    const auto sign_sigma = suanpan::sign(trial_stress(0));

    double jacobian;

    auto counter = 0u;
    while(true) {
        if(max_iteration == ++counter) return SUANPAN_FAIL;

        const auto backbone = compute_compression_backbone(kc);
        const auto residual = fabs(trial_stress(0)) - backbone(0);

        if(1u == counter && residual <= 0.) {
            const auto damage = compute_compression_damage(kc);

            trial_stress *= 1. - damage(0);
            trial_stiffness *= 1. - damage(0);

            return SUANPAN_SUCCESS;
        }

        jacobian = elastic_modulus + backbone(1);
        const auto incre = residual / jacobian;
        const auto error = fabs(incre);
        suanpan_debug("Local compression iteration error: {:.5E}.\n", error);
        if(error < tolerance || fabs(residual) < tolerance) break;

        const auto incre_e = incre * sign_sigma;

        kc += incre;
        plastic_strain += incre_e;
        trial_stress -= elastic_modulus * incre_e;
    }

    trial_stiffness -= elastic_modulus / jacobian * elastic_modulus;

    const auto damage = compute_compression_damage(kc);

    trial_stiffness *= 1. - damage(0);
    trial_stiffness -= trial_stress * damage(1) * elastic_modulus / jacobian;

    trial_stress *= 1. - damage(0);

    return SUANPAN_SUCCESS;
}

int ConcreteK4::compute_crack_close_branch() {
    auto& plastic_strain = trial_history(0);
    const auto& kt = trial_history(1);
    auto& kk = trial_history(3);

    const auto jacobian = elastic_modulus + hardening_k;
    auto incre = fabs(incre_strain(0)) * elastic_modulus / jacobian;

    if(incre > kt - kk) incre = kt - kk;
    else trial_stiffness -= elastic_modulus / jacobian * elastic_modulus;

    const auto incre_e = incre * suanpan::sign(trial_stress(0));

    kk += incre;
    plastic_strain += incre_e;
    trial_stress -= elastic_modulus * incre_e;

    return SUANPAN_SUCCESS;
}

ConcreteK4::ConcreteK4(const unsigned T, const double E, const double R)
    : DataConcreteK4{E, .2 * E}
    , Material1D(T, R) {}

int ConcreteK4::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    initialize_history(4);

    return SUANPAN_SUCCESS;
}

double ConcreteK4::get_parameter(const ParameterType P) const {
    if(ParameterType::DENSITY == P) return density;
    if(ParameterType::ELASTICMODULUS == P || ParameterType::YOUNGSMODULUS == P || ParameterType::E == P) return initial_stiffness(0);
    return 0.;
}

int ConcreteK4::update_trial_status(const vec& n_strain) {
    incre_strain = (trial_strain = n_strain) - current_strain;

    if(fabs(incre_strain(0)) <= tolerance) return SUANPAN_SUCCESS;

    trial_history = current_history;
    const auto& plastic_strain = trial_history(0);
    const auto& current_kt = current_history(1);
    const auto& current_kk = current_history(3);

    trial_stress = elastic_modulus * (trial_strain - plastic_strain);

    if(trial_stress(0) < 0. && incre_strain(0) < 0. && current_kt > current_kk && SUANPAN_SUCCESS != compute_crack_close_branch()) return SUANPAN_FAIL;

    return trial_stress(0) > 0. ? compute_tension_branch() : compute_compression_branch();
}

int ConcreteK4::clear_status() {
    current_strain.zeros();
    current_stress.zeros();
    current_history = initial_history;
    current_stiffness = initial_stiffness;
    return reset_status();
}

int ConcreteK4::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return SUANPAN_SUCCESS;
}

int ConcreteK4::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return SUANPAN_SUCCESS;
}

void ConcreteK4::print() {
    suanpan_info("A concrete model. doi: 10.1061/(ASCE)ST.1943-541X.000259\n");
    Material1D::print();
}

vec2 LinearK4::compute_tension_backbone(const double k) const { return vec2{f_t + hardening_t * k, hardening_t}; }

vec2 LinearK4::compute_compression_backbone(const double k) const {
    if(k < k_peak) return vec2{f_y + hardening_c * k, hardening_c};

    return vec2{f_c + hardening_d * (k - k_peak), hardening_d};
}

vec2 LinearK4::compute_tension_damage(const double k) const {
    const auto factor = exp(-k / e_ft);
    return vec2{1. - factor, factor / e_ft};
}

vec2 LinearK4::compute_compression_damage(double k) const {
    if(k < k_peak) return vec2{0., 0.};

    k -= k_peak;

    const auto factor = exp(-k / e_fc);
    return vec2{1. - factor, factor / e_fc};
}

LinearK4::LinearK4(const unsigned T, const double E, const double R)
    : DataLinearK4{.15 * E, .6 * E}
    , ConcreteK4(T, E, R) {}

unique_ptr<Material> LinearK4::get_copy() { return make_unique<LinearK4>(*this); }
