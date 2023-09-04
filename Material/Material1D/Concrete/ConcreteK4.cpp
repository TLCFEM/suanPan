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

void ConcreteK4::compute_tension_branch() {
    auto& plastic_strain = trial_history(0);
    auto& kt = trial_history(1);
    auto& kc = trial_history(2);
    auto& kk = trial_history(3);

    const auto sigma_t = f_t + hardening_t * kt;

    const auto f_trial = trial_stress(0) - sigma_t;

    if(f_trial <= 0.) return;

    const auto gamma = f_trial / (elastic_modulus + hardening_t);
    plastic_strain += gamma;
    kt += gamma;
    trial_stress -= elastic_modulus * gamma;
    trial_stiffness -= elastic_modulus / (elastic_modulus + hardening_t) * elastic_modulus;
}

void ConcreteK4::compute_compression_branch() {
    auto& plastic_strain = trial_history(0);
    auto& kt = trial_history(1);
    auto& kc = trial_history(2);
    auto& kk = trial_history(3);

    const bool first_segment = kc <= k_peak ? trial_stress(0) + k_peak * (elastic_modulus + hardening_c) - kc * elastic_modulus + f_y > 0. : false;

    double h;
    auto sigma_c = f_y;
    if(first_segment) {
        h = hardening_c;
        sigma_c += hardening_c * kc;
    }
    else {
        h = hardening_d;
        sigma_c += hardening_c * k_peak + (kc - k_peak) * hardening_d;
    }

    const auto f_trial = -trial_stress(0) - sigma_c;

    if(f_trial <= 0.) return;

    const auto gamma = f_trial / (elastic_modulus + h);
    plastic_strain -= gamma;
    kc += gamma;
    trial_stress += elastic_modulus * gamma;
    trial_stiffness -= elastic_modulus / (elastic_modulus + h) * elastic_modulus;
}

void ConcreteK4::compute_crack_close_branch() {
    auto& plastic_strain = trial_history(0);
    auto& kt = trial_history(1);
    auto& kc = trial_history(2);
    auto& kk = trial_history(3);

    const auto f_trial = fabs(trial_stress(0)) - hardening_k * kk;

    if(f_trial <= 0.) return;

    if(auto gamma = f_trial / (elastic_modulus + hardening_k); gamma < kt - kk) {
        plastic_strain -= gamma;
        kk += gamma;
        trial_stress += elastic_modulus * gamma;
        trial_stiffness -= elastic_modulus / (elastic_modulus + hardening_k) * elastic_modulus;
    }
    else {
        gamma = kt - kk;
        plastic_strain -= gamma;
        kk += gamma;
        trial_stress += elastic_modulus * gamma;
    }
}

ConcreteK4::ConcreteK4(const unsigned T, const double E, const double R)
    : DataConcreteK4{E, .1 * E, .3 * E, .1 * E, .1 * E}
    , Material1D(T, R) {}

int ConcreteK4::initialize(const shared_ptr<DomainBase>&) {
    trial_stiffness = current_stiffness = initial_stiffness = elastic_modulus;

    initialize_history(4);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> ConcreteK4::get_copy() { return make_unique<ConcreteK4>(*this); }

double ConcreteK4::get_parameter(const ParameterType P) const {
    if(ParameterType::DENSITY == P) return density;
    if(ParameterType::ELASTICMODULUS == P || ParameterType::YOUNGSMODULUS == P || ParameterType::E == P) return initial_stiffness(0);
    return 0.;
}

int ConcreteK4::update_trial_status(const vec& n_strain) {
    incre_strain = (trial_strain = n_strain) - current_strain;

    if(fabs(incre_strain(0)) <= 1E-15) return SUANPAN_SUCCESS;

    trial_history = current_history;
    const auto& plastic_strain = trial_history(0);
    const auto& current_kt = current_history(1);
    const auto& current_kc = current_history(2);
    const auto& current_kk = current_history(3);

    trial_stress = elastic_modulus * (trial_strain - plastic_strain);

    if(trial_stress(0) < 0. && incre_strain(0) < 0. && current_kk < current_kt) compute_crack_close_branch();

    trial_stress(0) > 0. ? compute_tension_branch() : compute_compression_branch();

    return SUANPAN_SUCCESS;
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
