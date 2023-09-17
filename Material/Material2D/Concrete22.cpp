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

#include "Concrete22.h"
#include <Domain/DomainBase.h>
#include <Toolbox/tensor.h>

Concrete22::Concrete22(const unsigned T, const double CS, const double TS, const double MCC, const double NCC, const double MTT, const double NTT, const double MP, const double CE, const double TE, const double SS, const double SR, const double R)
    : Material2D(T, PlaneType::N, R)
    , concrete_major(0, CS, TS, MCC, NCC, MTT, NTT, MP, CE, TE, R)
    , concrete_minor(0, CS, TS, MCC, NCC, MTT, NTT, MP, CE, TE, R)
    , shear_stress(SS)
    , shear_retention(SR) {}

int Concrete22::initialize(const shared_ptr<DomainBase>&) {
    if(SUANPAN_SUCCESS != concrete_major.initialize_base(nullptr) || SUANPAN_SUCCESS != concrete_major.initialize(nullptr) || SUANPAN_SUCCESS != concrete_minor.initialize_base(nullptr) || SUANPAN_SUCCESS != concrete_minor.initialize(nullptr)) return SUANPAN_FAIL;

    initial_stiffness.zeros(3, 3);
    initial_stiffness(2, 2) = shear_modulus = .5 * (initial_stiffness(0, 0) = initial_stiffness(1, 1) = concrete_major.get_parameter(ParameterType::ELASTICMODULUS));

    trial_stiffness = current_stiffness = initial_stiffness;

    shear_strain = shear_stress / shear_modulus;

    initialize_history(3);

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Concrete22::get_copy() { return make_unique<Concrete22>(*this); }

double Concrete22::get_parameter(const ParameterType P) const { return concrete_major.get_parameter(P); }

int Concrete22::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_history = current_history;
    auto& principal_angle = trial_history(0);
    auto& max_t_strain = trial_history(1);
    auto& max_c_strain = trial_history(2);

    bool if_rotating = false;

    if(max_t_strain < get_parameter(ParameterType::CRACKSTRAIN) && max_c_strain > get_parameter(ParameterType::PEAKSTRAIN)) {
        if_rotating = true;
        principal_angle = transform::strain::angle(trial_strain);
    }

    const auto trans = transform::strain::trans(principal_angle);

    const vec principal_strain = trans * trial_strain;
    const auto& strain_11 = principal_strain(0);
    const auto& strain_22 = principal_strain(1);

    if(strain_11 > max_t_strain) max_t_strain = strain_11;
    if(strain_22 < max_c_strain) max_c_strain = strain_22;

    // update status
    if(concrete_major.Material::update_trial_status(strain_11) + concrete_minor.Material::update_trial_status(strain_22) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

    // collect principal stiffness components
    trial_stiffness(0, 0) = concrete_major.get_trial_stiffness().at(0);
    trial_stiffness(1, 1) = concrete_minor.get_trial_stiffness().at(0);

    // collect principal stress components
    vec principal_stress(3);
    principal_stress(0) = concrete_major.get_trial_stress().at(0);
    principal_stress(1) = concrete_minor.get_trial_stress().at(0);

    if(if_rotating) {
        // doi.org/10.1061/(ASCE)0733-9399(1989)115:3(578)
        trial_stiffness(2, 2) = .5 * (principal_stress(0) - principal_stress(1)) / (principal_strain(0) - principal_strain(1)); // equation (9)
        if(!std::isfinite(trial_stiffness(2, 2))) trial_stiffness(2, 2) = shear_modulus;
        principal_stress(2) = 0.;
    }
    else if(principal_strain(2) > shear_strain) {
        trial_stiffness(2, 2) = shear_retention * shear_modulus;
        principal_stress(2) = shear_stress + (principal_strain(2) - shear_strain) * trial_stiffness(2, 2);
    }
    else if(principal_strain(2) < -shear_strain) {
        trial_stiffness(2, 2) = shear_retention * shear_modulus;
        principal_stress(2) = -shear_stress + (principal_strain(2) + shear_strain) * trial_stiffness(2, 2);
    }
    else {
        trial_stiffness(2, 2) = shear_modulus;
        principal_stress(2) = principal_strain(2) * shear_modulus;
    }

    // transform back to nominal direction
    trial_stress = trans.t() * principal_stress;
    trial_stiffness = trans.t() * diagmat(trial_stiffness) * trans;

    // BFGS type update
    // const vec elastic_stress = current_stiffness * trial_strain;
    // trial_stiffness = current_stiffness + trial_stress * trial_stress.t() / dot(trial_stress, trial_strain) - elastic_stress * elastic_stress.t() / dot(trial_strain, elastic_stress);

    suanpan_assert([&] { if(!trial_stress.is_finite() || !trial_stiffness.is_finite()) throw invalid_argument("infinite number detected"); });

    return SUANPAN_SUCCESS;
}

int Concrete22::clear_status() {
    trial_strain = current_strain.zeros();
    trial_stress = current_stress.zeros();
    trial_history = current_history = initial_history;
    trial_stiffness = current_stiffness = initial_stiffness;
    return concrete_major.clear_status() + concrete_minor.clear_status();
}

int Concrete22::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    current_history = trial_history;
    return concrete_major.commit_status() + concrete_minor.commit_status();
}

int Concrete22::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    trial_history = current_history;
    return concrete_major.reset_status() + concrete_minor.reset_status();
}

void Concrete22::print() {
    suanpan_info("Strain:", current_strain);
    suanpan_info("Stress:", current_stress);
}
