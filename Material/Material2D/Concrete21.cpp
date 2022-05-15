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

#include "Concrete21.h"
#include <Domain/DomainBase.h>
#include <Toolbox/tensorToolbox.h>

Concrete21::Concrete21(const unsigned T, const double CS, const double TS, const double MCC, const double NCC, const double MTT, const double NTT, const double MP, const double CE, const double TE, const double R)
    : Material2D(T, PlaneType::N, R)
    , concrete_major(0, CS, TS, MCC, NCC, MTT, NTT, MP, CE, TE, R)
    , concrete_minor(0, CS, TS, MCC, NCC, MTT, NTT, MP, CE, TE, R) {}

int Concrete21::initialize(const shared_ptr<DomainBase>&) {
    if(SUANPAN_SUCCESS != concrete_major.initialize_base(nullptr) || SUANPAN_SUCCESS != concrete_major.initialize(nullptr) || SUANPAN_SUCCESS != concrete_minor.initialize_base(nullptr) || SUANPAN_SUCCESS != concrete_minor.initialize(nullptr)) return SUANPAN_FAIL;

    initial_stiffness.zeros(3, 3);
    initial_stiffness(2, 2) = shear_modulus = .5 * (initial_stiffness(0, 0) = initial_stiffness(1, 1) = concrete_major.get_parameter(ParameterType::ELASTICMODULUS));

    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Concrete21::get_copy() { return make_unique<Concrete21>(*this); }

double Concrete21::get_parameter(const ParameterType P) const {
    if(ParameterType::PLANETYPE == P) return static_cast<double>(plane_type);
    return concrete_major.get_parameter(P);
}

int Concrete21::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

    const auto principal_strain = transform::strain::principal(trial_strain);

    // update status
    concrete_major.Material::update_trial_status(principal_strain(0));
    concrete_minor.Material::update_trial_status(principal_strain(1));

    vec principal_stress(3);

    // collect principal stress components
    principal_stress(0) = concrete_major.get_trial_stress().at(0);
    principal_stress(1) = concrete_minor.get_trial_stress().at(0);
    principal_stress(2) = 0.;

    // collect principal stiffness components
    trial_stiffness(0, 0) = concrete_major.get_trial_stiffness().at(0);
    trial_stiffness(1, 1) = concrete_minor.get_trial_stiffness().at(0);

    // doi.org/10.1061/(ASCE)0733-9399(1989)115:3(578)
    trial_stiffness(2, 2) = .5 * (principal_stress(0) - principal_stress(1)) / (principal_strain(0) - principal_strain(1)); // equation (9)

    if(!std::isfinite(trial_stiffness(2, 2))) trial_stiffness(2, 2) = shear_modulus;

    const auto trans = transform::strain::trans(transform::strain::angle(trial_strain));

    // transform back to nominal direction
    trial_stress = trans.t() * principal_stress;
    trial_stiffness = trans.t() * diagmat(trial_stiffness) * trans;

    return SUANPAN_SUCCESS;
}

int Concrete21::clear_status() {
    current_strain = trial_strain.zeros();
    current_stress = trial_stress.zeros();
    trial_stiffness = current_stiffness = initial_stiffness;
    return concrete_major.clear_status() + concrete_minor.clear_status();
}

int Concrete21::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    return concrete_major.commit_status() + concrete_minor.commit_status();
}

int Concrete21::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    return concrete_major.reset_status() + concrete_minor.reset_status();
}

void Concrete21::print() {
    suanpan_info("A planar concrete model: \n");
    current_strain.t().print("Strain:");
    current_stress.t().print("Stress:");
}
