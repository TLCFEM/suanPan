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

#include "UniversalOS.h"
#include <Domain/DomainBase.h>

UniversalOS::UniversalOS(const unsigned T, const unsigned BT, const unsigned MI, uvec&& FA, uvec&& FB)
    : StressWrapper(T, BT, MI, std::move(FA), std::move(FB), MaterialType::OS) {}

void UniversalOS::print() {
    suanpan_info("An open section material wrapper.\n");
    suanpan_info("Strain:", current_strain);
    suanpan_info("Stress:", current_stress);
    if(base) base->print();
}

OS146::OS146(const unsigned T, const unsigned BT, const unsigned MI)
    : UniversalOS(T, BT, MI, uvec{0, 3, 5}, uvec{1, 2, 4}) {}

unique_ptr<Material> OS146::get_copy() { return make_unique<OS146>(*this); }

OS146S::OS146S(const unsigned T, const unsigned BT, const double G)
    : Material(T, MaterialType::OS)
    , base_tag(BT)
    , shear_modulus(G) {}

int OS146S::initialize(const shared_ptr<DomainBase>& D) {
    base = D->initialized_material_copy(base_tag);

    if(nullptr == base || base->get_material_type() != MaterialType::D1) {
        suanpan_error("A valid 1D host material is required.\n");
        return SUANPAN_FAIL;
    }

    access::rw(density) = base->get_density();

    trial_stiffness = current_stiffness = initial_stiffness = diagmat(vec{base->get_initial_stiffness()(0), shear_modulus, shear_modulus});

    return SUANPAN_SUCCESS;
}

double OS146S::get_parameter(const ParameterType P) const { return base->get_parameter(P); }

unique_ptr<Material> OS146S::get_copy() { return make_unique<OS146S>(*this); }

int OS146S::update_trial_status(const vec& t_strain) {
    trial_strain = t_strain;
    if(SUANPAN_SUCCESS != base->update_trial_status(t_strain(0))) return SUANPAN_FAIL;

    trial_stress = shear_modulus * trial_strain;
    trial_stress(0) = base->get_trial_stress()(0);
    trial_stiffness = diagmat(vec{base->get_trial_stiffness()(0), shear_modulus, shear_modulus});

    return SUANPAN_SUCCESS;
}

int OS146S::clear_status() {
    trial_strain = current_strain.zeros();
    trial_stress = current_stress.zeros();
    trial_stiffness = current_stiffness = initial_stiffness;
    return base->clear_status();
}

int OS146S::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    return base->commit_status();
}

int OS146S::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    return base->reset_status();
}

std::vector<vec> OS146S::record(const OutputType P) { return base->record(P); }

void OS146S::print() { if(base) base->print(); }
