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

#include "Kelvin.h"

#include <Domain/DomainBase.h>
#include <Recorder/OutputType.h>

Kelvin::Kelvin(const unsigned T, const unsigned DT, const unsigned ST)
    : Material1D(T, 0.)
    , damper_tag(DT)
    , spring_tag(ST) {}

int Kelvin::initialize(const shared_ptr<DomainBase>& D) {
    damper = D->initialized_material_copy(damper_tag);
    spring = D->initialized_material_copy(spring_tag);

    if(nullptr == damper || nullptr == spring) return SUANPAN_FAIL;

    trial_strain_rate = current_strain_rate = incre_strain_rate.zeros(1);

    trial_damping = current_damping = initial_damping = damper->get_initial_damping();
    trial_stiffness = current_stiffness = initial_stiffness = spring->get_initial_stiffness();

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Kelvin::get_copy() { return std::make_unique<Kelvin>(*this); }

int Kelvin::update_trial_status(const vec&) {
    suanpan_error("Receives strain only from the associated element.\n");
    return SUANPAN_FAIL;
}

int Kelvin::update_trial_status(const vec& t_strain, const vec& t_strain_rate) {
    incre_strain = (trial_strain = t_strain) - current_strain;
    incre_strain_rate = (trial_strain_rate = t_strain_rate) - current_strain_rate;

    if(fabs(incre_strain(0) + fabs(incre_strain_rate(0))) <= datum::eps) return SUANPAN_SUCCESS;

    spring->update_trial_status(trial_strain);
    damper->update_trial_status(trial_strain, trial_strain_rate);

    trial_stiffness = spring->get_trial_stiffness();
    trial_damping = damper->get_trial_damping();

    if(!damper->get_trial_stiffness().empty()) trial_stiffness += damper->get_trial_stiffness();

    trial_stress = damper->get_trial_stress() + spring->get_trial_stress();

    return SUANPAN_SUCCESS;
}

int Kelvin::clear_status() {
    trial_strain = current_strain.zeros();
    trial_stress = current_stress.zeros();
    trial_strain_rate = current_strain_rate.zeros();
    trial_damping = current_damping = initial_damping;
    trial_stiffness = current_stiffness = initial_stiffness;
    return spring->clear_status() + damper->clear_status();
}

int Kelvin::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_strain_rate = trial_strain_rate;
    current_damping = trial_damping;
    current_stiffness = trial_stiffness;
    return spring->commit_status() + damper->commit_status();
}

int Kelvin::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_strain_rate = current_strain_rate;
    trial_damping = current_damping;
    trial_stiffness = current_stiffness;
    return spring->reset_status() + damper->reset_status();
}

std::vector<vec> Kelvin::record(const OutputType P) {
    if(OutputType::S == P) return {current_stress};
    if(OutputType::E == P) return {current_strain};
    if(OutputType::V == P) return {current_strain_rate};
    if(OutputType::SD == P) return {damper->get_current_stress()};
    if(OutputType::ED == P) return {damper->get_current_strain()};
    if(OutputType::VD == P) return {damper->get_current_strain_rate()};
    if(OutputType::SS == P) return {spring->get_current_stress()};
    if(OutputType::ES == P) return {spring->get_current_strain()};
    if(OutputType::VS == P) return {spring->get_current_strain_rate()};

    return {};
}

void Kelvin::print() {
    suanpan_info("A Kelvin material model.\n");
}
