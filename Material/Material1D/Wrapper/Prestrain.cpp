﻿/*******************************************************************************
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

#include "Prestrain.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Load/Amplitude/Ramp.h>

double Prestrain::get_prestrain() const { return amplitude->get_amplitude(*analysis_time) * magnitude; }

Prestrain::Prestrain(const unsigned T, const unsigned BT, const unsigned AT, const double M)
    : Material1D(T, 0.)
    , base_tag(BT)
    , amplitude_tag(AT)
    , magnitude(M) {}

int Prestrain::initialize(const shared_ptr<DomainBase>& D) {
    amplitude = D->get<Amplitude>(amplitude_tag);
    if(nullptr == amplitude || !amplitude->is_active()) {
        if(0u != amplitude_tag) suanpan_warning("The provided amplitude {} is not usable, using a default one instead.\n", amplitude_tag);
        amplitude = Ramp(0);
    }
    base = D->initialized_material_copy(base_tag);
    if(nullptr == base || base->get_material_type() != MaterialType::D1) {
        suanpan_error("A valid uniaxial host material is required.\n");
        return SUANPAN_FAIL;
    }

    analysis_time = &D->get_factory()->modify_trial_time();

    access::rw(density) = base->get_density();

    trial_stiffness = current_stiffness = base->get_initial_stiffness();
    trial_damping = current_damping = base->get_initial_damping();

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Prestrain::get_copy() { return std::make_unique<Prestrain>(*this); }

int Prestrain::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain + get_prestrain()) - current_strain;

    if(base->update_trial_status(trial_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

    trial_stress = base->get_trial_stress();
    trial_stiffness = base->get_trial_stiffness();

    return SUANPAN_SUCCESS;
}

int Prestrain::update_trial_status(const vec& t_strain, const vec& t_strain_rate) {
    incre_strain = (trial_strain = t_strain + get_prestrain()) - current_strain;
    incre_strain_rate = (trial_strain_rate = t_strain_rate) - current_strain_rate;

    if(base->update_trial_status(trial_strain, trial_strain_rate) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

    trial_stress = base->get_trial_stress();
    trial_stiffness = base->get_trial_stiffness();
    trial_damping = base->get_trial_damping();

    return SUANPAN_SUCCESS;
}

int Prestrain::clear_status() {
    trial_strain = current_strain.zeros();
    trial_strain_rate = current_strain_rate.zeros();
    trial_stress = current_stress.zeros();
    trial_stiffness = current_stiffness = initial_stiffness;
    trial_damping = current_damping = initial_damping;
    return base->clear_status();
}

int Prestrain::commit_status() {
    current_strain = trial_strain;
    current_strain_rate = trial_strain_rate;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    current_damping = trial_damping;
    return base->commit_status();
}

int Prestrain::reset_status() {
    trial_strain = current_strain;
    trial_strain_rate = current_strain_rate;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    trial_damping = current_damping;
    return base->reset_status();
}

std::vector<vec> Prestrain::record(const OutputType P) { return base->record(P); }

void Prestrain::print() {
    suanpan_info("A Prestrain container that holds the following material.\n");
    base->print();
}
