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
        amplitude = std::make_shared<Ramp>(0);
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
    incre_strain = (trial_strain = t_strain) - current_strain;

    const auto prestrain = get_prestrain();

    if(std::fabs(incre_strain(0) + prestrain) <= datum::eps) return SUANPAN_SUCCESS;

    if(base->update_trial_status(trial_strain + prestrain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

    trial_stress = base->get_trial_stress();
    trial_stiffness = base->get_trial_stiffness();

    return SUANPAN_SUCCESS;
}

int Prestrain::update_trial_status(const vec& t_strain, const vec& t_strain_rate) {
    incre_strain = (trial_strain = t_strain) - current_strain;
    incre_strain_rate = (trial_strain_rate = t_strain_rate) - current_strain_rate;

    const auto prestrain = get_prestrain();

    if(std::fabs(incre_strain(0) + incre_strain_rate(0) + prestrain) <= datum::eps) return SUANPAN_SUCCESS;

    if(base->update_trial_status(trial_strain + prestrain, trial_strain_rate) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

    trial_stress = base->get_trial_stress();
    trial_stiffness = base->get_trial_stiffness();
    trial_damping = base->get_trial_damping();

    return SUANPAN_SUCCESS;
}

int Prestrain::clear_status() { return base->clear_status(); }

int Prestrain::commit_status() { return base->commit_status(); }

int Prestrain::reset_status() { return base->reset_status(); }

std::vector<vec> Prestrain::record(const OutputType P) { return base->record(P); }

void Prestrain::print() {
    suanpan_info("A Prestrain container that holds the following material.\n");
    base->print();
}
