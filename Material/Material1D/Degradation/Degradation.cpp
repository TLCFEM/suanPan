/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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

#include "Degradation.h"
#include <Domain/DomainBase.h>
#include <Recorder/OutputType.h>

Degradation::Degradation(const unsigned T, const unsigned MT)
    : Material1D(T, 0.)
    , mat_tag(MT) {}

Degradation::Degradation(const Degradation& old_obj)
    : Material1D(old_obj)
    , mat_tag(old_obj.mat_tag)
    , base(suanpan::make_copy(old_obj.base)) {}

int Degradation::initialize(const shared_ptr<DomainBase>& D) {
    base = D->initialized_material_copy(mat_tag);

    if(nullptr == base) return SUANPAN_FAIL;

    access::rw(density) = base->get_density();

    const auto degrade = compute_positive_degradation(0.);
    const auto& d = degrade(0);

    trial_stiffness = current_stiffness = initial_stiffness = d * base->get_initial_stiffness();

    return SUANPAN_SUCCESS;
}

int Degradation::clear_status() {
    trial_strain = current_strain.zeros();
    trial_stress = current_stress.zeros();
    trial_history = current_history = initial_history;
    trial_stiffness = current_stiffness = initial_stiffness;
    return base->clear_status();
}

int Degradation::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_history = trial_history;
    current_stiffness = trial_stiffness;
    return base->commit_status();
}

int Degradation::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_history = current_history;
    trial_stiffness = current_stiffness;
    return base->reset_status();
}

int StrainDegradation::initialize(const shared_ptr<DomainBase>& D) {
    initialize_history(2);

    return Degradation::initialize(D);
}

int StrainDegradation::update_trial_status(const vec& t_strain) {
    if(SUANPAN_SUCCESS != base->update_trial_status(trial_strain = t_strain)) return SUANPAN_FAIL;

    trial_stress = base->get_trial_stress();

    trial_history = current_history;
    auto& max_strain = trial_history(0);
    auto& min_strain = trial_history(1);

    if(trial_stress(0) > 0.)
        // tension
        if(trial_strain(0) > max_strain) {
            max_strain = trial_strain(0);
            const auto degradation = compute_positive_degradation(max_strain);
            const auto& d = degradation(0);
            const auto& dd = degradation(1);

            trial_stress = d * base->get_trial_stress();

            trial_stiffness = d * base->get_trial_stiffness() + dd * base->get_trial_stress();
        }
        else {
            const auto degrade = compute_positive_degradation(max_strain);
            const auto& d = degrade(0);

            trial_stress = d * base->get_trial_stress();

            trial_stiffness = d * base->get_trial_stiffness();
        }
    else
        // compression
        if(trial_strain(0) < min_strain) {
            min_strain = trial_strain(0);
            const auto degradation = compute_negative_degradation(min_strain);
            const auto& d = degradation(0);
            const auto& dd = degradation(1);

            trial_stress = d * base->get_trial_stress();

            trial_stiffness = d * base->get_trial_stiffness() + dd * base->get_trial_stress();
        }
        else {
            const auto degrade = compute_negative_degradation(min_strain);
            const auto& d = degrade(0);

            trial_stress = d * base->get_trial_stress();

            trial_stiffness = d * base->get_trial_stiffness();
        }

    return SUANPAN_SUCCESS;
}

vector<vec> StrainDegradation::record(const OutputType P) {
    if(OutputType::DT == P) return {vec{compute_positive_degradation(current_history(0))(0)}};
    if(OutputType::DC == P) return {vec{compute_negative_degradation(current_history(1))(0)}};

    return Degradation::record(P);
}

int StressDegradation::initialize(const shared_ptr<DomainBase>& D) {
    initialize_history(2);

    return Degradation::initialize(D);
}

int StressDegradation::update_trial_status(const vec& t_strain) {
    if(SUANPAN_SUCCESS != base->update_trial_status(trial_strain = t_strain)) return SUANPAN_FAIL;

    trial_stress = base->get_trial_stress();

    trial_history = current_history;
    auto& max_stress = trial_history(0);
    auto& min_stress = trial_history(1);

    if(trial_stress(0) > 0.)
        // tension
        if(trial_stress(0) > max_stress) {
            max_stress = trial_stress(0);
            const auto degradation = compute_positive_degradation(max_stress);
            const auto& d = degradation(0);
            const auto& dd = degradation(1);

            trial_stress = d * base->get_trial_stress();

            trial_stiffness = (d + dd * base->get_trial_stress()) * base->get_trial_stiffness();
        }
        else {
            const auto degradation = compute_positive_degradation(max_stress);
            const auto& d = degradation(0);

            trial_stress = d * base->get_trial_stress();

            trial_stiffness = d * base->get_trial_stiffness();
        }
    else
        // compression
        if(trial_stress(0) < min_stress) {
            min_stress = trial_stress(0);
            const auto degradation = compute_negative_degradation(min_stress);
            const auto& d = degradation(0);
            const auto& dd = degradation(1);

            trial_stress = d * base->get_trial_stress();

            trial_stiffness = (d + dd * base->get_trial_stress()) * base->get_trial_stiffness();
        }
        else {
            const auto degradation = compute_negative_degradation(min_stress);
            const auto& d = degradation(0);

            trial_stress = d * base->get_trial_stress();

            trial_stiffness = d * base->get_trial_stiffness();
        }

    return SUANPAN_SUCCESS;
}

vector<vec> StressDegradation::record(const OutputType P) {
    if(OutputType::DT == P) return {vec{compute_positive_degradation(current_history(0))(0)}};
    if(OutputType::DC == P) return {vec{compute_negative_degradation(current_history(1))(0)}};

    return Degradation::record(P);
}
