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

#include "StressWrapper.h"
#include <Domain/DomainBase.h>

mat StressWrapper::form_stiffness(const mat& full_stiffness) const { return full_stiffness(F1, F1) - full_stiffness(F1, F2) * solve(full_stiffness(F2, F2), full_stiffness(F2, F1)); }

StressWrapper::StressWrapper(const unsigned T, const unsigned BT, const unsigned MI, uvec&& FA, uvec&& FB, const MaterialType MT)
    : Material(T, MT, 0.)
    , F1(std::forward<uvec>(FA))
    , F2(std::forward<uvec>(FB))
    , base_tag(BT)
    , max_iteration(MI) {}

int StressWrapper::initialize(const shared_ptr<DomainBase>& D) {
    base = suanpan::initialized_material_copy(D, base_tag);

    if(nullptr == base || base->get_material_type() != MaterialType::D3) {
        suanpan_error("A valid 3D host material is required.\n");
        return SUANPAN_FAIL;
    }

    trial_full_strain = current_full_strain.zeros(6);

    trial_stiffness = current_stiffness = initial_stiffness = form_stiffness(base->get_initial_stiffness());

    return SUANPAN_SUCCESS;
}

double StressWrapper::get_parameter(const ParameterType P) const { return base->get_parameter(P); }

int StressWrapper::update_trial_status(const vec& t_strain) {
    auto& t_stress = base->get_trial_stress();
    auto& t_stiffness = base->get_trial_stiffness();

    if(norm(incre_strain = t_strain - trial_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_full_strain(F1) = trial_strain = t_strain;

    if(1 == max_iteration) {
        trial_full_strain(F2) -= solve(t_stiffness(F2, F2), t_stress(F2) + t_stiffness(F2, F1) * incre_strain);

        if(SUANPAN_SUCCESS != base->update_trial_status(trial_full_strain)) return SUANPAN_FAIL;
    }
    else {
        auto counter = 0u;
        auto ref_error = 1.;
        while(true) {
            // do not fail the analysis here
            // some material models may have large tolerance
            if(max_iteration == ++counter) break;

            if(SUANPAN_SUCCESS != base->update_trial_status(trial_full_strain)) return SUANPAN_FAIL;

            const auto error = inf_norm(t_stress(F2));
            if(1u == counter) ref_error = error;
            suanpan_debug("Local iteration error: {:.5E}.\n", error);
            if(error < tolerance * ref_error || error < datum::eps) break;

            trial_full_strain(F2) -= solve(t_stiffness(F2, F2), t_stress(F2));
        }
    }

    trial_stress = t_stress(F1) - t_stiffness(F1, F2) * solve(t_stiffness(F2, F2), t_stress(F2));

    trial_stiffness = form_stiffness(t_stiffness);

    return SUANPAN_SUCCESS;
}

int StressWrapper::clear_status() {
    trial_full_strain = current_full_strain.zeros(6);
    trial_strain = current_strain.zeros();
    trial_stress = current_stress.zeros();
    trial_stiffness = current_stiffness = initial_stiffness;
    return base->clear_status();
}

int StressWrapper::commit_status() {
    current_full_strain = trial_full_strain;
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    return base->commit_status();
}

int StressWrapper::reset_status() {
    trial_full_strain = current_full_strain;
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    return base->reset_status();
}

std::vector<vec> StressWrapper::record(const OutputType P) { return base->record(P); }
