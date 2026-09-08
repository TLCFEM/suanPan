/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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
#include <set>

int StressWrapper::form_stiffness(mat& condensed_mat, const mat& full_stiffness) const {
    if(mat aux_mat; solve(aux_mat, full_stiffness(F2, F2), full_stiffness(F2, F1))) {
        condensed_mat = full_stiffness(F1, F1) - full_stiffness(F1, F2) * aux_mat;
        return SUANPAN_SUCCESS;
    }

    return SUANPAN_FAIL;
}

StressWrapper::StressWrapper(const unsigned T, const unsigned BT, const unsigned MI, uvec&& FA, const MaterialType MT)
    : Material(T, MT, 0.)
    , F1(std::move(FA))
    , base_tag(BT)
    , max_iteration(MI) {}

int StressWrapper::initialize_base(const shared_ptr<DomainBase>& D) {
    base = D->initialized_material_copy(base_tag);

    if(nullptr == base || base->get_material_type() != MaterialType::D3) {
        suanpan_error("A valid 3D host material is required.\n");
        return SUANPAN_FAIL;
    }

    return Material::initialize_base(D);
}

int StressWrapper::initialize(const shared_ptr<DomainBase>&) {
    access::rw(density) = base->get_density();

    const std::set nontrivial(F1.begin(), F1.end());

    std::vector workspace(F1.begin(), F1.end());
    const auto total_size = base->nonlocal_size() + 6u;
    for(auto I = 6u; I < total_size; ++I) workspace.emplace_back(I);
    access::rw(F1) = workspace;

    workspace.clear();
    for(auto I = 0u; I < 6u; ++I)
        if(!nontrivial.contains(I)) workspace.emplace_back(I);
    access::rw(F2) = workspace;

    trial_full_strain = current_full_strain.zeros(total_size);

    if(form_stiffness(initial_stiffness, base->get_initial_stiffness()) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

unsigned StressWrapper::nonlocal_size() const { return base->nonlocal_size(); }

double StressWrapper::get(const Parameter P) const { return base->get(P); }

int StressWrapper::update_trial_status(const vec& t_strain) {
    auto& t_stress = base->get_trial_stress();
    auto& t_stiffness = base->get_trial_stiffness();

    if(norm(incre_strain = t_strain - trial_strain) <= datum::eps) return SUANPAN_SUCCESS;

    trial_full_strain(F1) = trial_strain = t_strain;

    vec t_incre;

    if(1u == max_iteration) {
        if(!solve(t_incre, t_stiffness(F2, F2), t_stress(F2) + t_stiffness(F2, F1) * incre_strain)) return SUANPAN_FAIL;

        trial_full_strain(F2) -= t_incre;

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

            if(!solve(t_incre, t_stiffness(F2, F2), t_stress(F2))) return SUANPAN_FAIL;

            const auto error = suanpan::inf_norm(t_incre);

            if(1u == counter) ref_error = error;
            suanpan_debug("Local iteration error: {:.5E}.\n", error);
            if(error < tolerance * ref_error || (suanpan::inf_norm(t_stress(F2)) < tolerance && counter > 5u)) break;

            trial_full_strain(F2) -= t_incre;
        }
    }

    if(!solve(t_incre, t_stiffness(F2, F2), t_stress(F2))) return SUANPAN_FAIL;

    trial_stress = t_stress(F1) - t_stiffness(F1, F2) * t_incre;

    return form_stiffness(trial_stiffness, t_stiffness);
}

int StressWrapper::clear_status() {
    trial_full_strain = current_full_strain.zeros();
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

std::vector<vec> StressWrapper::record(const OutputType P) const { return base->record(P); }
