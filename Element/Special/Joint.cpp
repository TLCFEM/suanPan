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

#include "Joint.h"
#include <Domain/DomainBase.h>
#include <Material/Material.h>

Joint::Joint(const unsigned T, uvec&& NT, uvec&& MT)
    : MaterialElement1D(T, j_node, static_cast<unsigned>(MT.n_elem), std::forward<uvec>(NT), std::forward<uvec>(MT), false, {})
    , j_dof(Element::get_dof_number()) {}

int Joint::initialize(const shared_ptr<DomainBase>& D) {
    if(j_dof != material_tag.n_elem) return SUANPAN_FAIL;

    j_material.clear();
    j_material.reserve(j_dof);
    for(auto& I : material_tag) j_material.emplace_back(D->get<Material>(I)->get_copy());

    initial_stiffness.zeros(j_size, j_size);
    for(size_t I = 0, J = j_dof; I < j_dof; ++I, ++J) initial_stiffness(I, J) = initial_stiffness(J, I) = -(initial_stiffness(I, I) = initial_stiffness(J, J) = as_scalar(j_material[I]->get_initial_stiffness()));

    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}

int Joint::update_status() {
    const auto t_disp = get_trial_displacement();

    trial_resistance.zeros(j_size);
    trial_stiffness.zeros(j_size, j_size);
    for(unsigned I = 0, J = j_dof; I < j_material.size(); ++I, ++J) {
        if(SUANPAN_SUCCESS != j_material[I]->update_trial_status(t_disp(I) - t_disp(J))) return SUANPAN_FAIL;
        trial_resistance(J) = -(trial_resistance(I) = as_scalar(j_material[I]->get_trial_stress()));
        trial_stiffness(I, J) = trial_stiffness(J, I) = -(trial_stiffness(J, J) = trial_stiffness(I, I) = as_scalar(j_material[I]->get_trial_stiffness()));
    }

    return SUANPAN_SUCCESS;
}

int Joint::commit_status() {
    auto code = 0;
    for(const auto& I : j_material) code += I->commit_status();
    return code;
}

int Joint::clear_status() {
    auto code = 0;
    for(const auto& I : j_material) code += I->clear_status();
    return code;
}

int Joint::reset_status() {
    auto code = 0;
    for(const auto& I : j_material) code += I->reset_status();
    return code;
}

vector<vec> Joint::record(const OutputType P) {
    vector<vec> data;

    for(const auto& I : j_material) for(auto& J : I->record(P)) data.emplace_back(J);

    return data;
}

void Joint::print() { suanpan_info("A joint element that uses displacement as basic quantity. The material model used shall be based on displacement--force relationship.\n"); }
