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

#include "Parallel.h"
#include <Domain/DomainBase.h>

Parallel::Parallel(const unsigned T, uvec&& MT)
    : Material1D(T, 0.)
    , mat_tag(std::forward<uvec>(MT)) {}

Parallel::Parallel(const Parallel& old_obj)
    : Material1D(old_obj)
    , mat_tag(old_obj.mat_tag) { for(const auto& I : old_obj.mat_pool) mat_pool.emplace_back(I->get_copy()); }

int Parallel::initialize(const shared_ptr<DomainBase>& D) {
    mat_pool.clear();
    mat_pool.reserve(mat_tag.n_elem);
    for(const auto I : mat_tag) {
        mat_pool.emplace_back(suanpan::initialized_material_copy(D, I));
        if(nullptr == mat_pool.back() || mat_pool.back()->get_material_type() != MaterialType::D1) {
            suanpan_error("Parallel %u requires 1D host material models.\n", get_tag());
            return SUANPAN_FAIL;
        }
        access::rw(density) += mat_pool.back()->get_parameter(ParameterType::DENSITY);
    }

    initial_stiffness.zeros(1);
    initial_damping.zeros(1);
    for(const auto& I : mat_pool) {
        if(!I->get_initial_stiffness().empty()) initial_stiffness += I->get_initial_stiffness();
        if(!I->get_initial_damping().empty()) initial_damping += I->get_initial_damping();
    }
    trial_stiffness = current_stiffness = initial_stiffness;
    trial_damping = current_damping = initial_damping;

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Parallel::get_copy() { return make_unique<Parallel>(*this); }

int Parallel::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress.zeros();
    trial_stiffness.zeros();
    for(const auto& I : mat_pool) {
        if(I->update_trial_status(trial_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
        trial_stress += I->get_trial_stress();
        if(!I->get_trial_stiffness().empty()) trial_stiffness += I->get_trial_stiffness();
    }

    return SUANPAN_SUCCESS;
}

int Parallel::update_trial_status(const vec& t_strain, const vec& t_strain_rate) {
    incre_strain = (trial_strain = t_strain) - current_strain;
    incre_strain_rate = (trial_strain_rate = t_strain_rate) - current_strain_rate;

    if(fabs(incre_strain(0) + incre_strain_rate(0)) <= datum::eps) return SUANPAN_SUCCESS;

    trial_stress.zeros();
    trial_stiffness.zeros();
    trial_damping.zeros();
    for(const auto& I : mat_pool) {
        if(I->update_trial_status(trial_strain, trial_strain_rate) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
        trial_stress += I->get_trial_stress();
        if(!I->get_trial_stiffness().empty()) trial_stiffness += I->get_trial_stiffness();
        if(!I->get_trial_damping().empty()) trial_damping += I->get_trial_damping();
    }

    return SUANPAN_SUCCESS;
}

int Parallel::clear_status() {
    current_strain.zeros();
    trial_strain.zeros();
    current_stress.zeros();
    trial_stress.zeros();
    trial_stiffness = current_stiffness = initial_stiffness;
    auto code = 0;
    for(const auto& I : mat_pool) code += I->clear_status();
    return code;
}

int Parallel::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    auto code = 0;
    for(const auto& I : mat_pool) code += I->commit_status();
    return code;
}

int Parallel::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    auto code = 0;
    for(const auto& I : mat_pool) code += I->reset_status();
    return code;
}

vector<vec> Parallel::record(const OutputType P) {
    vector<vec> data;

    auto max_size = 0llu;
    for(const auto& I : mat_pool)
        for(const auto& J : I->record(P)) {
            if(J.n_elem > max_size) max_size = J.n_elem;
            data.emplace_back(J);
        }

    for(auto&& I : data) I.resize(max_size);

    return data;
}

void Parallel::print() {
    mat_tag.t().print("A Parallel container that holds following models:");
    for(const auto& I : mat_pool) I->print();
}
