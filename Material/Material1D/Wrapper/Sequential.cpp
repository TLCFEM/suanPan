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

#include "Sequential.h"
#include <Domain/DomainBase.h>

constexpr unsigned Sequential::max_iteration = 20;

Sequential::Sequential(const unsigned T, uvec&& MT)
    : Material1D(T, 0.)
    , mat_size(MT.n_elem - 1)
    , mat_tag(std::forward<uvec>(MT)) {}

Sequential::Sequential(const Sequential& old_obj)
    : Material1D(old_obj)
    , mat_size(old_obj.mat_size)
    , mat_tag(old_obj.mat_tag)
    , jacobian(old_obj.jacobian) { for(const auto& I : old_obj.mat_pool) mat_pool.emplace_back(I->get_copy()); }

int Sequential::initialize(const shared_ptr<DomainBase>& D) {
    mat_pool.clear();
    mat_pool.reserve(mat_tag.n_elem);
    for(const auto I : mat_tag) {
        mat_pool.emplace_back(suanpan::initialized_material_copy(D, I));
        if(nullptr == mat_pool.back() || mat_pool.back()->get_material_type() != MaterialType::D1) {
            suanpan_error("A valid 1D host material is required.\n");
            return SUANPAN_FAIL;
        }
        access::rw(density) += mat_pool.back()->get_parameter(ParameterType::DENSITY);
    }

    jacobian.zeros(mat_tag.n_elem, mat_tag.n_elem);
    jacobian.row(0).fill(1.);

    initial_stiffness.zeros(1);
    for(const auto& I : mat_pool) initial_stiffness += 1. / I->get_initial_stiffness();
    trial_stiffness = current_stiffness = initial_stiffness = 1. / initial_stiffness;

    return SUANPAN_SUCCESS;
}

unique_ptr<Material> Sequential::get_copy() { return make_unique<Sequential>(*this); }

int Sequential::update_trial_status(const vec& t_strain) {
    incre_strain = (trial_strain = t_strain) - current_strain;

    if(fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

    unsigned counter = 0;
    while(++counter < max_iteration) {
        vec residual(mat_tag.n_elem, fill::zeros);
        residual(0) = trial_strain(0) - mat_pool.front()->get_trial_strain().at(0) - mat_pool.back()->get_trial_strain().at(0);
        residual(1) = mat_pool.front()->get_trial_stress().at(0);
        residual(mat_size) -= mat_pool.back()->get_trial_stress().at(0);

        jacobian(1, 0) = -mat_pool.front()->get_trial_stiffness().at(0);
        jacobian(mat_size, mat_size) = mat_pool.back()->get_trial_stiffness().at(0);

        for(uword I = 1; I < mat_size; ++I) {
            residual(0) -= mat_pool[I]->get_trial_strain().at(0);
            residual(I) -= mat_pool[I]->get_trial_stress().at(0);
            residual(I + 1) += mat_pool[I]->get_trial_stress().at(0);
            jacobian(I + 1, I) = -(jacobian(I, I) = mat_pool[I]->get_trial_stiffness().at(0));
        }

        const vec i_strain = solve(jacobian, residual);

        for(size_t I = 0; I < mat_pool.size(); ++I) mat_pool[I]->update_trial_status(mat_pool[I]->get_trial_strain() + i_strain[I]);

        const auto error = norm(i_strain);
        suanpan_debug("Local iteration error: {:.5E}.\n", error);
        if(error <= tolerance) break;
    }

    if(max_iteration == counter) {
        suanpan_error("Cannot converge within {} iterations.\n", max_iteration);
        return SUANPAN_FAIL;
    }

    trial_stress = mat_pool.front()->get_trial_stress();

    auto t_stiff = 0.;
    for(const auto& I : mat_pool) t_stiff += 1. / I->get_trial_stiffness().at(0);

    trial_stiffness = 1. / t_stiff;

    return SUANPAN_SUCCESS;
}

int Sequential::clear_status() {
    current_strain.zeros();
    trial_strain.zeros();
    current_stress.zeros();
    trial_stress.zeros();
    trial_stiffness = current_stiffness = initial_stiffness;
    auto code = 0;
    for(const auto& I : mat_pool) code += I->clear_status();
    return code;
}

int Sequential::commit_status() {
    current_strain = trial_strain;
    current_stress = trial_stress;
    current_stiffness = trial_stiffness;
    auto code = 0;
    for(const auto& I : mat_pool) code += I->commit_status();
    return code;
}

int Sequential::reset_status() {
    trial_strain = current_strain;
    trial_stress = current_stress;
    trial_stiffness = current_stiffness;
    auto code = 0;
    for(const auto& I : mat_pool) code += I->reset_status();
    return code;
}

vector<vec> Sequential::record(const OutputType P) {
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

void Sequential::print() {
    suanpan_info("A wrapper of several material models.\n");
    for(const auto& I : mat_pool) I->print();
}
