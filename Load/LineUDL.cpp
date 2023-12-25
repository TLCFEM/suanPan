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

#include "LineUDL.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>
#include <Load/Amplitude/Amplitude.h>

LineUDL::LineUDL(const unsigned T, const unsigned S, const double L, uvec&& N, const unsigned DT, const unsigned AT, const uword D)
    : Load(T, S, AT, std::move(N), uvec{DT}, L)
    , dimension(D) {}

int LineUDL::initialize(const shared_ptr<DomainBase>& D) {
    if(node_encoding.n_elem % 2 != 0) return SUANPAN_FAIL;

    for(const auto I : node_encoding) if(const auto& t_node = D->get<Node>(I); t_node->get_reordered_dof().size() < dimension || t_node->get_coordinate().size() < dimension) return SUANPAN_FAIL;

    return Load::initialize(D);
}

LineUDL2D::LineUDL2D(const unsigned T, const unsigned S, const double L, uvec&& N, const unsigned DT, const unsigned AT)
    : LineUDL(T, S, L, std::move(N), DT, AT, 2llu) {}

int LineUDL2D::process(const shared_ptr<DomainBase>& D) {
    const auto& W = D->get_factory();

    trial_load.zeros(W->get_size());

    const auto ref_load = pattern * magnitude->get_amplitude(W->get_trial_time());

    for(auto I = 0llu, J = 1llu; J < node_encoding.n_elem; ++I, ++J) {
        const auto& node_i = D->get<Node>(node_encoding(I));
        const auto& node_j = D->get<Node>(node_encoding(J));
        const auto& dof_i = node_i->get_reordered_dof();
        const auto& dof_j = node_j->get_reordered_dof();

        const vec diff_coor = node_j->get_coordinate().head(dimension) - node_i->get_coordinate().head(dimension);

        if(0llu == dof_reference(0)) {
            trial_load(dof_i(0)) = trial_load(dof_j(0)) = -.5 * diff_coor(1) * ref_load;
            D->insert_loaded_dof(dof_i(0));
            D->insert_loaded_dof(dof_j(0));
        }
        else if(1llu == dof_reference(0)) {
            trial_load(dof_i(1)) = trial_load(dof_j(1)) = -.5 * diff_coor(0) * ref_load;
            D->insert_loaded_dof(dof_i(1));
            D->insert_loaded_dof(dof_j(1));
        }
    }

    return SUANPAN_SUCCESS;
}

LineUDL3D::LineUDL3D(const unsigned T, const unsigned S, const double L, uvec&& N, const unsigned DT, const unsigned AT)
    : LineUDL(T, S, L, std::move(N), DT, AT, 3llu) {}

int LineUDL3D::process(const shared_ptr<DomainBase>& D) {
    const auto& W = D->get_factory();

    trial_load.zeros(W->get_size());

    const auto ref_load = pattern * magnitude->get_amplitude(W->get_trial_time());

    for(auto I = 0llu, J = 1llu; J < node_encoding.n_elem; ++I, ++J) {
        const auto& node_i = D->get<Node>(node_encoding(I));
        const auto& node_j = D->get<Node>(node_encoding(J));
        const auto& dof_i = node_i->get_reordered_dof();
        const auto& dof_j = node_j->get_reordered_dof();

        const vec diff_coor = node_j->get_coordinate().head(dimension) - node_i->get_coordinate().head(dimension);

        if(0llu == dof_reference(0)) {
            trial_load(dof_i(0)) = trial_load(dof_j(0)) = -.5 * norm(diff_coor(uvec{1, 2})) * ref_load;
            D->insert_loaded_dof(dof_i(0));
            D->insert_loaded_dof(dof_j(0));
        }
        else if(1llu == dof_reference(0)) {
            trial_load(dof_i(1)) = trial_load(dof_j(1)) = -.5 * norm(diff_coor(uvec{0, 2})) * ref_load;
            D->insert_loaded_dof(dof_i(1));
            D->insert_loaded_dof(dof_j(1));
        }
        else if(2llu == dof_reference(0)) {
            trial_load(dof_i(2)) = trial_load(dof_j(2)) = -.5 * norm(diff_coor(uvec{0, 1})) * ref_load;
            D->insert_loaded_dof(dof_i(2));
            D->insert_loaded_dof(dof_j(2));
        }
    }

    return SUANPAN_SUCCESS;
}
