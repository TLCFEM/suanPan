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

#include "LineUDL.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

LineUDL::LineUDL(const unsigned T, const double L, uvec&& N, std::vector<Node::DOF>&& DT, const unsigned AT, const unsigned D)
    : Load(T, AT, std::move(N), 2u == D ? std::vector{Node::DOF::U1, Node::DOF::U2} : std::vector{Node::DOF::U1, Node::DOF::U2, Node::DOF::U3}, std::move(DT), L)
    , dimension(D) {}

LineUDL2D::LineUDL2D(const unsigned T, const double L, uvec&& N, std::vector<Node::DOF>&& DT, const unsigned AT)
    : LineUDL(T, L, std::move(N), std::move(DT), AT, 2u) {}

int LineUDL2D::process(const shared_ptr<DomainBase>& D) {
    const auto& W = D->get_factory();

    trial_load.zeros(W->get_size());

    const auto ref_load = magnitude * get_amplitude(D);

    for(auto I = 0llu, J = 1llu; J < target_node.n_elem; ++I, ++J) {
        const auto& node_i = D->get<Node>(target_node(I));
        const auto& node_j = D->get<Node>(target_node(J));
        const auto& dof_i = node_i->get_reordered_dof();
        const auto& dof_j = node_j->get_reordered_dof();

        const vec diff_coor = node_j->initial_position(dimension) - node_i->initial_position(dimension);

        if(1u == static_cast<unsigned>(get_dof_component()[0])) {
            trial_load(dof_i(0)) = trial_load(dof_j(0)) = -.5 * diff_coor(1) * ref_load;
            D->insert_loaded_dof(dof_i(0));
            D->insert_loaded_dof(dof_j(0));
        }
        else if(2u == static_cast<unsigned>(get_dof_component()[0])) {
            trial_load(dof_i(1)) = trial_load(dof_j(1)) = -.5 * diff_coor(0) * ref_load;
            D->insert_loaded_dof(dof_i(1));
            D->insert_loaded_dof(dof_j(1));
        }
    }

    return SUANPAN_SUCCESS;
}

LineUDL3D::LineUDL3D(const unsigned T, const double L, uvec&& N, std::vector<Node::DOF>&& DT, const unsigned AT)
    : LineUDL(T, L, std::move(N), std::move(DT), AT, 3u) {}

int LineUDL3D::process(const shared_ptr<DomainBase>& D) {
    const auto& W = D->get_factory();

    trial_load.zeros(W->get_size());

    const auto ref_load = magnitude * get_amplitude(D);

    for(auto I = 0llu, J = 1llu; J < target_node.n_elem; ++I, ++J) {
        const auto& node_i = D->get<Node>(target_node(I));
        const auto& node_j = D->get<Node>(target_node(J));
        const auto& dof_i = node_i->get_reordered_dof();
        const auto& dof_j = node_j->get_reordered_dof();

        const vec diff_coor = node_j->initial_position(dimension) - node_i->initial_position(dimension);

        if(1u == static_cast<unsigned>(get_dof_component()[0])) {
            trial_load(dof_i(0)) = trial_load(dof_j(0)) = -.5 * norm(diff_coor(uvec{1, 2})) * ref_load;
            D->insert_loaded_dof(dof_i(0));
            D->insert_loaded_dof(dof_j(0));
        }
        else if(2u == static_cast<unsigned>(get_dof_component()[0])) {
            trial_load(dof_i(1)) = trial_load(dof_j(1)) = -.5 * norm(diff_coor(uvec{0, 2})) * ref_load;
            D->insert_loaded_dof(dof_i(1));
            D->insert_loaded_dof(dof_j(1));
        }
        else if(3u == static_cast<unsigned>(get_dof_component()[0])) {
            trial_load(dof_i(2)) = trial_load(dof_j(2)) = -.5 * norm(diff_coor(uvec{0, 1})) * ref_load;
            D->insert_loaded_dof(dof_i(2));
            D->insert_loaded_dof(dof_j(2));
        }
    }

    return SUANPAN_SUCCESS;
}
