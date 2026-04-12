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

#include "NodalForce.h"

#include <Domain/Factory.hpp>

NodalForce::NodalForce(const unsigned T, const double L, uvec&& N, std::vector<Node::DOF>&& D, const unsigned AT)
    : Load(T, AT, {}, std::move(D), L) { target_node = std::move(N); }

int NodalForce::initialize(const shared_ptr<DomainBase>& D) {
    if(SUANPAN_SUCCESS != Load::initialize(D)) return SUANPAN_FAIL;

    target_node_dof = collect_node_dof(D);

    return SUANPAN_SUCCESS;
}

int NodalForce::process(const shared_ptr<DomainBase>& D) {
    D->insert_loaded_dof(target_node_dof);

    if(target_node_dof.is_empty()) trial_load.reset();
    else trial_load.zeros(D->get_factory()->get_size())(target_node_dof).fill(magnitude * get_amplitude(D));

    return SUANPAN_SUCCESS;
}

GroupNodalForce::GroupNodalForce(const unsigned T, const double L, uvec&& N, std::vector<Node::DOF>&& D, const unsigned AT)
    : GroupModifier(std::move(N))
    , NodalForce(T, L, {}, std::move(D), AT) {}

int GroupNodalForce::initialize(const shared_ptr<DomainBase>& D) {
    target_node = update_object_tag(D);

    return NodalForce::initialize(D);
}
