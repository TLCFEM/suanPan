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

#include "NodalDisplacement.h"

#include <Domain/Factory.hpp>

NodalDisplacement::NodalDisplacement(const unsigned T, const double L, uvec&& N, std::vector<Node::DOF>&& D, const unsigned AT)
    : Load(T, AT, {}, std::move(D), L) { target_node = std::move(N); }

int NodalDisplacement::initialize(const shared_ptr<DomainBase>& D) {
    if(SUANPAN_SUCCESS != Load::initialize(D)) return SUANPAN_FAIL;

    set_end_step(start_step + 1);

    D->get_factory()->update_reference_dof(target_node_dof);

    return SUANPAN_SUCCESS;
}

int NodalDisplacement::process(const shared_ptr<DomainBase>& D) {
    if(target_node_dof.empty()) trial_settlement.reset();
    else trial_settlement.zeros(D->get_factory()->get_size())(target_node_dof).fill(magnitude * get_amplitude(D));

    return SUANPAN_SUCCESS;
}

GroupNodalDisplacement::GroupNodalDisplacement(const unsigned T, const double L, uvec&& N, std::vector<Node::DOF>&& D, const unsigned AT)
    : GroupModifier(std::move(N))
    , NodalDisplacement(T, L, {}, std::move(D), AT) {}

int GroupNodalDisplacement::initialize(const shared_ptr<DomainBase>& D) {
    target_node = update_object_tag(D);

    return NodalDisplacement::initialize(D);
}
