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

#include "NodalAcceleration.h"

#include <Domain/Factory.hpp>

NodalAcceleration::NodalAcceleration(const unsigned T, const double L, uvec&& NT, std::vector<Node::DOF>&& DT, const unsigned AT)
    : Load(T, AT, {}, std::move(DT), L) { target_node = std::move(NT); }

int NodalAcceleration::initialize(const shared_ptr<DomainBase>& D) {
    if(SUANPAN_SUCCESS != Load::initialize(D)) return SUANPAN_FAIL;

    target_node_dof = collect_node_dof(D);

    return SUANPAN_SUCCESS;
}

int NodalAcceleration::process(const shared_ptr<DomainBase>& D) {
    auto& W = D->get_factory();

    if(auto& t_mass = W->get_mass(); t_mass && !target_node_dof.is_empty()) {
        trial_load.zeros(W->get_size())(target_node_dof).fill(magnitude * get_amplitude(D));
        trial_load = t_mass * trial_load;
    }
    else trial_load.reset();

    return SUANPAN_SUCCESS;
}
