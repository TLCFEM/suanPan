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

#include "NodalAcceleration.h"

#include <Domain/Factory.hpp>

NodalAcceleration::NodalAcceleration(const unsigned T, const double L, uvec&& NT, std::vector<Node::DOF>&& DT, const unsigned AT)
    : Load(T, AT, std::move(NT), {}, std::move(DT), L) {}

int NodalAcceleration::process(const shared_ptr<DomainBase>& D) {
    auto& W = D->get_factory();

    trial_load.reset();

    if(auto& t_mass = W->get_mass()) {
        trial_load.zeros(W->get_size())(target_dof).fill(magnitude * get_amplitude(D));

        trial_load = t_mass * trial_load;
    }

    return SUANPAN_SUCCESS;
}
