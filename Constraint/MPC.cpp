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

#include "MPC.h"

#include <Domain/Factory.hpp>

MPC::MPC(const unsigned T, const unsigned A, const double L, Pack&& P)
    : Constraint(T, A, {}, {}, 1)
    , magnitude(L)
    , pool(std::move(P)) {}

int MPC::initialize(const shared_ptr<DomainBase>& D) {
    auxiliary_stiffness.zeros(D->get_factory()->get_size(), lagrangian_size);

    for(auto [tag, target, weight] : pool) {
        if(auto& node = D->get<Node>(tag); node && node->is_active())
            if(const auto global = node->get_dof({target}); !global.empty()) {
                auxiliary_stiffness(global.front()) = weight;
                continue;
            }
        auxiliary_stiffness.reset();
        return SUANPAN_FAIL;
    }

    return Constraint::initialize(D);
}

int MPC::process(const shared_ptr<DomainBase>& D) {
    auxiliary_load = magnitude * get_amplitude(D);

    auxiliary_resistance = auxiliary_stiffness.t() * D->get_factory()->get_trial_displacement();

    return SUANPAN_SUCCESS;
}
