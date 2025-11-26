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

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>

MPC::MPC(const unsigned T, const unsigned A, const double L, std::vector<std::tuple<uword, uword, double>>&& P)
    : Constraint(T, A, {}, {}, {}, 1)
    , magnitude(L)
    , pool(std::move(P)) {}

int MPC::initialize(const shared_ptr<DomainBase>& D) {
    auto& W = D->get_factory();

    auxiliary_stiffness.zeros(W->get_size(), lagrangian_size);

    for(auto [tag, dof, weight] : pool) {
        auto& node = D->get<Node>(tag);
        if(!node || !node->is_active() || node->get_reordered_dof().size() <= dof) {
            auxiliary_stiffness.reset();
            return SUANPAN_FAIL;
        }
        auxiliary_stiffness(node->get_reordered_dof()(dof - 1)) = weight;
    }

    return Constraint::initialize(D);
}

int MPC::process(const shared_ptr<DomainBase>& D) {
    auxiliary_load = magnitude * get_amplitude(D);

    auxiliary_resistance = auxiliary_stiffness.t() * D->get_factory()->get_trial_displacement();

    return SUANPAN_SUCCESS;
}
