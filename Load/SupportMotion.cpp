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

#include "SupportMotion.h"

#include <Domain/Factory.hpp>
#include <Solver/Integrator/Integrator.h>
#include <Step/Step.h>

SupportMotion::SupportMotion(const unsigned T, const double L, uvec&& N, std::vector<Node::DOF>&& D, const unsigned AT)
    : Load(T, AT, {}, std::move(D), L) { target_node = std::move(N); }

int SupportMotion::initialize(const shared_ptr<DomainBase>& D) {
    if(SUANPAN_SUCCESS != Load::initialize(D)) return SUANPAN_FAIL;

    set_end_step(start_step + 1);

    D->get_factory()->update_reference_dof(target_node_dof);

    return SUANPAN_SUCCESS;
}

int SupportDisplacement::process(const shared_ptr<DomainBase>& D) {
    if(target_node_dof.is_empty()) trial_settlement.reset();
    else trial_settlement.zeros(D->get_factory()->get_size())(target_node_dof).fill(magnitude * get_amplitude(D));

    return SUANPAN_SUCCESS;
}

int SupportVelocity::process(const shared_ptr<DomainBase>& D) {
    const auto& W = D->get_factory();
    const auto& G = D->get_current_step()->get_integrator();

    if(target_node_dof.is_empty()) trial_settlement.reset();
    else trial_settlement.zeros(W->get_size())(target_node_dof) = G->from_total_velocity(magnitude * get_amplitude(D), target_node_dof);

    return SUANPAN_SUCCESS;
}

int SupportAcceleration::process(const shared_ptr<DomainBase>& D) {
    const auto& W = D->get_factory();
    const auto& G = D->get_current_step()->get_integrator();

    if(target_node_dof.is_empty()) trial_settlement.reset();
    else trial_settlement.zeros(W->get_size())(target_node_dof) = G->from_total_acceleration(magnitude * get_amplitude(D), target_node_dof);

    return SUANPAN_SUCCESS;
}
