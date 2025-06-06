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

#include "SupportMotion.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Load/Amplitude/Amplitude.h>
#include <Solver/Integrator/Integrator.h>
#include <Step/Step.h>

SupportMotion::SupportMotion(const unsigned T, const double L, uvec&& N, uvec&& D, const unsigned AT)
    : Load(T, AT, std::move(N), std::move(D), L) { enable_displacement_control(); }

int SupportMotion::initialize(const shared_ptr<DomainBase>& D) {
    set_end_step(start_step + 1);

    D->get_factory()->update_reference_dof(encoding = get_nodal_active_dof(D));

    return Load::initialize(D);
}

int SupportDisplacement::process(const shared_ptr<DomainBase>& D) {
    const auto& W = D->get_factory();

    trial_settlement.zeros(W->get_size());

    trial_settlement(encoding).fill(pattern * amplitude->get_amplitude(W->get_trial_time()));

    return SUANPAN_SUCCESS;
}

int SupportVelocity::process(const shared_ptr<DomainBase>& D) {
    const auto& W = D->get_factory();
    const auto& G = D->get_current_step()->get_integrator();

    trial_settlement.zeros(W->get_size());

    trial_settlement(encoding) = G->from_total_velocity(pattern * amplitude->get_amplitude(W->get_trial_time()), encoding);

    return SUANPAN_SUCCESS;
}

int SupportAcceleration::process(const shared_ptr<DomainBase>& D) {
    const auto& W = D->get_factory();
    const auto& G = D->get_current_step()->get_integrator();

    trial_settlement.zeros(W->get_size());

    trial_settlement(encoding) = G->from_total_acceleration(pattern * amplitude->get_amplitude(W->get_trial_time()), encoding);

    return SUANPAN_SUCCESS;
}
