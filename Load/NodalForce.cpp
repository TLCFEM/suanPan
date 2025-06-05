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

#include "NodalForce.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Load/Amplitude/Amplitude.h>

NodalForce::NodalForce(const unsigned T, const unsigned S, const double L, uvec&& N, const unsigned D, const unsigned AT)
    : Load(T, S, AT, std::move(N), uvec{D}, L) {}

NodalForce::NodalForce(const unsigned T, const unsigned S, const double L, uvec&& N, uvec&& D, const unsigned AT)
    : Load(T, S, AT, std::move(N), std::move(D), L) {}

int NodalForce::process(const shared_ptr<DomainBase>& D) {
    const auto& W = D->get_factory();

    const auto active_dof = get_nodal_active_dof(D);

    D->insert_loaded_dof(active_dof);

    trial_load.zeros(W->get_size());

    trial_load(active_dof).fill(pattern * amplitude->get_amplitude(W->get_trial_time()));

    return SUANPAN_SUCCESS;
}
