/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
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
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Load/Amplitude/Amplitude.h>

NodalAcceleration::NodalAcceleration(const unsigned T, const unsigned ST, const double L, const unsigned DT, const unsigned AT)
    : Load(T, ST, AT, {}, uvec{DT}, L) {}

NodalAcceleration::NodalAcceleration(const unsigned T, const unsigned ST, const double L, uvec&& DT, const unsigned AT)
    : Load(T, ST, AT, {}, std::forward<uvec>(DT), L) {}

NodalAcceleration::NodalAcceleration(const unsigned T, const unsigned ST, const double L, uvec&& NT, const unsigned DT, const unsigned AT)
    : Load(T, ST, AT, std::forward<uvec>(NT), uvec{DT}, L) {}

int NodalAcceleration::process(const shared_ptr<DomainBase>& D) {
    auto& W = D->get_factory();

    trial_load.reset();

    if(nullptr == W->get_mass()) return SUANPAN_SUCCESS;

    trial_load.zeros(W->get_size());

    trial_load(node_encoding.is_empty() ? get_all_nodal_active_dof(D) : get_nodal_active_dof(D)).fill(1.);

    trial_load = W->get_mass() * trial_load * pattern * magnitude->get_amplitude(W->get_trial_time());

    return SUANPAN_SUCCESS;
}
