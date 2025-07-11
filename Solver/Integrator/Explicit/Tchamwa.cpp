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

#include "Tchamwa.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>

Tchamwa::Tchamwa(const unsigned T, const double R)
    : ExplicitIntegrator(T)
    , PHI(2. / (1. + std::max(0., std::min(1., R)))) {}

int Tchamwa::update_trial_status(bool) {
    const auto D = get_domain();
    auto& W = D->get_factory();

    W->update_incre_displacement(DT * W->get_current_velocity() + PHI * DT * DT * W->get_current_acceleration());
    W->update_incre_velocity(DT * W->get_current_acceleration());

    return D->update_trial_status();
}

void Tchamwa::update_parameter(const double NT) { DT = NT; }

void Tchamwa::print() {
    suanpan_info("A Tchamwa solver.\n");
}
