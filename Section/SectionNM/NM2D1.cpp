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

#include "NM2D1.h"

#include <Domain/DomainBase.h>

NM2D1::NM2D1(const unsigned T, const double EEA, const double EEIS, const double LD)
    : SectionNM2D(T, EEA, EEIS, LD) {}

unique_ptr<Section> NM2D1::get_copy() { return make_unique<NM2D1>(*this); }

int NM2D1::update_trial_status(const vec& t_deformation) {
    trial_resistance = trial_stiffness * (trial_deformation = t_deformation);

    return SUANPAN_SUCCESS;
}

void NM2D1::print() {
    suanpan_info("A NM2D1 elastic section.\n");
}
