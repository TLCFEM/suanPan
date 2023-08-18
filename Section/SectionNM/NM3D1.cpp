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

#include "NM3D1.h"
#include <Domain/DomainBase.h>

NM3D1::NM3D1(const unsigned T, const double EEA, const double EEIS, const double EEIW, const double LD)
    : SectionNM3D(T, EEA, EEIS, EEIW, LD) {}

unique_ptr<Section> NM3D1::get_copy() { return make_unique<NM3D1>(*this); }

int NM3D1::update_trial_status(const vec& t_deformation) {
    trial_resistance = trial_stiffness * (trial_deformation = t_deformation);

    return SUANPAN_SUCCESS;
}

void NM3D1::print() {
    suanpan_info("A NM3D1 elastic section.\n");
}
