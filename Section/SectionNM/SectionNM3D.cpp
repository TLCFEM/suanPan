/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
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

#include "SectionNM3D.h"
#include <Domain/DomainBase.h>

SectionNM3D::SectionNM3D(const unsigned T, const double EEA, const double EEIS, const double EEIW, const double LD)
    : DataSectionNM3D{EEA, EEIS, EEIW}
    , SectionNM(T, SectionType::NM3D) { access::rw(linear_density) = LD; }

int SectionNM3D::initialize(const shared_ptr<DomainBase>&) {
    initial_stiffness.zeros(6, 6);

    initial_stiffness(0, 0) = EA;
    initial_stiffness(5, 5) = ET;

    initial_stiffness(1, 1) = initial_stiffness(2, 2) = 2. * (initial_stiffness(1, 2) = initial_stiffness(2, 1) = 2. * EIS);
    initial_stiffness(3, 3) = initial_stiffness(4, 4) = 2. * (initial_stiffness(3, 4) = initial_stiffness(4, 3) = 2. * EIW);

    trial_stiffness = current_stiffness = initial_stiffness;

    return SUANPAN_SUCCESS;
}
