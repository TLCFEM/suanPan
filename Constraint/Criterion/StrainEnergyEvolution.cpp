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

#include "StrainEnergyEvolution.h"
#include <Element/Element.h>

StrainEnergyEvolution::StrainEnergyEvolution(const unsigned T, const unsigned ST, const unsigned IL, const unsigned FL, const double WT, const unsigned IT, const unsigned RR, const double PW, const double TL)
    : EnergyEvolution(T, ST, IL, FL, WT, IT, RR, PW, TL) { get_energy = &Element::get_strain_energy; }

unique_ptr<Criterion> StrainEnergyEvolution::get_copy() { return make_unique<StrainEnergyEvolution>(*this); }
