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

#include "TrussSection.h"
#include <Material/Material1D/Material1D.h>

TrussSection::TrussSection(const unsigned T, const double A, const unsigned M)
    : Section1D(T, M, A) {}

unique_ptr<Section> TrussSection::get_copy() { return make_unique<TrussSection>(*this); }

void TrussSection::print() { suanpan_info("A 1D generalized section with an area of %.3E.\n", area); }
