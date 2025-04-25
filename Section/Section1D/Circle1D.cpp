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

#include "Circle1D.h"

#include <Material/Material.h>

Circle1D::Circle1D(const unsigned T, const double B, const unsigned M)
    : Section1D(T, M, B * B * datum::pi)
    , radius(B) {}

unique_ptr<Section> Circle1D::get_copy() { return make_unique<Circle1D>(*this); }

void Circle1D::print() {
    suanpan_info("A uniaxial circular section with radius of {:.3E}.\nMaterial: ", radius);
    s_material->print();
}
