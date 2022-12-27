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

#include "Rectangle1D.h"
#include <Material/Material.h>

Rectangle1D::Rectangle1D(const unsigned T, const double B, const double H, const unsigned M)
    : Section1D(T, M, B * H)
    , width(B)
    , height(H) {}

unique_ptr<Section> Rectangle1D::get_copy() { return make_unique<Rectangle1D>(*this); }

void Rectangle1D::print() {
    suanpan_info("A 1D rectangle section with width %.3E and height %.3E.\n", width, height);
    suanpan_info("Material: ");
    s_material->print();
}
