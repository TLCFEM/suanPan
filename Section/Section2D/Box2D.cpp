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

#include "Box2D.h"
#include <Material/Material1D/Material1D.h>

Box2D::Box2D(const unsigned T, const double B, const double H, const double TH, const unsigned M, const unsigned S, const double EC)
    : ISection2D(T, B, TH, B, TH, H - 2. * TH, 2. * TH, M, S, EC) {}

Box2D::Box2D(const unsigned T, vec&& D, const unsigned M, const unsigned S, const double EC)
    : ISection2D(T, D(0), D(2), D(0), D(2), D(1) - 2. * D(2), 2. * D(2), M, S, EC) {}

unique_ptr<Section> Box2D::get_copy() { return make_unique<Box2D>(*this); }

void Box2D::print() {
    suanpan_info("A 2D box section.\n");
}
