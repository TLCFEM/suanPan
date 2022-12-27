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

#include "TSection2D.h"
#include <Material/Material1D/Material1D.h>

TSection2D::TSection2D(const unsigned T, const double TFW, const double TFT, const double WH, const double WT, const unsigned MT, const unsigned IP, const double EC)
    : ISection2D(T, TFW, TFT, 0., 0., WH, WT, MT, IP, EC) {}

TSection2D::TSection2D(const unsigned T, vec&& D, const unsigned MT, const unsigned IP, const double EC)
    : ISection2D(T, D(0), D(1), 0., 0., D(2), D(3), MT, IP, EC) {}

unique_ptr<Section> TSection2D::get_copy() { return make_unique<TSection2D>(*this); }
