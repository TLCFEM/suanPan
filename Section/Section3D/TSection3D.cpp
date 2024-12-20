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

#include "TSection3D.h"
#include <Material/Material1D/Material1D.h>

TSection3D::TSection3D(const unsigned T, const double TFW, const double TFT, const double WH, const double WT, const unsigned MT, const unsigned IP, vec&& EC)
    : ISection3D(T, TFW, TFT, 0., 0., WH, WT, MT, IP, std::move(EC)) {}

TSection3D::TSection3D(const unsigned T, vec&& D, const unsigned MT, const unsigned IP, vec&& EC)
    : ISection3D(T, D(0), D(1), 0., 0., D(2), D(3), MT, IP, std::move(EC)) {}

unique_ptr<Section> TSection3D::get_copy() { return make_unique<TSection3D>(*this); }
