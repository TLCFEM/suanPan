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

#include "Material2D.h"
#include <Recorder/OutputType.h>
#include <Toolbox/tensor.h>

Material2D::Material2D(const unsigned T, const PlaneType PT, const double R)
    : Material(T, MaterialType::D2, R)
    , plane_type(PT) {}

vector<vec> Material2D::record(const OutputType P) {
    vector<vec> data;

    if(P == OutputType::S) data.emplace_back(current_stress);
    else if(P == OutputType::SP) data.emplace_back(transform::stress::principal(current_stress));
    else if(P == OutputType::S11) data.emplace_back(vec{current_stress(0)});
    else if(P == OutputType::S22) data.emplace_back(vec{current_stress(1)});
    else if(P == OutputType::S12) data.emplace_back(vec{current_stress(2)});
    else if(P == OutputType::E) data.emplace_back(current_strain);
    else if(P == OutputType::EP) data.emplace_back(transform::strain::principal(current_strain));
    else if(P == OutputType::E11) data.emplace_back(vec{current_strain(0)});
    else if(P == OutputType::E22) data.emplace_back(vec{current_strain(1)});
    else if(P == OutputType::E12) data.emplace_back(vec{current_strain(2)});

    return data;
}
