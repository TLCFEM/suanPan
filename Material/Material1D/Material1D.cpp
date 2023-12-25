/*******************************************************************************
 * Copyright (C) 2017-2024 Theodore Chang
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

#include "Material1D.h"
#include <Recorder/OutputType.h>

Material1D::Material1D(const unsigned T, const double D)
    : Material(T, MaterialType::D1, D) { set_symmetric(true); }

vector<vec> Material1D::record(const OutputType P) {
    if(P == OutputType::S11) return {current_stress};
    if(P == OutputType::E11) return {current_strain};

    return Material::record(P);
}

void Material1D::print() {
    suanpan_info("Current Strain: {:.3E}\tCurrent Stress: {:.3E}\n", current_strain(0), current_stress(0));
}
