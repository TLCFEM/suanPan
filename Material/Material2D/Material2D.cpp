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

#include "Material2D.h"

#include <Recorder/OutputType.h>
#include <Toolbox/tensor.h>

Material2D::Material2D(const unsigned T, const PlaneType PT, const double R)
    : Material(T, MaterialType::D2, R) { access::rw(plane_type) = PT; }

std::vector<vec> Material2D::record(const OutputType P) const {
    if(P == OutputType::HIST) return {current_history};
    if(P == OutputType::YF) return {vec{any(current_history != 0.) ? 1. : 0.}};

    const auto remap = [](const vec& in) {
        vec out(6, fill::zeros);
        out(uvec{0, 1, 3}) = in;
        return out;
    };

    if(plane_type == PlaneType::S) {
        if(P == OutputType::S) return {remap(current_stress)};
        if(P == OutputType::SP) return {transform::stress::principal(current_stress)};
    }
    else if(plane_type == PlaneType::E) {
        if(P == OutputType::E) return {remap(current_strain)};
        if(P == OutputType::EP) return {transform::strain::principal(current_strain)};
    }

    return {};
}
