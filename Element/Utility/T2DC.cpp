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

#include "T2DC.h"
#include <Element/Element.h>
#include <Toolbox/tensor.h>

unique_ptr<Orientation> T2DC::get_copy() { return make_unique<T2DC>(*this); }

void T2DC::update_transformation() {
    const auto coord = get_coordinate(element_ptr, 2);
    const auto t_disp = get_trial_displacement(element_ptr);

    vec x_axis(2);
    x_axis(0) = coord(1, 0) - coord(0, 0) + t_disp(2) - t_disp(0);
    x_axis(1) = coord(1, 1) - coord(0, 1) + t_disp(3) - t_disp(1);

    length = norm(x_axis);
    direction_cosine = x_axis / length;
    inclination = transform::atan2(direction_cosine);
}

T2DC::T2DC(const unsigned T)
    : T2DL(T) {}

bool T2DC::is_nlgeom() const { return true; }

mat T2DC::to_global_geometry_mat(const mat& l_mat) const {
    auto g_mat = to_global_stiffness_mat(-l_mat);

    g_mat.diag() += l_mat(0);
    for(unsigned I = 0; I < 2; ++I) g_mat(I + 2llu, I) = g_mat(I, I + 2llu) -= l_mat(0);

    return g_mat;
}
